"""
Oracle V8 EquationSet Abstraction
==================================

The EquationSet is the contract that defines which dynamical-core
equations V8 is integrating. It captures the assumptions about the
continuity constraint, the buoyancy formulation, and the elliptic
operator that the pressure projection inverts.

Concrete subclasses for V8:
    - LH82AnelasticEquationSet: Lipps & Hemler (1982) anelastic system.
      Continuity constraint is ∇·(ρ̄(z) u) = 0.
      Standard reference for nonhydrostatic deep-convection modeling
      under modest perturbation amplitudes. The default for V8.0.

    - PseudoIncompressibleEquationSet: Durran (1989) pseudo-incompressible
      system. Continuity constraint is ∇·(ρ̄(z) θ̄(z) u) = 0 (or
      equivalent pseudo-density form). Permits larger thermodynamic
      perturbations than LH82; planned for V8.x to characterize where
      LH82's small-perturbation assumption affects TC eyewall results.

The point of the abstraction: every component that depends on the
governing equations (advection, buoyancy forcing, the projection
operator, the constraint diagnostic) interrogates an EquationSet object
rather than hard-coding the equation form. Switching equation sets is a
single object swap; the rest of the solver is equation-set-agnostic.

References:
    Lipps, F. B., and R. S. Hemler (1982): A scale analysis of deep moist
        convection and some related numerical calculations. J. Atmos.
        Sci., 39, 2192-2210.
    Durran, D. R. (1989): Improving the anelastic approximation. J.
        Atmos. Sci., 46, 1453-1461.
    Klein, R., U. Achatz, D. Bresch, O. M. Knio, and P. K. Smolarkiewicz
        (2010): Regime of validity of soundproof atmospheric flow models.
        J. Atmos. Sci., 67, 3226-3237.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

# Standard gravity (m/s²).  Module-level rather than buried in the buoyancy
# inner loop.  NOTE: the same literal 9.81 is still duplicated across ~16 other
# modules (vortex_init, tendency, production_config, …); a project-wide
# constants module is the proper home — see the follow-up cleanup task.
GRAVITY = 9.81

if TYPE_CHECKING:
    # Forward-declare to avoid circular imports; these will exist
    # in reference_state/ and grid/ when the modules are written.
    from oracle_v8.reference_state.base_states import DryBaseState
    from oracle_v8.grid.staggering import GridStaggering


@dataclass
class ConstraintResidual:
    """The output of evaluating the continuity constraint on a velocity field.

    For LH82: ∇·(ρ̄ u). Should be ~0 to machine precision after projection.
    For PI:   ∇·(ρ̄ θ̄ u). Same expectation.

    Attributes:
        max_abs_residual: max |residual| over the domain
        l2_residual: discrete L2 norm of the residual
        residual_field: the full 3D residual (for inspection)
        constraint_name: human-readable name of the constraint, e.g.
            "div(rho_bar * u)" for LH82 or "div(rho_bar * theta_bar * u)" for PI
    """
    max_abs_residual: float
    l2_residual: float
    residual_field: np.ndarray
    constraint_name: str


class EquationSet(abc.ABC):
    """
    Abstract base class for the V8 dynamical-core governing equations.

    A concrete EquationSet captures three things:

        1. The form of the continuity constraint that the projection
           must enforce (and the diagnostic that checks it).

        2. The coefficient that appears in the variable-coefficient
           Poisson equation: ∇·(C(z) ∇φ) = f.
           For LH82: C(z) = ρ̄(z).
           For PI:   C(z) = ρ̄(z) θ̄(z) (or equivalent pseudo-density).

        3. The buoyancy formulation: how potential-temperature
           perturbations enter the vertical momentum equation.

    The abstraction is intentionally narrow. It does NOT encapsulate
    advection schemes, time integration, or surface forcing — those are
    TendencyComponent concerns. EquationSet captures only what changes
    when you swap LH82 for PI.
    """

    name: str = "abstract"

    @abc.abstractmethod
    def poisson_coefficient(
        self, base: "DryBaseState", staggering: "GridStaggering"
    ) -> np.ndarray:
        """
        Return the coefficient C(z) that multiplies ∇φ inside the
        divergence in the elliptic equation ∇·(C(z) ∇φ) = f.

        Returned shape and staggering must match what the Poisson solver
        expects. The staggering object handles where C(z) lives on the
        grid (typically half-levels for the divergence stencil).

        For LH82: returns ρ̄(z) on the half-levels where w lives.
        For PI:   returns ρ̄(z) θ̄(z) on the same half-levels.
        """

    @abc.abstractmethod
    def compute_constraint_residual(
        self,
        u: np.ndarray, v: np.ndarray, w: np.ndarray,
        base: "DryBaseState",
        staggering: "GridStaggering",
        dx: float, dy: float,
    ) -> ConstraintResidual:
        """
        Evaluate the continuity constraint ∇·(C(z) u) on the velocity
        field. Returns a ConstraintResidual.

        For LH82: ∇·(ρ̄ u). Should be ~0 after projection.
        For PI:   ∇·(ρ̄ θ̄ u). Same expectation.

        This is the diagnostic the validation suite calls every step to
        verify the projection is doing its job.
        """

    @abc.abstractmethod
    def compute_buoyancy_tendency(
        self,
        theta_prime: np.ndarray,
        base: "DryBaseState",
        staggering: "GridStaggering",
    ) -> np.ndarray:
        """
        Compute the buoyancy contribution to the vertical momentum
        tendency, given the potential-temperature perturbation θ′.

        For LH82 and PI, the leading-order form is:
            b = g · θ′ / θ̄(z)
        evaluated on the w-staggering levels (half-levels in Lorenz,
        or wherever w lives in the chosen staggering).

        Returns the buoyancy field on the w-staggering, in m/s².
        """

    @abc.abstractmethod
    def base_state_compatibility(self, base: "DryBaseState") -> tuple[bool, str]:
        """
        Verify that a given base state is compatible with this equation
        set's assumptions.

        For LH82: returns (False, reason) if θ̄(z) variation is so steep
        that small-perturbation assumption is questionable, or if base
        state is missing required fields.

        For PI: returns (False, reason) if base state lacks the
        thermodynamic closure needed for the pseudo-density formulation.

        The validation suite calls this before any test to give an
        early, clear error rather than a downstream numerical surprise.
        """

    def __repr__(self) -> str:
        return f"<EquationSet: {self.name}>"


class LH82AnelasticEquationSet(EquationSet):
    """
    Lipps & Hemler (1982) anelastic equation set.

    Continuity constraint: ∇·(ρ̄(z) u) = 0
    Poisson coefficient:   ρ̄(z)
    Buoyancy:              b = g · θ′ / θ̄(z)

    Assumptions:
        - Small thermodynamic perturbations: |θ′| << θ̄(z).
        - Hydrostatically balanced base state.
        - Sound waves filtered analytically by the divergence constraint.

    Limitations:
        - The small-perturbation assumption is questionable in TC
          eyewalls where θ′ can reach 10-20 K against θ̄ ≈ 300 K.
          This is a 3-7% perturbation; whether LH82's leading-order
          neglect of higher-order terms is acceptable at this magnitude
          is a question V8 is designed to investigate empirically.
    """

    name = "LH82_anelastic"

    def poisson_coefficient(
        self, base: "DryBaseState", staggering: "GridStaggering"
    ) -> np.ndarray:
        # LH82 uses ρ̄(z) as the divergence-form coefficient.
        # Implementation deferred until staggering API exists.
        raise NotImplementedError(
            "Will return ρ̄(z) on the staggering's half-levels. "
            "Implementation pending the GridStaggering interface."
        )

    def compute_constraint_residual(
        self,
        u: np.ndarray, v: np.ndarray, w: np.ndarray,
        base: "DryBaseState",
        staggering: "GridStaggering",
        dx: float, dy: float,
    ) -> ConstraintResidual:
        raise NotImplementedError(
            "Will compute ∇·(ρ̄u) using staggering's derivative operators. "
            "The existing compute_anelastic_residual in test_harness.py "
            "is a candidate implementation pending the staggering API."
        )

    def compute_buoyancy_tendency(
        self,
        theta_prime: np.ndarray,
        base: "DryBaseState",
        staggering: "GridStaggering",
    ) -> np.ndarray:
        """
        Compute the LH82 buoyancy tendency b = g · θ′ / θ̄(z) on the
        w-staggering levels.

        Algorithm:
            1. Evaluate b_full = g · θ′ / θ̄(z) on full levels, where
               θ′ and θ̄ both naturally live in the Lorenz convention.
            2. Map b_full to half levels via the staggering's
               symmetric 0.5-averaging interpolation. This produces a
               (..., nz+1) array with boundaries pre-zeroed by the
               generic interpolation.
            3. Explicitly assert the rigid-w kinematic constraint at
               surface and lid: b_half[..., 0] = b_half[..., -1] = 0.
               This is redundant given the staggering's default-zero
               boundary, but it documents the *physical* reason
               (rigid w forces dw/dt=0 at boundaries) at the equation-
               set level rather than relying on interpolation defaults.

        Returns
        -------
        b_half : np.ndarray, shape (..., nz+1)
            Buoyancy tendency on the w-staggering levels, in m/s².
            Suitable for direct assignment to a Tendency's dw_dt field.
        """
        # Step 1: physics on full levels
        # theta_prime has shape (..., nz); base.theta0 has shape (nz,)
        # numpy broadcasts the 1-D theta0 along the trailing axis of theta_prime
        b_full = GRAVITY * theta_prime / base.theta0

        # Step 2: symmetric interpolation to half levels
        b_half = staggering.interpolate_full_to_half(b_full)

        # Step 3: assert the rigid-w kinematic constraint at boundaries
        # (redundant with staggering's default-zero, but documents intent)
        b_half[..., 0] = 0.0
        b_half[..., -1] = 0.0

        return b_half

    def base_state_compatibility(self, base: "DryBaseState") -> tuple[bool, str]:
        """
        LH82 requires a hydrostatically balanced dry base state with
        finite, positive ρ̄(z) and θ̄(z) throughout, plus a discrete
        hydrostatic operator that V8 trusts.

        Hardening per Five's P1.5: check for None, missing required
        fields, non-finite values, and the load-bearing-readiness flag.
        """
        if base is None:
            return (False, "base state is None")

        # Required fields. Use hasattr rather than dataclass introspection
        # so we tolerate base states constructed from anywhere (test
        # harness, real ref state, future moist V8.x).
        required_fields = ("rho0", "theta0", "z", "integration_scheme")
        for fname in required_fields:
            if not hasattr(base, fname):
                return (False, f"base state is missing required field: {fname!r}")
            if getattr(base, fname) is None:
                return (False, f"base state field {fname!r} is None")

        # The is_load_bearing_ready() flag is the contract that the
        # base state was built using V8's discrete vertical operator,
        # not the trapezoidal placeholder. LH82 must reject placeholder
        # base states.
        ready_method = getattr(base, "is_load_bearing_ready", None)
        if not callable(ready_method) or not ready_method():
            return (False,
                    "base state is not load-bearing ready "
                    f"(integration_scheme={base.integration_scheme!r}); "
                    "construct using V8's discrete vertical operator before "
                    "treating LH82 results as load-bearing")

        # Finite-value checks. NaN or inf in the base state will silently
        # propagate through the solver and produce confusing failures
        # downstream; catch them here.
        if not np.all(np.isfinite(base.rho0)):
            return (False,
                    f"non-finite ρ̄(z) detected; "
                    f"min={np.nanmin(base.rho0)}, max={np.nanmax(base.rho0)}")
        if not np.all(np.isfinite(base.theta0)):
            return (False,
                    f"non-finite θ̄(z) detected; "
                    f"min={np.nanmin(base.theta0)}, max={np.nanmax(base.theta0)}")

        # Positivity checks (anelastic systems require positive ρ̄, θ̄).
        if np.any(base.rho0 <= 0):
            return (False, f"non-positive ρ̄(z); min = {base.rho0.min()}")
        if np.any(base.theta0 <= 0):
            return (False, f"non-positive θ̄(z); min = {base.theta0.min()}")

        return (True, "compatible")


class PseudoIncompressibleEquationSet(EquationSet):
    """
    Durran (1989) pseudo-incompressible equation set.

    The exact functional form of the constraint and the Poisson
    coefficient is TBD pending derivation from Durran (1989) and
    Klein et al. (2010). The continuity-constraint form ∇·(ρ̄θ̄ u) = 0
    appears commonly in the literature as the leading-order pseudo-
    density variant, but Durran's full PI system involves a slightly
    more involved closure on the divergence operator and the exact
    pseudo-density product depends on whether one chooses ρ̄θ̄,
    ρ̄(p̄/p_ref)^(c_v/c_p), or a related thermodynamic combination.

    PROVISIONAL form (pending the V8.x derivation pass):
        Constraint:           ∇·(<pseudo-density>(z) u) = 0
        Poisson coefficient:  <pseudo-density>(z)
        Buoyancy:             b = g · θ′ / θ̄(z) (same leading order as LH82)

    Advantages over LH82 (independent of the exact form):
        - No small-perturbation assumption on θ′.
        - Thermodynamically consistent under large diabatic heating
          (Klein et al. 2010, JAS).
        - Extends the validity regime of soundproof models by orders of
          magnitude in perturbation amplitude.

    Implementation note:
        The constraint involves a thermodynamic-product pseudo-density
        rather than ρ̄, so the Poisson operator's coefficient changes
        but the structure is identical. The hybrid FFT-tridiagonal
        solver works the same way; only the C(z) coefficient changes.
        This is exactly what the EquationSet abstraction is for.

    Status: planned for V8.x. V8.0 ships LH82 as the validated baseline.
    The exact PI form is a derivation task to perform when V8.0 is
    stable and we have a working comparison baseline.
    """

    name = "pseudo_incompressible"

    def poisson_coefficient(
        self, base: "DryBaseState", staggering: "GridStaggering"
    ) -> np.ndarray:
        raise NotImplementedError(
            "Will return the pseudo-density on the staggering's half-levels. "
            "Exact form pending Durran 1989 derivation; not implemented in V8.0."
        )

    def compute_constraint_residual(
        self,
        u: np.ndarray, v: np.ndarray, w: np.ndarray,
        base: "DryBaseState",
        staggering: "GridStaggering",
        dx: float, dy: float,
    ) -> ConstraintResidual:
        raise NotImplementedError(
            "Planned for V8.x; exact constraint form pending Durran derivation."
        )

    def compute_buoyancy_tendency(
        self,
        theta_prime: np.ndarray,
        base: "DryBaseState",
        staggering: "GridStaggering",
    ) -> np.ndarray:
        raise NotImplementedError("Planned for V8.x.")

    def base_state_compatibility(self, base: "DryBaseState") -> tuple[bool, str]:
        return (False, "PseudoIncompressibleEquationSet not implemented in V8.0; "
                "planned for V8.x with derivation pending")
