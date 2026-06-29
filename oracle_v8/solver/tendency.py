"""
Oracle V8 Step Components
==========================

The unit of physics in V8 is a StepComponent — anything that
participates in the step loop. Step components fall into three
categories with different interfaces:

    TendencyComponent — produces a Tendency that is summed into the
        provisional state update. Most physics lives here: buoyancy,
        pressure-gradient, advection, Coriolis, surface drag, sponge
        damping.

    ProjectionComponent — enforces a constraint by transforming the
        State directly, not by adding a tendency. The anelastic
        projection lives here: it solves the variable-coefficient
        Poisson equation and corrects velocities to satisfy the
        continuity constraint after the provisional step. Projection
        is NOT additive; it cannot be folded into the tendency sum.

    DiagnosticComponent — read-only consumer of state. Produces records
        (mass budgets, energy budgets, constraint residuals) without
        modifying state. Logged at every step.

Each StepComponent declares:
    - name (for diagnostics and logging)
    - stage (PRE_PROJECTION, PROJECTION, POST_PROJECTION, SLOW, DIAGNOSTIC)
    - reads() and writes() (state variables it touches; type-checked
      via the StateVar enum so typos cannot silently pass)

Design rationale: Five's P0.1 review correctly identified that treating
projection as a normal additive tendency component would force the
solver to special-case it. The StepComponent abstraction with stages
makes the special-case explicit and the solver loop straightforward.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

from oracle_v8.backend import xp as np
from oracle_v8.backend import to_numpy, asarray, xp
import os as _os
import numpy as _np_cpu      # always CPU numpy — used in Helmholtz FFT

if TYPE_CHECKING:
    from oracle_v8.solver.equation_set import EquationSet
    from oracle_v8.grid.staggering import GridStaggering
    from oracle_v8.reference_state.base_states import DryBaseState


# -------------------------------------------------------------------------
# State variable enumeration (Five's P1.8)
# -------------------------------------------------------------------------


class StateVar(str, Enum):
    """
    Canonical names for V8 prognostic state variables.

    Use this enum (not bare strings) when declaring component reads()
    and writes(). This catches the class of bug where a component's
    metadata says reads=("theta",) but the State field is named
    "theta_prime" — the enum makes such typos impossible at runtime.

    Adding a new prognostic variable means adding to this enum AND
    extending the State dataclass below.
    """
    U = "u"
    V = "v"
    W = "w"
    THETA_PRIME = "theta_prime"
    PROJECTION_POTENTIAL = "projection_potential"


# -------------------------------------------------------------------------
# Step stages (Five's P1.4)
# -------------------------------------------------------------------------


class StepStage(Enum):
    """
    Where in the step sequence a StepComponent fires.

    PRE_PROJECTION:
        Fast-mode tendency components that produce contributions to
        u*, the provisional velocity field, before projection enforces
        the continuity constraint. Examples: buoyancy, pressure
        gradient. In RK3 split-explicit time integration these are
        sub-stepped within each outer RK3 step.

    PROJECTION:
        The constraint-enforcement step. ProjectionComponent subclasses
        live here. Applied after PRE_PROJECTION tendencies have been
        accumulated into u*.

    POST_PROJECTION:
        Rare. Components that need to fire after projection but before
        the next sub-step. Reserved.

    SLOW:
        Components evaluated at the outer RK3 timestep. Examples:
        advection, Coriolis, surface drag, sponge damping. These do
        not need acoustic-mode-rate sub-stepping.

    DIAGNOSTIC:
        Read-only components that produce records, not state changes.
    """
    PRE_PROJECTION = "pre_projection"
    PROJECTION = "projection"
    POST_PROJECTION = "post_projection"
    SLOW = "slow"
    DIAGNOSTIC = "diagnostic"


# -------------------------------------------------------------------------
# Half-level convention (Five's P0.3)
# -------------------------------------------------------------------------
#
# V8 LOCKS THE FOLLOWING CONVENTION FOR LORENZ STAGGERING:
#
#   Full-level fields (u, v, theta_prime, projection_potential):
#       shape = (nx, ny, nz)
#       indexed by k = 0, 1, ..., nz-1 (cell centers)
#
#   Half-level field (w):
#       shape = (nx, ny, nz+1)
#       indexed by k = 0, 1, ..., nz (cell faces including boundaries)
#       w[:, :, 0]    = surface boundary value (= 0 for rigid surface)
#       w[:, :, nz]   = top boundary value (= 0 for rigid lid)
#       w[:, :, k]    for k = 1, ..., nz-1 = interior half-levels
#
# The boundary values are stored EXPLICITLY rather than implied by BC
# logic in operator code. This means the divergence operator can compute
# d(rho_bar w)/dz without special-casing the boundaries: the boundary
# w values are simply zero in the array. Documenting and enforcing this
# convention now (P0.3) prevents the off-by-one and shape-mismatch bugs
# that historically plague variable-staggering atmospheric models.
# -------------------------------------------------------------------------


@dataclass
class State:
    """
    The full prognostic state of the V8 system at a given time.

    Field-by-field shape conventions:

        u, v: shape (nx, ny, nz) — full levels, horizontal momentum
        w: shape (nx, ny, nz+1) — half levels INCLUDING boundaries.
            w[:, :, 0] is the surface (= 0 for rigid surface).
            w[:, :, nz] is the top (= 0 for rigid lid).
        theta_prime: shape (nx, ny, nz) — full levels
        projection_potential: shape (nx, ny, nz) — full levels.
            DIAGNOSTIC, not prognostic. Recomputed each step by the
            ProjectionComponent. Has units of pressure-like quantity
            but is a Lagrange multiplier enforcing ∇·(C(z) u) = 0,
            NOT a meteorological pressure perturbation.

    `t` is the current simulation time in seconds.

    Future moist V8.x will extend this dataclass with q_v, q_l, q_i.
    """
    u: np.ndarray
    v: np.ndarray
    w: np.ndarray
    theta_prime: np.ndarray
    projection_potential: np.ndarray
    t: float


@dataclass
class Tendency:
    """
    A contribution to the time derivative of state, returned by
    TendencyComponents.

    Note: there is intentionally NO d_projection_potential_dt field.
    The projection potential is diagnostic — it is recomputed by the
    ProjectionComponent each step rather than evolved by tendency
    accumulation. If a future equation set prognostically evolves
    pressure (e.g. a fully compressible V9 with Exner perturbation),
    that will be a separate prognostic variable, not a redefinition
    of projection_potential.

    A component returns a Tendency with non-zero entries only for the
    variables it modifies. The solver sums Tendency objects across
    components before stepping.
    """
    du_dt: np.ndarray | None = None
    dv_dt: np.ndarray | None = None
    dw_dt: np.ndarray | None = None
    dtheta_prime_dt: np.ndarray | None = None

    @classmethod
    def zeros_like(cls, state: State) -> "Tendency":
        """
        Return a Tendency with zero arrays matching the shapes of the
        given state. Useful as the initial value for tendency
        accumulation: total = Tendency.zeros_like(state).
        """
        return cls(
            du_dt=np.zeros_like(state.u),
            dv_dt=np.zeros_like(state.v),
            dw_dt=np.zeros_like(state.w),
            dtheta_prime_dt=np.zeros_like(state.theta_prime),
        )

    def add_(self, other: "Tendency") -> "Tendency":
        """
        In-place additive merge of another Tendency. Treats None entries
        as zero. Returns self for chaining.

        Used by the solver to accumulate tendencies from multiple
        components: total = Tendency.zeros_like(state); for c in
        components: total.add_(c.compute_tendency(...)).
        """
        for fname in ("du_dt", "dv_dt", "dw_dt", "dtheta_prime_dt"):
            other_val = getattr(other, fname)
            if other_val is None:
                continue
            current_val = getattr(self, fname)
            if current_val is None:
                setattr(self, fname, other_val.copy())
            else:
                current_val += other_val  # in-place
        return self

    def validate_against_state(self, state: State) -> None:
        """
        Verify all non-None tendency fields have shapes matching their
        corresponding state fields. Raises ValueError on mismatch.

        Use in debugging or test mode; the solver inner loop should
        call validate_against_state at most once per outer step.
        """
        pairs = [
            ("du_dt", "u"),
            ("dv_dt", "v"),
            ("dw_dt", "w"),
            ("dtheta_prime_dt", "theta_prime"),
        ]
        for tend_name, state_name in pairs:
            tend_val = getattr(self, tend_name)
            if tend_val is None:
                continue
            state_val = getattr(state, state_name)
            if tend_val.shape != state_val.shape:
                raise ValueError(
                    f"Tendency.{tend_name} shape {tend_val.shape} does not "
                    f"match State.{state_name} shape {state_val.shape}"
                )


# -------------------------------------------------------------------------
# Step component hierarchy (Five's P0.1 + P1.4)
# -------------------------------------------------------------------------


class StepComponent(abc.ABC):
    """
    Abstract base for anything that participates in the V8 step loop.

    Three concrete subclasses serve different semantics:
        TendencyComponent — produces a Tendency to be summed
        ProjectionComponent — transforms State directly (constraint
            enforcement)
        DiagnosticComponent — produces records, does not modify state

    All step components declare their stage, name, and read/write
    contracts. The solver dispatches by stage and by type.
    """

    name: str = "abstract"
    stage: StepStage = StepStage.SLOW

    @abc.abstractmethod
    def reads(self) -> tuple[StateVar, ...]:
        """State variables this component reads."""

    @abc.abstractmethod
    def writes(self) -> tuple[StateVar, ...]:
        """State variables this component writes (or modifies indirectly)."""

    def to_log_dict(self) -> dict:
        """Serializable description for validation logging."""
        return {
            "name": self.name,
            "stage": self.stage.value,
            "reads": [v.value for v in self.reads()],
            "writes": [v.value for v in self.writes()],
            "class": type(self).__name__,
            "kind": self._kind(),
        }

    def _kind(self) -> str:
        if isinstance(self, TendencyComponent):
            return "tendency"
        if isinstance(self, ProjectionComponent):
            return "projection"
        if isinstance(self, DiagnosticComponent):
            return "diagnostic"
        return "unknown"

    def __repr__(self) -> str:
        return f"<{type(self).__name__}: {self.name} ({self.stage.value})>"


class TendencyComponent(StepComponent):
    """
    StepComponents that produce a Tendency. Most physics lives here.
    """

    @abc.abstractmethod
    def compute_tendency(
        self,
        state: State,
        equation_set: "EquationSet",
        staggering: "GridStaggering",
        base: "DryBaseState",
        dt: float,
    ) -> Tendency:
        """
        Compute the tendency this component contributes given the
        current state. Pure function: no hidden component state, fully
        determined by inputs.
        """


class ProjectionComponent(StepComponent):
    """
    StepComponents that enforce constraints by transforming State
    directly. CANNOT be folded into tendency accumulation — Five's
    P0.1 finding.

    The canonical example is the anelastic projection: solve the
    variable-coefficient Poisson equation, correct velocities to
    satisfy the continuity constraint, return the corrected State.
    """

    stage: StepStage = StepStage.PROJECTION

    @abc.abstractmethod
    def apply_projection(
        self,
        state: State,
        equation_set: "EquationSet",
        staggering: "GridStaggering",
        base: "DryBaseState",
        dt: float,
    ) -> State:
        """
        Transform the state to satisfy the constraint. Returns a new
        State (or modifies in-place and returns self; concrete
        subclasses document which).
        """


class DiagnosticComponent(StepComponent):
    """
    Read-only StepComponents that produce records without modifying
    state. Mass budget, energy budget, constraint residual, vorticity
    statistics — all live here.
    """

    stage: StepStage = StepStage.DIAGNOSTIC

    @abc.abstractmethod
    def compute_diagnostic(
        self,
        state: State,
        equation_set: "EquationSet",
        staggering: "GridStaggering",
        base: "DryBaseState",
    ) -> dict:
        """
        Compute and return a serializable diagnostic record.
        """

    def writes(self) -> tuple[StateVar, ...]:
        # Diagnostics never write state, so the default is empty.
        return ()


# -------------------------------------------------------------------------
# Shared spatial helpers
# -------------------------------------------------------------------------


def anelastic_divergence(u, v, w, rho_bar_full, rho_bar_half, kx, ky, dz):
    """
    Discrete d = ∇·(ρ̄ u) on the cell-centered Lorenz grid — the EXACT operator
    the AnelasticProjection inverts and corrects against.

    Single source of truth: both AnelasticProjection (via
    _compute_anelastic_divergence) and the validation harness
    (validation/test_harness.compute_anelastic_residual) call this, so the
    constraint a test measures is the one the solver enforces — never a
    re-derived stencil.

    Discretization:
      - Horizontal ∂(ρ̄u)/∂x + ∂(ρ̄v)/∂y: spectral (multiply by ikx, iky),
        periodic in x, y.
      - Vertical ∂(ρ̄w)/∂z: flux form across the Δz half-level stencil,
        (ρ̄w)_{k+1/2} − (ρ̄w)_{k−1/2}, with w on half levels (length nz+1)
        and rigid w=0 stored explicitly at the surface/lid faces.

    Parameters
    ----------
    u, v : (nx, ny, nz) full-level horizontal velocity.
    w : (nx, ny, nz+1) half-level vertical velocity.
    rho_bar_full : (nz,) ρ̄ on full levels.
    rho_bar_half : (nz+1,) ρ̄ on half levels.
    kx, ky : wavenumber arrays (2π·fftfreq) matching the grid.
    dz : vertical spacing.

    Returns d on full levels, shape (nx, ny, nz).
    """
    rho_u = rho_bar_full[None, None, :] * u
    rho_v = rho_bar_full[None, None, :] * v

    ikx = 1j * kx[:, None, None]
    iky = 1j * ky[None, :, None]

    rho_u_hat = np.fft.fft2(rho_u, axes=(0, 1))
    rho_v_hat = np.fft.fft2(rho_v, axes=(0, 1))
    d_hat = ikx * rho_u_hat + iky * rho_v_hat
    horiz_div = np.real(np.fft.ifft2(d_hat, axes=(0, 1)))

    rho_w_lower = rho_bar_half[None, None, :-1] * w[:, :, :-1]
    rho_w_upper = rho_bar_half[None, None, 1:] * w[:, :, 1:]
    vert_div = (rho_w_upper - rho_w_lower) / dz

    return horiz_div + vert_div


def _phi_gradient(
    phi: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    dz: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Spectral-horizontal + finite-difference-vertical gradient of φ.

    Shared by AnelasticProjection and PressureGradientComponent so
    the stencil is defined exactly once.  Both components use
    identical grid geometry, so their kx/ky arrays are numerically
    equal; the results are therefore bit-for-bit identical for the
    same φ.

    Horizontal
    ----------
    FFT differentiation on the periodic x–y planes.  Output lives on
    full levels (same grid as φ), shape (nx, ny, nz).

    Vertical
    --------
    Centered finite differences, full levels → half levels:
        ∂φ/∂z|_{k+½} = (φ_{k+1} − φ_k) / dz    k = 0 … nz−2
    Output shape (nx, ny, nz+1).  Surface (index 0) and lid
    (index nz) entries are zero by Neumann-BC convention — the rigid
    boundary forces ∂φ/∂z = 0 at both faces, so the velocity
    correction there is zero and the rigid-wall w = 0 is preserved.

    Parameters
    ----------
    phi : ndarray, shape (nx, ny, nz)
        Projection potential on full levels.
    kx  : ndarray, shape (nx,)
        Zonal wavenumber array (rad/m).  Typically
        ``2π · np.fft.fftfreq(nx, d=Lx/nx)``.
    ky  : ndarray, shape (ny,)
        Meridional wavenumber array (rad/m).
    dz  : float
        Uniform vertical grid spacing (m).

    Returns
    -------
    dphi_dx : ndarray, shape (nx, ny, nz) — ∂φ/∂x on full levels
    dphi_dy : ndarray, shape (nx, ny, nz) — ∂φ/∂y on full levels
    dphi_dz : ndarray, shape (nx, ny, nz+1) — ∂φ/∂z on half levels
    """
    phi_hat = np.fft.fft2(phi, axes=(0, 1))
    ikx = 1j * kx[:, None, None]
    iky = 1j * ky[None, :, None]
    dphi_dx = np.real(np.fft.ifft2(ikx * phi_hat, axes=(0, 1)))
    dphi_dy = np.real(np.fft.ifft2(iky * phi_hat, axes=(0, 1)))

    nz = phi.shape[2]
    nx, ny = phi.shape[0], phi.shape[1]
    dphi_dz = np.zeros((nx, ny, nz + 1), dtype=phi.dtype)
    # Interior half-levels: centered difference between adjacent full levels
    dphi_dz[:, :, 1:-1] = (phi[:, :, 1:] - phi[:, :, :-1]) / dz
    # Boundary half-levels stay at zero (Neumann BC, rigid surface and lid)

    return dphi_dx, dphi_dy, dphi_dz


# -------------------------------------------------------------------------
# Concrete component declarations
# -------------------------------------------------------------------------


class BuoyancyComponent(TendencyComponent):
    """Buoyancy forcing: dw/dt += g · θ′ / θ̄(z)."""
    name = "buoyancy"
    stage = StepStage.PRE_PROJECTION

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        """
        Buoyancy tendency on the w-staggering, delegated to the equation
        set per the EquationSet contract.

        This component is a pure orchestrator. It:
          1. Asks the equation set for the buoyancy tendency on
             half-levels (the equation set internally computes physics
             on full levels and uses the staggering's interpolation).
          2. Validates that the returned shape matches state.w
             (i.e., the equation set honored its contract).
          3. Wraps the result in a sparse Tendency with only dw_dt
             populated.

        The physics formulation (LH82: g·θ′/θ̄; PI: same to leading
        order with higher-order corrections; etc.) lives in the
        equation set. The spatial discretization lives in the
        staggering. The component itself contains no physics and no
        spatial-mapping logic — it just composes those abstractions.
        """
        b_half = equation_set.compute_buoyancy_tendency(
            state.theta_prime, base, staggering,
        )
        # Validate the equation set honored the w-staggering contract.
        # This catches subtle bugs like an equation set returning b on
        # full levels instead of half levels.
        if b_half.shape != state.w.shape:
            raise ValueError(
                f"equation_set.compute_buoyancy_tendency returned shape "
                f"{b_half.shape}, expected state.w.shape = {state.w.shape}. "
                f"The equation set must return buoyancy on w-staggering "
                f"levels per the EquationSet contract."
            )
        return Tendency(dw_dt=b_half)

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.THETA_PRIME,)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.W,)


class PressureGradientComponent(TendencyComponent):
    """
    Pressure-gradient forcing on momentum:
        du/dt += −∂φ/∂x
        dv/dt += −∂φ/∂y
        dw/dt += −∂φ/∂z
    where φ = state.projection_potential (full levels, shape nx×ny×nz).

    Grid parameters must be supplied at construction to configure the
    wavenumber arrays used for the spectral horizontal derivative:

        pg = PressureGradientComponent(nx=64, ny=64, nz=32,
                                       Lx=1e5, Ly=1e5, Lz=1e4)

    These must match the grid parameters given to AnelasticProjection.
    The same kx/ky convention guarantees that ∇φ here and −∇φ inside
    the projection's velocity correction use numerically identical
    stencils — the two operators are adjoints of each other and bit-
    for-bit consistency matters for energy conservation diagnostics.

    Construction without arguments (PressureGradientComponent()) is
    valid for smoke tests that inspect only metadata (name, stage,
    reads, writes).  Calling compute_tendency on an unconfigured
    instance raises RuntimeError with an actionable message.

    Architecture note
    -----------------
    The pressure gradient −∇φ is equation-set-agnostic at LH82 and
    PI leading order: both use the projection potential in the same
    form.  The spatial gradient therefore lives here rather than in
    the EquationSet.  If PI higher-order pressure corrections are
    needed (V8.x), add equation_set.pressure_gradient_correction()
    as a residual term accumulated after this component's output.

    The gradient stencil is implemented in the module-level helper
    _phi_gradient, shared with AnelasticProjection._compute_phi_gradient
    so the discretization is defined exactly once.
    """
    name = "pressure_gradient"
    stage = StepStage.PRE_PROJECTION

    def __init__(
        self,
        nx: int | None = None,
        ny: int | None = None,
        nz: int | None = None,
        Lx: float | None = None,
        Ly: float | None = None,
        Lz: float | None = None,
    ):
        self._kx: np.ndarray | None = None
        self._ky: np.ndarray | None = None
        self._dz: float | None = None
        if nx is not None and Lx is not None:
            self._configure(nx, ny, nz, Lx, Ly, Lz)

    def _configure(
        self,
        nx: int,
        ny: int,
        nz: int,
        Lx: float,
        Ly: float,
        Lz: float,
    ) -> None:
        """Build and cache the wavenumber arrays and dz from grid params."""
        self._dz = Lz / nz
        self._kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=Lx / nx)
        self._ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=Ly / ny)

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        """
        Compute −∇φ from state.projection_potential and return as a
        Tendency with all three momentum components populated.

        The component is a pure orchestrator:
          1. Validates the instance is configured.
          2. Delegates gradient computation to _phi_gradient().
          3. Validates output shapes against state.
          4. Returns Tendency(du_dt=−∂φ/∂x, dv_dt=−∂φ/∂y, dw_dt=−∂φ/∂z).

        No physics or spatial-mapping logic lives here beyond the
        delegation and shape contracts.
        """
        if self._kx is None:
            raise RuntimeError(
                "PressureGradientComponent.compute_tendency called on an "
                "unconfigured instance.  Supply grid parameters at construction:\n"
                "    PressureGradientComponent("
                "nx=..., ny=..., nz=..., Lx=..., Ly=..., Lz=...)"
            )

        phi = state.projection_potential
        dphi_dx, dphi_dy, dphi_dz = _phi_gradient(
            phi, self._kx, self._ky, self._dz,
        )

        # Shape contracts: identical to what BuoyancyComponent checks for dw_dt.
        if dphi_dx.shape != state.u.shape:
            raise ValueError(
                f"∂φ/∂x shape {dphi_dx.shape} != state.u.shape {state.u.shape}; "
                f"check that grid parameters match the simulation grid."
            )
        if dphi_dz.shape != state.w.shape:
            raise ValueError(
                f"∂φ/∂z shape {dphi_dz.shape} != state.w.shape {state.w.shape}; "
                f"PressureGradientComponent nz must match AnelasticProjection nz."
            )

        return Tendency(
            du_dt=-dphi_dx,
            dv_dt=-dphi_dy,
            dw_dt=-dphi_dz,
        )

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.PROJECTION_POTENTIAL,)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V, StateVar.W)


class AnelasticProjection(ProjectionComponent):
    """
    Anelastic constraint enforcement.

    Given a provisional state with velocities (u*, v*, w*) that may
    not satisfy ∇·(ρ̄ u) = 0, solve the variable-coefficient Poisson
    equation ∇·(ρ̄ ∇φ) = ∇·(ρ̄ u*), then correct:
        u_final = u* - ∂φ/∂x
        v_final = v* - ∂φ/∂y
        w_final = w* - ∂φ/∂z
    and store φ as state.projection_potential.

    Implementation uses the hybrid FFT-x,y + tridiagonal-z solver.
    Diagnostics (compatibility residual, discrete operator residual)
    are exposed via the last_solve attribute for logging.

    Construction caches the solver to avoid recomputing the wavenumber
    grids each step. If the grid resolution changes mid-run, a new
    AnelasticProjection instance must be constructed.
    """
    name = "anelastic_projection"
    # stage inherited from ProjectionComponent (= PROJECTION)

    def __init__(self, nx: int, ny: int, nz: int,
                 Lx: float, Ly: float, Lz: float):
        # Lazy import to avoid circular dependency: poisson.py imports
        # nothing from tendency, but we want the dependency direction
        # to be tendency → poisson, not the reverse.
        from oracle_v8.solver.poisson import VariableCoefficientPoissonSolver
        self._solver = VariableCoefficientPoissonSolver(
            nx, ny, nz, Lx, Ly, Lz,
        )
        self._dx = Lx / nx
        self._dy = Ly / ny
        self._dz = Lz / nz
        self.last_solve = None  # latest PoissonSolveResult for diagnostics
        self._rho_cache   = None  # (rho_full_ref, rho_half) — stable ids for
                                  # the solver's Thomas cache (V8.6.3)
        self._solve_count = 0

    def apply_projection(self, state, equation_set, staggering, base, dt):
        """
        Project the provisional velocity field onto the divergence-free
        manifold. Modifies state in-place: corrects u, v, w and updates
        projection_potential.

        Parameters
        ----------
        state : State
            Provisional state with possibly non-divergence-free velocity.
        equation_set, staggering, base : context objects
            Currently used for ρ̄(z) lookup; future moist V8.x may use
            equation_set for ρ̄θ̄ pseudo-density.
        dt : float
            Timestep, currently unused (the projection itself is not
            timestep-dependent in the LH82 formulation).

        Returns
        -------
        State (same instance, modified in-place).
        """
        # Pull ρ̄ on full and half levels from the base state.
        # base.rho0 is on full levels; we construct half-level values
        # by averaging adjacent full levels, with surface/top ghost
        # values set to the boundary cell value (cell-centered Neumann
        # convention).
        nz = self._solver.nz
        rho_bar_full = base.rho0
        # V8.6.3: build the half-level profile ONCE per base state.  Beyond
        # the (small) saving, the stable array identities let the solver's
        # Thomas precompute cache hit on every call.
        if self._rho_cache is None or self._rho_cache[0] is not rho_bar_full:
            _rbh = np.zeros(nz + 1)
            _rbh[1:-1] = 0.5 * (rho_bar_full[:-1] + rho_bar_full[1:])
            _rbh[0] = rho_bar_full[0]
            _rbh[-1] = rho_bar_full[-1]
            self._rho_cache = (rho_bar_full, _rbh)
        rho_bar_half = self._rho_cache[1]

        # Step 1: compute the source d = ∇·(ρ̄ u*)
        # u, v on full levels (nx, ny, nz); w on half levels (nx, ny, nz+1)
        # with explicit boundary values.
        d = self._compute_anelastic_divergence(
            state.u, state.v, state.w, rho_bar_full, rho_bar_half,
        )

        # Step 2: solve ∇·(ρ̄ ∇φ) = d
        # V8.6.3: the residual diagnostics cost a FULL second operator
        # application (two extra FFTs) plus two host syncs per solve.  Sample
        # them every ORACLE_POISSON_DIAG_EVERY solves (default 60) instead of
        # every solve; set the env var to 1 to restore historical behavior.
        # Skipped solves report residuals of 0.0 in last_solve.
        self._solve_count += 1
        _de = max(1, int(_os.environ.get("ORACLE_POISSON_DIAG_EVERY", "60")))
        result = self._solver.solve(
            d, rho_bar_full, rho_bar_half,
            log_diagnostics=(_de == 1 or self._solve_count % _de == 1),
        )
        self.last_solve = result

        # Step 3: correct velocities by subtracting ∇φ
        u_correction, v_correction, w_correction = self._compute_phi_gradient(
            result.phi,
        )
        state.u = state.u - u_correction
        state.v = state.v - v_correction
        # w correction needs to fit the (nz+1) half-level shape; the
        # gradient routine returns it sized correctly.
        state.w = state.w - w_correction

        # Step 4: store φ as the new projection_potential
        state.projection_potential = result.phi

        return state

    def _compute_anelastic_divergence(self, u, v, w, rho_bar_full, rho_bar_half):
        """
        Compute d = ∇·(ρ̄ u) on the cell-centered Lorenz grid.

        Thin wrapper around the module-level :func:`anelastic_divergence` so
        the projection and the validation harness share ONE operator (the
        test suite must never re-derive this stencil).  Returns d on full
        levels, shape (nx, ny, nz).
        """
        return anelastic_divergence(
            u, v, w, rho_bar_full, rho_bar_half,
            self._solver.kx, self._solver.ky, self._dz,
        )

    def _compute_phi_gradient(self, phi):
        """
        Compute ∇φ on the velocity grid.

        Delegates to the module-level _phi_gradient helper, which is also
        used by PressureGradientComponent — guaranteeing that the stencil
        applied during projection correction and the one applied during the
        pressure-gradient tendency step are numerically identical.

        Returns (du_correction, dv_correction, dw_correction) where
        du, dv have full-level shape (nx, ny, nz) and dw has half-level
        shape (nx, ny, nz+1).  Surface and lid entries of dw are zero,
        preserving the rigid-boundary condition on w.
        """
        return _phi_gradient(
            phi,
            self._solver.kx,
            self._solver.ky,
            self._dz,
        )

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V, StateVar.W)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V, StateVar.W,
                StateVar.PROJECTION_POTENTIAL)


class AdvectionComponent(TendencyComponent):
    """
    Full advection of all prognostic variables (V8.1).

    Advects θ′ AND momentum (u, v, w):

        ∂θ′/∂t = −u·∇θ′ − w_c·∂θ̄/∂z
        ∂u/∂t  = −u·∇u  (momentum advection)
        ∂v/∂t  = −u·∇v
        ∂w/∂t  = −u_h·∇w  (on half levels)

    The base-state term −w_c·∂θ̄/∂z in the θ′ equation captures the
    adiabatic adjustment as the parcel is displaced into an environment
    with different θ̄(z).  It is essential for physical secondary
    circulation spin-up.

    Spatial scheme: second-order centred differences throughout.
        Horizontal (x, y): periodic BCs via np.roll.
        Vertical (z): centred interior, one-sided at surface and lid.

    Staggering:
        u, v, θ′ live on full levels (nx, ny, nz).
        w lives on half levels (nx, ny, nz+1).
        For u, v advection: w is interpolated full←half (w_c).
        For w advection: u, v are interpolated half←full (u_h, v_h).

    Because this component now writes momentum fields, the integrator's
    _apply_slow will call apply_projection after both Strang half-steps
    to restore ∇·(ρ̄u) = 0.

    Grid parameters must be supplied at construction:
        adv = AdvectionComponent(nx=64, ny=64, nz=32,
                                  Lx=1e5, Ly=1e5, Lz=1e4)
    Construction without arguments is valid for smoke tests.
    """
    name = "advection"
    stage = StepStage.SLOW

    def __init__(
        self,
        nx: int | None = None,
        ny: int | None = None,
        nz: int | None = None,
        Lx: float | None = None,
        Ly: float | None = None,
        Lz: float | None = None,
        scheme: str = "centered2",
    ) -> None:
        # scheme selects the HORIZONTAL advection discretization:
        #   "centered2"  — 2nd-order centred (V8.0 default; bit-identical to the
        #                  validated runs; no built-in dissipation).
        #   "upwind5h"   — 5th-order upwind-biased (Wicker-Skamarock 2002) in x,y;
        #                  carries 4th-order numerical dissipation that damps the
        #                  nonlinear 2Δx cascade.  Vertical advection stays centred
        #                  in both modes (w is O(0.1 m/s); the aliasing instability
        #                  is horizontal).  Advective form — flux form for strict
        #                  conservation is a follow-up.
        if scheme not in ("centered2", "upwind5h"):
            raise ValueError(
                f"AdvectionComponent scheme must be 'centered2' or 'upwind5h', "
                f"got {scheme!r}"
            )
        self.scheme = scheme
        self._dx: float | None = None
        self._dy: float | None = None
        self._dz: float | None = None
        self._nz: int | None = nz
        if nx is not None and Lx is not None:
            self._dx = Lx / nx
            self._dy = Ly / ny
            self._dz = Lz / nz
            self._nz = nz

    # ------------------------------------------------------------------
    # Private stencil helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _cx(f: np.ndarray, dx: float) -> np.ndarray:
        """∂f/∂x, second-order centred, periodic in x."""
        return (np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0)) / (2.0 * dx)

    @staticmethod
    def _cy(f: np.ndarray, dy: float) -> np.ndarray:
        """∂f/∂y, second-order centred, periodic in y."""
        return (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2.0 * dy)

    @staticmethod
    def _cz_full(f: np.ndarray, dz: float) -> np.ndarray:
        """∂f/∂z for f on full levels; one-sided at surface and lid."""
        g = np.empty_like(f)
        g[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2.0 * dz)
        g[:, :,  0]   = (f[:, :,  1] - f[:, :,  0])  / dz
        g[:, :, -1]   = (f[:, :, -1] - f[:, :, -2])  / dz
        return g

    @staticmethod
    def _cz_half(f: np.ndarray, dz: float) -> np.ndarray:
        """∂f/∂z for f on half levels (nz+1); one-sided at boundaries."""
        g = np.empty_like(f)
        g[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2.0 * dz)
        g[:, :,  0]   = (f[:, :,  1] - f[:, :,  0])  / dz
        g[:, :, -1]   = (f[:, :, -1] - f[:, :, -2])  / dz
        return g

    @staticmethod
    def _upwind5(q: np.ndarray, a: np.ndarray, d: float, axis: int) -> np.ndarray:
        """
        5th-order upwind-biased first derivative ∂q/∂x_axis (Wicker-Skamarock
        2002), periodic via np.roll, upwind direction set by the advecting
        velocity `a` (same shape as q).

        Decomposition: 6th-order centred derivative minus |a|-signed 6th-difference
        dissipation, so the truncation error is a 5th-order (∝ ∂⁶q) dissipative
        term — exactly the built-in scale-selective damping that lets the scheme
        control the nonlinear 2Δx cascade without an external filter.

            a >= 0:  (-2 q_{i-3} +15 q_{i-2} -60 q_{i-1} +20 q_i +30 q_{i+1} -3 q_{i+2})/(60 d)
            a <  0:  ( 2 q_{i+3} -15 q_{i+2} +60 q_{i+1} -20 q_i -30 q_{i-1} +3 q_{i-2})/(60 d)
        """
        qm3 = np.roll(q,  3, axis); qm2 = np.roll(q,  2, axis); qm1 = np.roll(q,  1, axis)
        qp1 = np.roll(q, -1, axis); qp2 = np.roll(q, -2, axis); qp3 = np.roll(q, -3, axis)
        d_pos = (-2.0*qm3 + 15.0*qm2 - 60.0*qm1 + 20.0*q + 30.0*qp1 - 3.0*qp2) / (60.0 * d)
        d_neg = ( 2.0*qp3 - 15.0*qp2 + 60.0*qp1 - 20.0*q - 30.0*qm1 + 3.0*qm2) / (60.0 * d)
        return np.where(a >= 0.0, d_pos, d_neg)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        if self._dx is None:
            raise RuntimeError(
                "AdvectionComponent.compute_tendency called on an "
                "unconfigured instance.  Provide grid parameters:\n"
                "    AdvectionComponent("
                "nx=..., ny=..., nz=..., Lx=..., Ly=..., Lz=...)"
            )

        u  = state.u            # (nx, ny, nz)
        v  = state.v            # (nx, ny, nz)
        w  = state.w            # (nx, ny, nz+1)
        θp = state.theta_prime  # (nx, ny, nz)

        dx, dy, dz = self._dx, self._dy, self._dz

        # ------------------------------------------------------------------
        # Level interpolations
        # ------------------------------------------------------------------

        # w half → full (for u, v, θ′ advection)
        w_c = 0.5 * (w[:, :, :-1] + w[:, :, 1:])        # (nx, ny, nz)

        # u, v full → half (for w advection); extrapolate at boundaries
        u_h = np.empty_like(w)
        u_h[:, :, 1:-1] = 0.5 * (u[:, :, :-1] + u[:, :, 1:])
        u_h[:, :,  0]   = u[:, :,  0]
        u_h[:, :, -1]   = u[:, :, -1]

        v_h = np.empty_like(w)
        v_h[:, :, 1:-1] = 0.5 * (v[:, :, :-1] + v[:, :, 1:])
        v_h[:, :,  0]   = v[:, :,  0]
        v_h[:, :, -1]   = v[:, :, -1]

        # Horizontal advective derivative operators.  Centered ignores the
        # advecting velocity `a`; upwind5h uses its sign to bias the stencil.
        if self.scheme == "upwind5h":
            ddx = lambda q, a: self._upwind5(q, a, dx, 0)
            ddy = lambda q, a: self._upwind5(q, a, dy, 1)
        else:
            ddx = lambda q, a: self._cx(q, dx)
            ddy = lambda q, a: self._cy(q, dy)

        # ------------------------------------------------------------------
        # θ′ advection (full levels)
        #   ∂θ′/∂t = −u·∇θ′ − w_c·∂θ̄/∂z
        # ------------------------------------------------------------------
        θ0 = base.theta0                    # (nz,)
        dθ0_dz        = np.empty(self._nz)
        dθ0_dz[1:-1]  = (θ0[2:] - θ0[:-2]) / (2.0 * dz)
        dθ0_dz[0]     = (θ0[1]  - θ0[0])   / dz
        dθ0_dz[-1]    = (θ0[-1] - θ0[-2])  / dz

        dθp_dt = -(
            u   * ddx(θp, u)
            + v   * ddy(θp, v)
            + w_c * (self._cz_full(θp, dz) + dθ0_dz[None, None, :])
        )

        # ------------------------------------------------------------------
        # Momentum advection — u (full levels)
        #   ∂u/∂t = −(u ∂u/∂x + v ∂u/∂y + w_c ∂u/∂z)
        # ------------------------------------------------------------------
        du_dt = -(
            u   * ddx(u, u)
            + v   * ddy(u, v)
            + w_c * self._cz_full(u, dz)
        )

        # ------------------------------------------------------------------
        # Momentum advection — v (full levels)
        # ------------------------------------------------------------------
        dv_dt = -(
            u   * ddx(v, u)
            + v   * ddy(v, v)
            + w_c * self._cz_full(v, dz)
        )

        # ------------------------------------------------------------------
        # Momentum advection — w (half levels)
        #   ∂w/∂t = −(u_h ∂w/∂x + v_h ∂w/∂y + w ∂w/∂z)
        # Rigid-BC: w = 0 at surface (k=0) and lid (k=nz), so
        # dw/dt = 0 there too.
        # ------------------------------------------------------------------
        dw_dt = -(
            u_h * ddx(w, u_h)
            + v_h * ddy(w, v_h)
            + w   * self._cz_half(w, dz)
        )
        dw_dt[:, :,  0] = 0.0    # surface rigid BC
        dw_dt[:, :, -1] = 0.0    # lid rigid BC

        return Tendency(
            du_dt           = du_dt,
            dv_dt           = dv_dt,
            dw_dt           = dw_dt,
            dtheta_prime_dt = dθp_dt,
        )

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V, StateVar.W, StateVar.THETA_PRIME)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V, StateVar.W, StateVar.THETA_PRIME)


class CoriolisComponent(TendencyComponent):
    """
    Coriolis forcing — f-plane and β-plane modes.

    f-plane (default):
        du/dt += +f₀ · v
        dv/dt += −f₀ · u

        Uniform Coriolis parameter everywhere in the domain.  Correct for
        warm-bubble tests, component validation, and short integrations
        where meridional displacement is small compared to the Rossby
        deformation radius.

    β-plane:
        f(y) = f₀ + β · (y − y₀)

        du/dt += +f(y) · v
        dv/dt += −f(y) · u

        Spatially varying Coriolis parameter.  Required for multi-day TC
        track simulations.  Beta drift — the poleward-and-westward
        self-propagation of a vortex on a β-plane — is the physical
        mechanism that distinguishes real TC tracks from pure steering.
        V7 had a β-plane module; V8 inherits it here.

    Parameters
    ----------
    f : float
        Coriolis parameter f₀ (s⁻¹).  Default 5×10⁻⁵ s⁻¹ ≈ 20°N.
        Also serves as f₀ for the β-plane.
    mode : str
        "f_plane" (default) or "beta_plane".
    beta : float, optional
        Meridional gradient of f, df/dy (m⁻¹ s⁻¹).  If None and
        mode="beta_plane", derived from f₀ via:
            β = 2Ω cos(lat) / R_Earth
        where lat = arcsin(f₀ / 2Ω).  Approximately 2.17×10⁻¹¹ m⁻¹s⁻¹
        at 20°N.
    Ly : float, optional
        Domain meridional extent (m).  Required for beta_plane.
    ny : int, optional
        Number of meridional grid cells.  Required for beta_plane.
    y0 : float, optional
        Reference y-coordinate where f = f₀ (m).  Defaults to domain
        centre (Ly / 2).

    Sign convention (NH, f > 0):
        Eastward inflow  (u < 0): dv/dt = −f·u > 0 → cyclonic ✓
        Northward inflow (v < 0): du/dt = +f·v < 0 → cyclonic ✓
    """
    name = "coriolis"
    stage = StepStage.SLOW

    # Earth constants
    _OMEGA     = 7.2921e-5   # rad s⁻¹
    _R_EARTH   = 6.371e6     # m

    def __init__(
        self,
        f: float = 5e-5,
        mode: str = "f_plane",
        beta: float = None,
        Ly: float = None,
        ny: int = None,
        y0: float = None,
        u_env: float = 0.0,
        v_env: float = 0.0,
        periodic_taper: bool = False,
    ) -> None:
        import numpy as _np

        self._f0   = f
        self._mode = mode
        # Geostrophic background winds: Coriolis acts only on the PERTURBATION
        # velocity (u - u_env, v - v_env).  Without this, the uniform steering
        # flow has no geostrophic balance in the periodic domain and the entire
        # airmass undergoes an inertial oscillation with period T = 2pi/f ≈ 26h,
        # causing the TC to loop rather than translate (diagnosed by Gemini,
        # ensemble review May 2026).
        self._u_env = float(u_env)
        self._v_env = float(v_env)

        if mode == "beta_plane":
            if Ly is None or ny is None:
                raise ValueError(
                    "CoriolisComponent beta_plane mode requires Ly and ny."
                )

            # Derive β from f₀ if not provided
            if beta is None:
                lat_rad = _np.arcsin(
                    _np.clip(f / (2.0 * self._OMEGA), -1.0, 1.0)
                )
                beta = 2.0 * self._OMEGA * _np.cos(lat_rad) / self._R_EARTH

            self._beta = float(beta)

            # Cell-centred y coordinates (shape (ny,))
            dy       = Ly / ny
            y_ref    = y0 if y0 is not None else Ly / 2.0
            y_ctrs   = (_np.arange(ny) + 0.5) * dy

            from oracle_v8.backend import asarray
            f_deviation = self._beta * (y_ctrs - y_ref)  # β·(y - y₀) raw

            if periodic_taper:
                # CRITICAL FIX: f(y) = f₀ + β·y is LINEAR and DISCONTINUOUS
                # at the N/S periodic boundary: Δf = β·Ly = 4.34×10⁻⁵ s⁻¹
                # (64% of f₀).  This delta-function of divergence floods the
                # FFT Poisson solver with Gibbs noise at every step, driving
                # the φ growth observed in Run 10.  (Diagnosed Gemini, second
                # ensemble review May 2026.)
                #
                # Fix: apply a cosine² taper to the β·(y-y₀) deviation in
                # the outer taper_frac of the domain at each N and S boundary,
                # so f(0) = f(Ly) = f₀ — exactly continuous for the FFT.
                # The interior (1 - 2*taper_frac) fraction is the true
                # beta-plane; the taper zones damp to f₀ smoothly.
                taper_frac  = 0.20          # 20% at each boundary = 400 km
                n_taper     = max(1, int(round(taper_frac * ny)))
                taper       = _np.ones(ny)
                for i in range(n_taper):
                    # sin²(π/2 × fraction) → 0 at boundary, 1 at interior edge
                    t = _np.sin(_np.pi * (i + 0.5) / (2.0 * n_taper)) ** 2
                    taper[i]          = t   # south boundary
                    taper[ny - 1 - i] = t   # north boundary
                f_deviation = f_deviation * taper

            f_y_1d      = _np.float64(f) + f_deviation
            self._f_y   = asarray(f_y_1d)[None, :, None]  # (1, ny, 1)
        else:
            self._beta  = None
            self._f_y   = None

    @property
    def f(self) -> float:
        """Reference Coriolis parameter f₀ (s⁻¹)."""
        return self._f0

    @property
    def mode(self) -> str:
        """Active approximation: "f_plane" or "beta_plane"."""
        return self._mode

    @property
    def beta(self):
        """β = df/dy (m⁻¹ s⁻¹), or None for f-plane."""
        return self._beta

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        # Coriolis acts on PERTURBATION velocity only (departure from
        # geostrophic background).  This assumes the background steering
        # flow is in geostrophic balance with an implied synoptic pressure
        # gradient — physically correct for a DLM steering flow.
        u_p = state.u - self._u_env
        v_p = state.v - self._v_env
        if self._mode == "f_plane":
            f = self._f0
            return Tendency(du_dt=+f * v_p, dv_dt=-f * u_p)
        else:  # beta_plane
            return Tendency(
                du_dt=+self._f_y * v_p,
                dv_dt=-self._f_y * u_p,
            )

    def set_env(self, u_env: float, v_env: float) -> None:
        """
        Update the geostrophic background reference (u_env, v_env) mid-run.

        Used by time-varying steering (mean-flow relaxation).  SAFE only when
        the caller simultaneously shifts the background flow ALREADY PRESENT in
        the state by the same increment — i.e. field and reference advance
        together.  Updating this reference ALONE (while the state's background
        stays frozen) decouples the geostrophic fix from the actual flow and
        drives a spurious Coriolis torque (the eastward drift seen in Runs 5–6,
        when only the reference was updated).  The β-plane f(y) profile does not
        depend on u_env, so nothing else needs recomputing.
        """
        self._u_env = float(u_env)
        self._v_env = float(v_env)

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)


class SurfaceDragComponent(TendencyComponent):
    """
    Bulk aerodynamic surface drag with boundary-layer depth profile.

    Acts on the PERTURBATION wind (u - u_env, v - v_env), consistent with
    CoriolisComponent and SpongeDampingComponent.

        |V'|     = √((u₀ - u_env)² + (v₀ - v_env)²)   at the lowest level (k=0)
        α(z)     = Cd · |V'| / dz · max(0, 1 − z / H_bl)
        du/dt|_k = −α_k · (u_k − u_env)
        dv/dt|_k = −α_k · (v_k − v_env)

    WHY PERTURBATION-RELATIVE (changed from V8.3):
    The V8.3 form damped the TOTAL wind toward zero, which spun down the
    background steering flow in the boundary layer over the run and
    manufactured spurious vertical shear between the dragged surface layers
    and the free-tropospheric steering flow — exactly the failure mode
    SpongeDampingComponent's docstring is written to avoid.  The steering
    flow is a maintained geostrophic background (its surface stress is
    implicitly balanced by the synoptic pressure gradient that holds it up),
    so only the vortex's DEPARTURE from it should be dragged.  Building the
    drag rate from |V'| as well makes the environment an exact fixed point:
    where the perturbation vanishes, the drag vanishes, regardless of the
    background magnitude.

    At the eyewall |V'| ≈ |V_total| (the vortex dominates the few-m/s
    background), so the angular-momentum sink that caps the runaway
    secondary intensification is unchanged:

        τ = dz / (Cd·|V'|) = 625 / (1.5e-3 · 150) ≈ 2 780 s ≈ 46 min

    e-folding at the eyewall surface — consistent with the brief-4 estimate
    that explicit drag is stable at these wind speeds (Δu per slow step well
    under 1%; no implicit treatment required).

    Parameters
    ----------
    Cd : float
        Bulk drag coefficient (dimensionless).  Default 1.5×10⁻³
        (typical open-ocean value used in TC models).
    H_bl : float
        Boundary-layer depth (m).  Drag is zero at and above this height.
        Default 1000 m.
    u_env, v_env : float
        Deep-layer-mean steering wind (m/s) toward which the dragged layers
        relax.  Pass the same values used for the vortex initialization and
        CoriolisComponent.  Default 0.0 reproduces the old total-wind
        behaviour (warm-bubble tests, component validation).
    """
    name = "surface_drag"
    stage = StepStage.SLOW

    def __init__(self, Cd: float = 1.5e-3, H_bl: float = 1000.0,
                 u_env: float = 0.0, v_env: float = 0.0) -> None:
        self._Cd    = float(Cd)
        self._H_bl  = float(H_bl)
        self._u_env = float(u_env)
        self._v_env = float(v_env)

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        u  = state.u    # (nx, ny, nz)
        v  = state.v
        z  = base.z     # (nz,) cell-centred heights

        # Grid spacing from the base-state z array (uniform grid assumed)
        dz = float(z[1] - z[0]) if len(z) > 1 else float(z[0])

        # Departure from the maintained steering flow
        u_p = u - self._u_env
        v_p = v - self._v_env

        # Bulk drag rate from the SURFACE perturbation speed: α₀ = Cd·|V'|/dz
        V_sfc     = np.sqrt(u_p[:, :, 0]**2 + v_p[:, :, 0]**2)   # (nx, ny)
        alpha_sfc = self._Cd * V_sfc / dz                        # (nx, ny) [s⁻¹]

        # Vertical weight: 1 at surface → 0 at H_bl, 0 above.  Vectorised
        # over z (the old per-k Python loop was nz separate GPU kernels/step).
        weight = np.maximum(0.0, 1.0 - z / self._H_bl)           # (nz,)
        alpha  = alpha_sfc[:, :, None] * weight[None, None, :]   # (nx, ny, nz)

        return Tendency(du_dt=-alpha * u_p, dv_dt=-alpha * v_p)

    def set_env(self, u_env: float, v_env: float) -> None:
        """Update the steering reference mid-run (time-varying steering).

        Same lockstep requirement as CoriolisComponent.set_env: only update
        this in concert with shifting the state's background flow by the same
        increment.  Since the drag rate and force are both perturbation-relative
        (u - u_env), advancing field and reference together leaves the eyewall
        drag unchanged and the environment an exact fixed point at the new DLM.
        """
        self._u_env = float(u_env)
        self._v_env = float(v_env)

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)


class IntensityCapComponent(TendencyComponent):
    """
    Perturbation-wind speed limiter (V8.4.2).

    Gently relaxes the perturbation wind |V'| = |(u-u_env, v-v_env)| toward a
    ceiling v_cap *only where it exceeds the ceiling*, preserving direction and
    leaving the background/steering flow untouched.  In magnitude:

        d|V'|/dt = -(|V'| - v_cap) / tau            for |V'| > v_cap
                 =  0                                otherwise

    realised as the direction-preserving tendency dV'/dt = -alpha·V' with
    alpha = (|V'| - v_cap) / (tau·|V'|).  So the broad outer circulation below
    v_cap is untouched and only the runaway eyewall is bled back toward v_cap.

    WHY: the f-plane translation test showed the vortex is numerically unstable
    above ~90 m/s (init 120 → peak 372–493 m/s, NaN-prone; ε and drag were
    acting as stability crutches).  Katrina's mid-run max|u|~150 rides that
    runaway, which (a) makes "intensity" unmeasurable and (b) inflates the outer
    circulation that drives the reversed-taper β-gyre.  Capping V' to a physical
    ceiling gives an honest intensity and an uncontaminated ventilation read,
    and — if the runaway was amplifying the taper deficit — may shrink the lag.

    This is a numerical SPEED LIMITER, not a buoyancy/physics replacement; label
    it as such in the paper.  Perturbation-relative ⇒ it carries a steering
    reference and MUST be advanced in lockstep with time-varying steering
    (set_env), exactly like SurfaceDragComponent.
    """

    def __init__(self, v_cap: float = 70.0, tau: float = 300.0,
                 u_env: float = 0.0, v_env: float = 0.0) -> None:
        self._v_cap = float(v_cap)
        self._tau = float(tau)
        self._u_env = float(u_env)
        self._v_env = float(v_env)

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        u_p = state.u - self._u_env
        v_p = state.v - self._v_env
        spd = np.sqrt(u_p ** 2 + v_p ** 2)
        excess = np.maximum(0.0, spd - self._v_cap)
        # alpha>0 only above the cap; preserve direction; guard tiny spd
        alpha = excess / (self._tau * np.maximum(spd, 1e-6))
        return Tendency(du_dt=-alpha * u_p, dv_dt=-alpha * v_p)

    def set_env(self, u_env: float, v_env: float) -> None:
        """Lockstep steering-reference update (see SurfaceDragComponent.set_env)."""
        self._u_env = float(u_env)
        self._v_env = float(v_env)

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)


class SpongeDampingComponent(TendencyComponent):
    """
    Rayleigh sponge damping near the top boundary.

    Damps all prognostic perturbations toward the BACKGROUND ENVIRONMENT,
    not toward zero.  This is critical when a steering flow (u_env, v_env)
    is present: damping toward zero would destroy the steering wind in the
    sponge layer, creating artificial vertical shear that shreds storm
    outflow.  With u_env and v_env set to the deep-layer mean steering
    wind, the sponge absorbs gravity waves and vortex outflow perturbations
    while leaving the environmental flow undisturbed.

    Damping profile (quadratic, zero-derivative at z_sponge):

        α(z) = α_max · [(z − z_sp) / (z_top − z_sp)]²   z > z_sp
             = 0                                          z ≤ z_sp

    Applied as perturbation-relative tendencies:

        du/dt       += −α(z) · (u − u_env)
        dv/dt       += −α(z) · (v − v_env)
        dw/dt       += −α_half(z) · w          (w_env = 0 by definition)
        dθ′/dt      += −α(z) · θ′              (θ_env already in base state)

    For warm-bubble tests and component validation, leave u_env = v_env = 0
    (default) — behaviour is identical to the previous zero-damping form.

    Parameters
    ----------
    Lz : float
        Domain height (m).  Default 10 000 m.
    alpha_max : float
        Maximum damping rate at the lid (s⁻¹).  Default 0.01 (100 s).
    sponge_fraction : float
        Fraction of domain height occupied by the sponge layer.
        Default 0.3 (top 30%, i.e. z_sp = 0.7 · Lz).
    u_env : float
        Background zonal wind (m/s) toward which u is damped.
        Set to the deep-layer mean steering wind for TC track runs.
        Default 0.0 (backward-compatible).
    v_env : float
        Background meridional wind (m/s) toward which v is damped.
        Default 0.0 (backward-compatible).
    """
    name = "sponge_damping"
    stage = StepStage.SLOW

    def __init__(
        self,
        Lz: float = 10_000.0,
        alpha_max: float = 0.01,
        sponge_fraction: float = 0.3,
        u_env: float = 0.0,
        v_env: float = 0.0,
    ) -> None:
        self._Lz          = Lz
        self._alpha_max   = alpha_max
        self._z_sponge    = Lz * (1.0 - sponge_fraction)
        self._u_env       = u_env
        self._v_env       = v_env

    def _alpha_full(self, z: np.ndarray) -> np.ndarray:
        """Damping profile on full levels (shape matches z)."""
        alpha = np.zeros_like(z)
        above = z > self._z_sponge
        if np.any(above):
            alpha[above] = self._alpha_max * (
                (z[above] - self._z_sponge)
                / (self._Lz - self._z_sponge)
            ) ** 2
        return alpha

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        z     = base.z
        alpha = self._alpha_full(z)       # (nz,)

        # Half-level alpha for w
        nz   = len(z)
        ah   = np.zeros(nz + 1)
        ah[1:-1] = 0.5 * (alpha[:-1] + alpha[1:])
        ah[0]    = alpha[0]
        ah[-1]   = alpha[-1]

        # Damp perturbations from the background environment.
        # For u and v: tendency drives (u - u_env) → 0.
        # For w and θ′: environment values are zero by construction.
        return Tendency(
            du_dt           = -alpha[None, None, :] * (state.u - self._u_env),
            dv_dt           = -alpha[None, None, :] * (state.v - self._v_env),
            dw_dt           = -ah[None, None, :]    *  state.w,
            dtheta_prime_dt = -alpha[None, None, :] *  state.theta_prime,
        )

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V, StateVar.W, StateVar.THETA_PRIME)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V, StateVar.W, StateVar.THETA_PRIME)


class HorizontalDiffusionComponent(TendencyComponent):
    """
    Second-order horizontal Laplacian diffusion for TC track simulations.

    The nonlinear self-advection of a strong vortex (u.grad(u)) on a Cartesian
    grid generates non-zero divergence every timestep.  Without diffusion this
    accumulates into an aliasing instability at ~t=1h for Cat 4/5 vortices on
    O(15 km) grids.  This component damps grid-scale noise:

      grid scale (15 km):  tau_diff = dx^2/nu_H ~  25 min   (fast, stabilising)
      vortex scale (75 km): tau_diff = Rmax^2/nu_H ~ 400 days (vortex preserved)

    Standard practice in WRF-ARW, CM1, and all operational TC models.
    Documented in the Oracle V8 methods paper.

    Parameters
    ----------
    nu_H : float  Horizontal diffusivity (m^2/s).
                  Recommended: 2e5 m^2/s for dx=15.6 km, dt=30s runs.
    Lx, Ly : float  Domain extent (m).
    nx, ny : int    Grid dimensions.
    """
    name  = "horizontal_diffusion"
    stage = StepStage.SLOW

    def __init__(self, nu_H: float, Lx: float, Ly: float,
                 nx: int, ny: int) -> None:
        self._nu_H = float(nu_H)
        self._dx2  = (Lx / nx) ** 2
        self._dy2  = (Ly / ny) ** 2

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        u  = state.u
        v  = state.v
        nu = self._nu_H
        lap_u = (
            (np.roll(u, -1, axis=0) - 2.0 * u + np.roll(u, +1, axis=0)) / self._dx2
            + (np.roll(u, -1, axis=1) - 2.0 * u + np.roll(u, +1, axis=1)) / self._dy2
        )
        lap_v = (
            (np.roll(v, -1, axis=0) - 2.0 * v + np.roll(v, +1, axis=0)) / self._dx2
            + (np.roll(v, -1, axis=1) - 2.0 * v + np.roll(v, +1, axis=1)) / self._dy2
        )
        return Tendency(du_dt=nu * lap_u, dv_dt=nu * lap_v)

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)


class HyperDiffusionComponent(TendencyComponent):
    """
    4th-order biharmonic (hyper)diffusion for TC track simulations.

    du/dt = -nu4 * nabla^4 u    (applied via two Laplacian applications)

    Scale selectivity vs 2nd-order Laplacian diffusion:

      Scale      2nd-order tau  4th-order tau
      2*dx        ~20 min         ~20 min (same by design)
      Rmax=75km   ~8 h            ~7 days  ← vortex preserved!
      2*Rmax      ~32 h           >> 1 year

    At the same grid-scale damping time, hyperdiffusion leaves the
    vortex structure essentially untouched over a 30h run.  This
    directly resolves the stability-intensity tradeoff that killed
    the vortex in runs with 2nd-order diffusion (diagnosed ensemble
    review May 2026).

    Stability criterion: nu4 < dx^4 / (64 * dt)
    For dx=15.625 km, dt=30s: nu4_max = 3.1e12 m^4/s

    Parameters
    ----------
    nu4 : float
        Biharmonic diffusivity (m^4/s).
        Recommended: 3e11 m^4/s for dx=15.6 km, dt=30 s.
    Lx, Ly : float  Domain extent (m).
    nx, ny : int    Grid dimensions.
    """
    name  = "hyper_diffusion"
    stage = StepStage.SLOW

    def __init__(self, nu4: float, Lx: float, Ly: float,
                 nx: int, ny: int) -> None:
        self._nu4  = float(nu4)
        self._dx2  = (Lx / nx) ** 2
        self._dy2  = (Ly / ny) ** 2

    def _laplacian(self, f):
        """Second-order isotropic Laplacian with periodic BCs."""
        return (
            (np.roll(f, -1, axis=0) - 2.0 * f + np.roll(f, +1, axis=0)) / self._dx2
            + (np.roll(f, -1, axis=1) - 2.0 * f + np.roll(f, +1, axis=1)) / self._dy2
        )

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        u, v = state.u, state.v
        # nabla^4 via two Laplacian applications
        bilap_u = self._laplacian(self._laplacian(u))
        bilap_v = self._laplacian(self._laplacian(v))
        # Negative: hyperdiffusion is dissipative (removes high-k energy)
        return Tendency(du_dt=-self._nu4 * bilap_u, dv_dt=-self._nu4 * bilap_v)

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)


class NewtonianCoolingComponent(TendencyComponent):
    """
    Newtonian cooling: relaxes theta' toward zero on timescale tau.

    Standard in all dry TC models (Emanuel 1986 and forward) to prevent
    runaway adiabatic cooling.  Without it, the AdvectionComponent's
    base-state lapse-rate term (-w * dtheta0/dz) continuously generates
    theta' anomalies from the secondary circulation.  In a dry model with
    no moisture to cap latent heating, this causes ~-88 K in the updraft
    column after ~8h, triggering explosive intensification via BuoyancyComponent.

    With tau=3600s (1h), theta' is bounded at:
      theta'_max ≈ w * dtheta0/dz * tau ≈ 1 m/s × 3e-3 K/m × 3600s ≈ 10 K

    This is a realistic TC warm-core amplitude and allows a modest secondary
    circulation to persist without runaway.  Physical interpretation: represents
    radiative cooling that limits warm-core growth in the absence of SST
    heat fluxes and moisture.

    Parameters
    ----------
    tau : float
        Relaxation timescale (s).  Default 3600s (1h) for TC track runs.
        Increase to weaken; decrease to enforce stricter barotropic behaviour.
    """
    name  = "newtonian_cooling"
    stage = StepStage.SLOW

    def __init__(self, tau: float = 3600.0) -> None:
        self._alpha = 1.0 / float(tau)

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        return Tendency(dtheta_prime_dt=-self._alpha * state.theta_prime)

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.THETA_PRIME,)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.THETA_PRIME,)


class DiabaticHeatingComponent(TendencyComponent):
    """
    Prescribed annular eyewall heating — an idealized θ′ source (K/s).

    A dry anelastic vortex has no diabatic driver, so a *balanced* warm core is
    inert: it sits in equilibrium and is eroded by Newtonian cooling, producing
    no secondary circulation (measured: max|w| ~ 0.1 m/s, θ′ relaxes on τ_cool).
    This component supplies the missing forcing — a fixed annular heating that
    mimics eyewall latent-heat release, sustaining θ′ against cooling and
    driving the overturning circulation (eyewall ascent, compensating
    subsidence, low-level inflow / upper-level outflow via the projection).

    This is the simplest possible driver: prescribed (state-independent), fixed
    in space, centered at the domain center.  It is NOT a moist parameterization
    and NOT in the production (barotropic-track) config — it exists for the
    buoyancy / eyewall study, where the question is how LH82's small-perturbation
    assumption fares once θ′ is held at eyewall amplitudes.

    Shape:
        Q(r, z) = Q_max · exp(−((r − r_eyewall)/w_r)²) · exp(−((z − z_peak)/w_z)²)
    annular in r (ring at the eyewall radius), mid-tropospheric in z.

    Steady-state estimate (heating ≈ Newtonian cooling): θ′_eq ≈ Q_max · τ_cool,
    so Q_max ≈ 5e-3 K/s with τ_cool = 1800 s targets θ′_eq ~ 9 K.

    Parameters
    ----------
    Q_max : float       peak heating rate (K/s).
    r_eyewall : float   radius of the heating annulus (m); ~ Rmax.
    width_r : float     radial half-width of the ring (m).
    z_peak : float      altitude of peak heating (m); ~ mid-troposphere.
    width_z : float     vertical half-width (m).
    nx, ny, nz, Lx, Ly, Lz : grid.
    x_c, y_c : annulus center (m); defaults to domain center.

    Limitation: the center is fixed.  For a translating vortex (steering on) the
    annulus would need to track the storm center — a follow-up.
    """
    name  = "diabatic_heating"
    stage = StepStage.SLOW

    def __init__(self, Q_max: float, r_eyewall: float, width_r: float,
                 z_peak: float, width_z: float,
                 nx: int, ny: int, nz: int,
                 Lx: float, Ly: float, Lz: float,
                 x_c: float | None = None, y_c: float | None = None) -> None:
        dx, dy, dz = Lx / nx, Ly / ny, Lz / nz
        x_c = x_c if x_c is not None else Lx / 2.0
        y_c = y_c if y_c is not None else Ly / 2.0
        x = (np.arange(nx) + 0.5) * dx
        y = (np.arange(ny) + 0.5) * dy
        z = (np.arange(nz) + 0.5) * dz
        X, Y = np.meshgrid(x, y, indexing="ij")
        r = np.sqrt((X - x_c) ** 2 + (Y - y_c) ** 2)             # (nx, ny)
        f_r = np.exp(-((r - r_eyewall) / width_r) ** 2)          # annular ring
        f_z = np.exp(-((z - z_peak) / width_z) ** 2)             # mid-trop
        # precompute the static source on the compute device
        self._Q = float(Q_max) * f_r[:, :, None] * f_z[None, None, :]
        self._Q_max = float(Q_max)

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        return Tendency(dtheta_prime_dt=self._Q)

    def reads(self) -> tuple[StateVar, ...]:
        return ()

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.THETA_PRIME,)


class HelmholtzDivergenceDampingComponent(TendencyComponent):
    """
    Divergence damping via Helmholtz decomposition — correct zero-vorticity version.

    PROBLEM WITH NAIVE DIVERGENCE DAMPING:
    The formula du/dt += gamma_d * d(D)/dx expands via the vector identity
    nabla(nabla.u) = nabla^2(u) + nabla x (nabla x u)
    to include a vorticity-coupling term nabla x zeta.  For the TC eyewall
    (dzeta/dy ~ 1.33e-8 m^-1 s^-1), this adds ~12 m/s/h of spurious tangential
    spin-up at gamma_d = 5e5, driving a barotropic vortex to 105 m/s (Run 11).

    CORRECT APPROACH — Helmholtz decomposition:
    Any velocity field decomposes uniquely into irrotational (divergent) and
    solenoidal (rotational) parts.  We isolate and damp ONLY the divergent part:

      1. Compute horizontal divergence: D = du/dx + dv/dy
      2. Solve the 2D horizontal Poisson equation: nabla^2_h(psi) = D  (via FFT)
      3. Apply: du/dt -= (epsilon/dt) * d(psi)/dx
                dv/dt -= (epsilon/dt) * d(psi)/dy

    Since nabla x nabla(psi) = 0 for any scalar psi, this adds EXACTLY zero
    vorticity.  The rotational vortex circulation is completely untouched.

    Parameters
    ----------
    epsilon : float
        Fraction of horizontal divergence removed per SLOW half-step.
        epsilon=0.5 removes 75% of D per full step (two half-steps of 0.5 each).
        epsilon=1.0 removes 100% of D per step (full correction each half-step).
        Recommended: 0.5 for stable TC track runs.
    Lx, Ly : float  Domain extent (m).
    nx, ny : int    Grid dimensions.
    """
    name  = "helmholtz_divergence_damping"
    stage = StepStage.SLOW

    def __init__(self, epsilon: float, Lx: float, Ly: float,
                 nx: int, ny: int) -> None:
        self._epsilon = float(epsilon)
        self._dx = Lx / nx
        self._dy = Ly / ny

        # Pre-compute wavenumbers for 2D horizontal FFT — DEVICE-resident
        # (V8.6.3).  This component previously lived entirely on the CPU
        # ("CPU for FFT"): every call pulled the full u,v fields to host,
        # ran 4 numpy FFTs, and shipped the result back — twice per step.
        # The profiler measured it at 48% of total wall (~96% of physics
        # time) and it is why the GPU backend barely beat the CPU all
        # campaign.  Identical float64 math, now on xp.
        kx = xp.fft.fftfreq(nx, d=self._dx) * (2.0 * xp.pi)
        ky = xp.fft.fftfreq(ny, d=self._dy) * (2.0 * xp.pi)
        Kx, Ky = xp.meshgrid(kx, ky, indexing="ij")                # (nx, ny)
        K2 = Kx ** 2 + Ky ** 2
        K2[0, 0] = 1.0          # avoid division by zero at k=0 (mean mode)

        self._K2_inv = xp.where(K2 > 0, -1.0 / K2, 0.0)           # (nx, ny)
        self._iKx    = 1j * Kx                                      # (nx, ny)
        self._iKy    = 1j * Ky                                      # (nx, ny)

    def compute_tendency(self, state, equation_set, staggering, base, dt):
        # FULLY SPECTRAL implementation (Five, brief 3 diagnosis):
        # The previous version computed divergence via centered FD then solved
        # spectrally.  FD and spectral operators are inconsistent at high-k,
        # so the correction was not a clean projection in the model's discrete
        # space — explaining why ε=1.0 was WORSE than ε=0.5.
        #
        # Here everything is spectral:
        #   D_hat = iKx·û + iKy·v̂  (spectral divergence — exact)
        #   psi_hat = D_hat / (-K²)  (Helmholtz solve)
        #   u_div = IFFT(iKx·psi_hat), v_div = IFFT(iKy·psi_hat)
        #
        # In spectral space ∇·(∇ψ) = k²ψ is exact, so zero vorticity holds
        # in the same discrete operators the model uses for everything else.
        u = state.u             # (nx, ny, nz) — stays on the compute device
        v = state.v

        # Spectral representation of u and v
        u_hat = xp.fft.fft2(u, axes=(0, 1))                    # (nx, ny, nz)
        v_hat = xp.fft.fft2(v, axes=(0, 1))

        # Spectral divergence: D_hat = iKx·û + iKy·v̂
        D_hat = (self._iKx[:, :, None] * u_hat
                 + self._iKy[:, :, None] * v_hat)               # (nx, ny, nz)

        # Helmholtz solve: ∇²ψ = D  →  ψ_hat = D_hat / (−K²)
        psi_hat = D_hat * self._K2_inv[:, :, None]

        # Spectral gradient of ψ → divergent velocity component
        u_div = xp.real(xp.fft.ifft2(
            self._iKx[:, :, None] * psi_hat, axes=(0, 1)))
        v_div = xp.real(xp.fft.ifft2(
            self._iKy[:, :, None] * psi_hat, axes=(0, 1)))

        # Tendency: u_new = u + du_dt*dt = u − ε·u_div
        # (already device-resident; no host round-trip — V8.6.3)
        return Tendency(
            du_dt=-(self._epsilon / dt) * u_div,
            dv_dt=-(self._epsilon / dt) * v_div,
        )

    def reads(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)

    def writes(self) -> tuple[StateVar, ...]:
        return (StateVar.U, StateVar.V)
