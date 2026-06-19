"""
Oracle V8 GridStaggering Abstraction
=====================================

GridStaggering encapsulates the choice of vertical grid staggering and
the interpolation/derivative operators that depend on it. Every operator
in V8 that touches z-derivatives or interpolates between fields uses
this abstraction rather than implementing its own staggering logic.

Concrete subclasses for V8:
    - LorenzStaggering: θ, ρ, u, v, φ on full levels (cell centers);
      w on half levels (cell faces). The default for V8.0.

    - CharneyPhillipsStaggering: θ co-located with w on half levels;
      ρ, u, v, φ on full levels. Planned for V8.x if hydrostatic
      adjustment tests reveal Lorenz computational-mode artifacts.

The point of the abstraction: switching staggering is a single object
swap. The advection, buoyancy, and projection components don't need to
know which staggering is in use — they just call interpolation and
derivative methods on the staggering object.

Key references:
    Lorenz, E. N. (1960): Energy and numerical weather prediction.
        Tellus, 12, 364-373. (Original Lorenz grid.)
    Charney, J. G., and N. A. Phillips (1953): Numerical integration
        of the quasi-geostrophic equations for barotropic and simple
        baroclinic flows. J. Meteor., 10, 71-99. (CP grid origin.)
    Arakawa, A., and Y.-J. G. Konor (1996): Vertical differencing of the
        primitive equations based on the Charney-Phillips grid in
        hybrid σ-p vertical coordinates. Mon. Wea. Rev., 124, 511-528.

V8 design choice (post-triangulation): start with Lorenz, retain the
GridStaggering abstraction so CP can replace it if hydrostatic
adjustment tests show 2Δz computational-mode artifacts in θ. The cost
of the abstraction is contained; the cost of switching after writing
non-abstracted code is large.
"""

from __future__ import annotations

import abc
from enum import Enum
from typing import Literal

import numpy as np


class LevelType(Enum):
    """
    Where a variable lives in the vertical staggering.

    FULL: cell-center levels, indexed k = 0, 1, ..., nz-1
    HALF: cell-face (interface) levels, indexed k+1/2 = 0+1/2, ..., (nz-1)+1/2
    SURFACE: at z=0 only (the bottom boundary)
    TOP: at z=Lz only (the top boundary)
    """
    FULL = "full"
    HALF = "half"
    SURFACE = "surface"
    TOP = "top"


VariableName = Literal["u", "v", "w", "theta", "theta_prime", "rho", "phi", "p_prime"]


class GridStaggering(abc.ABC):
    """
    Abstract base class for vertical staggering policies.

    A concrete staggering specifies:
        1. Where each prognostic variable lives (full vs half levels).
        2. Interpolation operators between full and half levels.
        3. Vertical derivative operators that respect the staggering.
        4. The discrete hydrostatic relation specific to this staggering.

    All operators that touch z derivatives interrogate the staggering
    object. Switching staggering means swapping this object; no other
    code changes.
    """

    name: str = "abstract"

    @abc.abstractmethod
    def level_type(self, variable: VariableName) -> LevelType:
        """Return where a variable lives in the staggering."""

    @abc.abstractmethod
    def interpolate_full_to_half(self, field_on_full: np.ndarray) -> np.ndarray:
        """
        Interpolate a field from full levels (nz values) to half levels
        (nz+1 values, the cell faces including top and bottom boundaries
        OR nz-1 internal half-levels — the concrete staggering decides).

        For Lorenz: interior half-levels are nz-1 simple averages of
        adjacent full levels. Boundary half-levels (surface and top)
        use one-sided extrapolation or boundary-condition values.

        The exact dimensionality (nz+1 vs nz-1) is staggering-specific
        and documented in each concrete subclass.
        """

    @abc.abstractmethod
    def interpolate_half_to_full(self, field_on_half: np.ndarray) -> np.ndarray:
        """
        Interpolate a field from half levels back to full levels
        (returns nz values).
        """

    @abc.abstractmethod
    def vertical_derivative(
        self,
        field: np.ndarray,
        from_level_type: LevelType,
        to_level_type: LevelType,
        dz: float,
    ) -> np.ndarray:
        """
        Compute the vertical derivative of a field, transferring it from
        one level type to another if needed.

        Common cases:
            ∂(field on full)/∂z → field on half: centered differences
                between adjacent full levels.
            ∂(field on half)/∂z → field on full: centered differences
                between adjacent half levels.
            ∂(field on full)/∂z → field on full: requires either a
                wider stencil or interpolation; staggering specifies.
        """

    @abc.abstractmethod
    def discrete_hydrostatic_relation(
        self,
        rho_bar: np.ndarray,
        gravity: float,
        dz: float,
    ) -> np.ndarray:
        """
        Return the discrete vertical pressure increment dp̄ that
        satisfies hydrostatic balance EXACTLY in this staggering.

        This is the operator the base-state construction uses to
        produce a base state in V8's notion of discrete hydrostatic
        balance, replacing the trapezoidal placeholder in
        base_states.py.

        For Lorenz with ρ on full levels and p on full levels:
            dp̄_k = -ρ̄_{k+1/2} g dz, integrated upward from p_surface
            where ρ̄_{k+1/2} is the half-level interpolation of ρ̄.

        For CP, the staggering of p and ρ may differ and the discrete
        relation has different terms.
        """

    def __repr__(self) -> str:
        return f"<GridStaggering: {self.name}>"


class LorenzStaggering(GridStaggering):
    """
    Lorenz vertical staggering.

    LOCKED ARRAY-SHAPE CONVENTION (Five's P0.3):

        Full-level fields (u, v, theta_prime, rho, p, projection_potential):
            shape = (nx, ny, nz)
            indexed by k = 0, 1, ..., nz-1 (cell centers)

        Half-level field (w):
            shape = (nx, ny, nz+1)
            w[:, :, 0]  = surface (= 0 for rigid surface)
            w[:, :, nz] = top (= 0 for rigid lid)
            w[:, :, k]  for k = 1, ..., nz-1 = interior half-levels

    The boundary values for w are stored EXPLICITLY (not implied by
    BC logic in operators). This allows operators that compute
    d(rho_bar w)/dz to apply uniform stencils across the entire array
    without special-casing boundaries.

    Known issue: Lorenz has a 2Δz stationary computational mode for
    potential temperature that interpolation cannot eliminate (Arakawa
    1972; well-documented). Our hydrostatic adjustment validation test
    is designed to detect this artifact. If it appears at threshold-
    relevant amplitude, swap to CharneyPhillipsStaggering.

    Reference: Lorenz (1960), Tellus 12, 364-373.
    """

    name = "lorenz"

    def level_type(self, variable: VariableName) -> LevelType:
        # Full levels: thermodynamic, density, horizontal momentum, projection potential
        # Half levels: vertical velocity only
        if variable == "w":
            return LevelType.HALF
        return LevelType.FULL

    def interpolate_full_to_half(self, field_on_full: np.ndarray) -> np.ndarray:
        """
        Map a field from full levels to half levels via symmetric arithmetic
        averaging at the interior faces.

        Input shape:  (..., nz)
        Output shape: (..., nz+1) — half levels including surface (index 0)
                      and lid (index nz)

        Interior faces (k = 1, ..., nz-1):
            out[..., k] = 0.5 * (field[..., k-1] + field[..., k])

        Surface (k=0) and lid (k=nz):
            Left at zero (from np.zeros initialization). Boundary handling
            is intentionally NOT baked into this generic interpolation —
            it is the caller's responsibility to assert the appropriate
            boundary condition. For LH82 with rigid surfaces, the equation
            set zeros these locations as an explicit assertion of the
            kinematic constraint w=0 → dw/dt=0 at boundaries.

        Symmetry of 0.5-averaging is mathematically required for energetic
        consistency on the Lorenz grid: it guarantees that the discrete
        buoyancy flux operator and the discrete continuity projection
        operator commute, eliminating spurious entropy production in
        long-time integrations.
        """
        # Allocate output with nz+1 in the last axis; rest of shape matches input
        out_shape = field_on_full.shape[:-1] + (field_on_full.shape[-1] + 1,)
        out = np.zeros(out_shape, dtype=field_on_full.dtype)
        # Interior faces: arithmetic mean of adjacent full levels
        out[..., 1:-1] = 0.5 * (field_on_full[..., :-1] + field_on_full[..., 1:])
        # Surface (out[..., 0]) and lid (out[..., -1]) stay at zero.
        # Callers asserting different boundary behavior must do so explicitly.
        return out

    def interpolate_half_to_full(self, field_on_half: np.ndarray) -> np.ndarray:
        # Inverse direction: average of adjacent half levels back to full.
        # Internal full levels k = 1, ..., nz-2 use 0.5*(half[k-1+1/2] + half[k+1/2]).
        # Boundary full levels (k=0, k=nz-1) use one-sided combinations or
        # boundary conditions, depending on what's being interpolated.
        raise NotImplementedError(
            "Lorenz half-to-full interpolation. Pending implementation; "
            "boundary handling will be staggering-specific."
        )

    def vertical_derivative(
        self,
        field: np.ndarray,
        from_level_type: LevelType,
        to_level_type: LevelType,
        dz: float,
    ) -> np.ndarray:
        # Centered differences across the appropriate level pairs.
        # FULL → HALF: ∂field/∂z at half-level k+1/2 = (field[k+1] - field[k]) / dz
        # HALF → FULL: ∂field/∂z at full-level k = (field[k+1/2] - field[k-1/2]) / dz
        raise NotImplementedError(
            "Lorenz vertical derivative. Pending implementation; "
            "boundary stencils will use one-sided differences."
        )

    def discrete_hydrostatic_relation(
        self,
        rho_bar: np.ndarray,
        gravity: float,
        dz: float,
    ) -> np.ndarray:
        # Lorenz with p, ρ on full levels:
        # dp̄_{k+1/2} = -ρ̄_{k+1/2} * g * dz where ρ̄_{k+1/2} = 0.5*(ρ̄_k + ρ̄_{k+1})
        # Integrating: p̄_{k+1} = p̄_k + dp̄_{k+1/2}
        raise NotImplementedError(
            "Lorenz discrete hydrostatic relation. Replaces the trapezoidal "
            "placeholder in base_states.py once this method is implemented."
        )


class CharneyPhillipsStaggering(GridStaggering):
    """
    Charney-Phillips vertical staggering.

    Variable placement:
        Full levels (k = 0, ..., nz-1):
            ρ, p, u, v, φ
        Half levels (k+1/2):
            θ, θ′, w (co-located, eliminating the 2Δz θ-w mode)
        Surface and top:
            boundary conditions for w, surface fluxes

    Advantages over Lorenz:
        - No 2Δz stationary computational mode for θ.
        - Better gravity wave dispersion at marginal grid scales.
        - Documented improvements in TC eyewall asymmetry (cited in
          Decision 2 architectural review).

    Status: planned for V8.x. V8.0 ships Lorenz as the simpler baseline,
    swap if hydrostatic adjustment validation shows artifacts.
    """

    name = "charney_phillips"

    def level_type(self, variable: VariableName) -> LevelType:
        if variable in ("w", "theta", "theta_prime"):
            return LevelType.HALF
        return LevelType.FULL

    def interpolate_full_to_half(self, field_on_full: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Planned for V8.x.")

    def interpolate_half_to_full(self, field_on_half: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Planned for V8.x.")

    def vertical_derivative(
        self,
        field: np.ndarray,
        from_level_type: LevelType,
        to_level_type: LevelType,
        dz: float,
    ) -> np.ndarray:
        raise NotImplementedError("Planned for V8.x.")

    def discrete_hydrostatic_relation(
        self,
        rho_bar: np.ndarray,
        gravity: float,
        dz: float,
    ) -> np.ndarray:
        raise NotImplementedError("Planned for V8.x.")
