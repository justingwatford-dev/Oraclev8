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

Scope of the abstraction in V8.0 (honest statement): V8.0 is Lorenz-only.
The abstraction currently captures exactly two things that are actually on
the execution path — vertical level placement (`level_type`) and the
full→half interpolation used by the buoyancy component
(`interpolate_full_to_half`). The vertical *derivative* stencils are
implemented directly inside AdvectionComponent and AnelasticProjection
(hardcoded for the Lorenz grid), NOT routed through this object. So a
staggering swap is NOT a drop-in today: CharneyPhillipsStaggering is a
forward-looking placeholder (only `level_type` is meaningful), and actually
switching to it would require refactoring those components to call the
staggering's derivative operators. The class is retained as the seam for
that future work, not as a working swap.

Key references:
    Lorenz, E. N. (1960): Energy and numerical weather prediction.
        Tellus, 12, 364-373. (Original Lorenz grid.)
    Charney, J. G., and N. A. Phillips (1953): Numerical integration
        of the quasi-geostrophic equations for barotropic and simple
        baroclinic flows. J. Meteor., 10, 71-99. (CP grid origin.)
    Arakawa, A., and Y.-J. G. Konor (1996): Vertical differencing of the
        primitive equations based on the Charney-Phillips grid in
        hybrid σ-p vertical coordinates. Mon. Wea. Rev., 124, 511-528.

V8 design choice (post-triangulation): ship Lorenz, retain a *minimal*
GridStaggering seam (level placement + full→half interpolation) so the
later move to CP — if hydrostatic-adjustment tests show 2Δz θ artifacts —
has a defined entry point. The unimplemented derivative/hydrostatic methods
that previously lived here were removed: they were never called (the
components compute Lorenz stencils inline), so keeping them as abstract
contract overstated what the abstraction delivers.
"""

from __future__ import annotations

import abc
from enum import Enum
from typing import Literal

import numpy as np


def _xp_like(arr):
    """
    Return the array module (cupy or numpy) that ``arr`` lives on, so output
    arrays are allocated on the same device as the input.  Buoyancy fields are
    device-resident under the GPU backend; allocating the interpolation output
    with bare numpy and then assigning a cupy slice into it raises
    "implicit conversion to a NumPy array is not allowed".
    """
    try:
        import cupy
        return cupy.get_array_module(arr)
    except ImportError:
        return np


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

    A concrete staggering specifies, in V8.0:
        1. Where each prognostic variable lives (full vs half levels)
           — `level_type`.
        2. The full→half interpolation used by the buoyancy component
           — `interpolate_full_to_half`.

    NOTE: vertical-derivative and discrete-hydrostatic operators are NOT
    part of this contract in V8.0 — they are implemented inline in the
    components (AdvectionComponent, AnelasticProjection) for the Lorenz
    grid. See the module docstring.
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
        # Allocate output with nz+1 in the last axis; rest of shape matches input.
        # Use the input's array module so the output lands on the same device
        # (GPU-safe: b_full is device-resident under the cupy backend).
        xp = _xp_like(field_on_full)
        out_shape = field_on_full.shape[:-1] + (field_on_full.shape[-1] + 1,)
        out = xp.zeros(out_shape, dtype=field_on_full.dtype)
        # Interior faces: arithmetic mean of adjacent full levels
        out[..., 1:-1] = 0.5 * (field_on_full[..., :-1] + field_on_full[..., 1:])
        # Surface (out[..., 0]) and lid (out[..., -1]) stay at zero.
        # Callers asserting different boundary behavior must do so explicitly.
        return out


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

    Status: forward-looking placeholder for V8.x. Only `level_type` is
    implemented (and exercised by the abstractions smoke test); the
    interpolation is not implemented. V8.0 ships Lorenz; moving to CP is a
    V8.x task that also requires routing the component vertical stencils
    through the staggering (see module docstring).
    """

    name = "charney_phillips"

    def level_type(self, variable: VariableName) -> LevelType:
        if variable in ("w", "theta", "theta_prime"):
            return LevelType.HALF
        return LevelType.FULL

    def interpolate_full_to_half(self, field_on_full: np.ndarray) -> np.ndarray:
        raise NotImplementedError("CharneyPhillipsStaggering is a V8.x placeholder.")
