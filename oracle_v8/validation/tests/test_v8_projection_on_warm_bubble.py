"""
Projection-Driven Pressure Structure Test
==========================================

The V7-vs-V8 distinguishing test, run with what we have so far.

V8's full hydrostatic adjustment requires buoyancy + pressure-gradient +
projection working together over many timesteps. Buoyancy and pressure-
gradient aren't implemented yet, so we can't do the full integration.
But we CAN test whether the projection step alone produces the right
pressure structure when handed a state with a buoyancy-induced velocity
divergence.

Procedure:

  1. Build a hydrostatically balanced base state (constant-N dry).
  2. Initialize a warm bubble (B&F 2002 setup).
  3. Manually compute a single buoyancy step: dw/dt = g·θ′/θ̄, integrate
     w forward by Δt = 1 second. This produces a divergent velocity
     field that V7 by construction would NOT correct (no hydrostatic
     mass evacuation).
  4. Hand this divergent state to V8's AnelasticProjection.
  5. Examine the resulting projection_potential φ.

What we expect:

  - V7 (Boussinesq): if we did this on V7, the projection would still
    enforce ∇·u = 0, but the corresponding φ would be the standard
    Boussinesq pressure, with NO contribution from the warm-column
    mass deficit. φ at the surface beneath the bubble: small, ~0 Pa.

  - V8 (LH82 anelastic): the projection enforces ∇·(ρ̄u) = 0. With a
    warm column having reduced effective density, the constraint
    forces hydrostatic mass evacuation, and φ at the surface beneath
    the bubble develops a NEGATIVE deficit on the order of tens to
    hundreds of Pa.

This is not the full hydrostatic adjustment test; that needs the time
loop with buoyancy + pressure-gradient working together. But it IS the
core mechanism the full adjustment depends on, and it's exactly the
piece V7's Boussinesq core was structurally incapable of producing.

If φ comes out negative, with reasonable magnitude, and centered under
the bubble — that's strong evidence V8's anelastic core is doing the
thing V7 couldn't.
"""

from __future__ import annotations

import sys

import numpy as np

from oracle_v8.solver.poisson import VariableCoefficientPoissonSolver
from oracle_v8.solver.tendency import AnelasticProjection, State
from oracle_v8.validation import base_states as bs
from oracle_v8.validation.tests.test_hydrostatic_adjustment import (
    WarmBubbleParams,
    warm_bubble_theta_perturbation,
    expected_pressure_deficit_at_center,
)


def main() -> int:
    print("=" * 70)
    print("V8 PROJECTION: PRESSURE-DEFICIT RESPONSE TEST")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Setup: small domain, fast test
    # ------------------------------------------------------------------
    Lx = Ly = 100_000.0    # 100 km
    Lz = 10_000.0          # 10 km vertical
    nx = ny = 64
    nz = 32
    dx = Lx / nx
    dy = Ly / ny
    dz = Lz / nz

    # Cell-centered z grid (matching the Poisson solver's convention)
    z_centers = (np.arange(nz) + 0.5) * dz

    # Build base state on the cell-centered z grid.
    # base_states.constant_N_dry_base_state expects z[0]=0 currently;
    # for this test we work around by giving it the cell-centered grid
    # and using its outputs for full-level fields. This is acceptable
    # because the test is about the projection's response, not about
    # base-state thermodynamic closure.
    base = bs.constant_N_dry_base_state(
        z=np.linspace(0, Lz, nz, endpoint=False),  # placeholder z grid
        N=0.01,
        theta_surface=300.0,
        staggering=bs.GridStaggering.UNSTAGGERED_PLACEHOLDER,
    )

    # Re-evaluate base-state fields on the actual cell-centered z to
    # match the Poisson solver's grid layout.
    theta0_cell = 300.0 * np.exp(0.01**2 * z_centers / 9.81)
    # Crude trapezoidal Π integration on the cell-centered grid:
    Pi = np.zeros(nz)
    Pi[0] = 1.0 - (9.81 / 1004.5) * (z_centers[0]) / theta0_cell[0]
    for k in range(nz - 1):
        dz_local = z_centers[k + 1] - z_centers[k]
        Pi[k + 1] = Pi[k] - (9.81 / 1004.5) * (dz_local / 2.0) * (
            1.0 / theta0_cell[k] + 1.0 / theta0_cell[k + 1]
        )
    p0_cell = 100_000.0 * Pi ** (1004.5 / 287.04)
    T0_cell = theta0_cell * Pi
    rho0_cell = p0_cell / (287.04 * T0_cell)

    # Wrap into a base-state-like object the projection can consume
    class CellCenteredBase:
        z = z_centers
        rho0 = rho0_cell
        theta0 = theta0_cell

    base_cc = CellCenteredBase()

    print(f"\nBase state: constant-N (N=0.01), θ_surface=300 K")
    print(f"  ρ̄(0)   = {rho0_cell[0]:.4f} kg/m³")
    print(f"  ρ̄(top) = {rho0_cell[-1]:.4f} kg/m³  (factor {rho0_cell[0]/rho0_cell[-1]:.2f}× density variation)")

    # ------------------------------------------------------------------
    # Warm bubble at default B&F parameters
    # ------------------------------------------------------------------
    bubble = WarmBubbleParams()  # 2 K, 10 km radius, 2 km altitude
    x = np.linspace(0, Lx, nx, endpoint=False)
    y = np.linspace(0, Ly, ny, endpoint=False)
    X, Y, Z = np.meshgrid(x, y, z_centers, indexing="ij")
    theta_prime = warm_bubble_theta_perturbation(X, Y, Z, Lx, Ly, bubble)
    print(f"\nWarm bubble: max θ′ = {theta_prime.max():.2f} K, "
          f"radius = {bubble.xr_m/1000:.0f} km, "
          f"center z = {bubble.zc_m/1000:.1f} km")

    expected_deficit = expected_pressure_deficit_at_center(base_cc, bubble)
    print(f"\nExpected hydrostatic-adjustment pressure deficit at surface")
    print(f"  beneath bubble (analytical column integral): {expected_deficit:.2f} Pa")

    # ------------------------------------------------------------------
    # Manually apply one buoyancy step: dw/dt = g · θ′/θ̄
    # ------------------------------------------------------------------
    # b is on full levels (where θ′ lives in our placeholder unstaggered
    # convention). For the locked Lorenz w-shape (nx, ny, nz+1), we
    # interpolate b to the half levels where w lives.
    # For this single-step test, we use a small Δt to get a velocity
    # divergence whose magnitude we can interpret.
    dt = 1.0  # 1 second of buoyancy acceleration
    g = 9.81

    buoyancy_full = g * theta_prime / theta0_cell[None, None, :]   # m/s² on full levels
    # Interpolate to half levels: simple averaging of adjacent full levels.
    # Half-level k+1/2 (k = 0..nz-2) sits between full levels k and k+1.
    # Boundary half-levels (k=0 surface, k=nz lid) keep w=0 by rigid BC.
    buoyancy_half = np.zeros((nx, ny, nz + 1))
    buoyancy_half[:, :, 1:-1] = 0.5 * (
        buoyancy_full[:, :, :-1] + buoyancy_full[:, :, 1:]
    )
    # Surface and lid: rigid BC, w stays zero, so we leave buoyancy_half
    # at 0 there (a buoyancy force on a zero-velocity rigid boundary
    # is suppressed by the BC).

    # Initial state: rest + θ′ + one buoyancy kick of duration dt
    state = State(
        u=np.zeros((nx, ny, nz)),
        v=np.zeros((nx, ny, nz)),
        w=buoyancy_half * dt,         # divergent velocity from buoyancy
        theta_prime=theta_prime,
        projection_potential=np.zeros((nx, ny, nz)),
        t=dt,
    )

    print(f"\nProvisional state after one buoyancy step (Δt = {dt} s):")
    print(f"  max w (interior) = {np.max(np.abs(state.w[:, :, 1:-1])):.4e} m/s")
    print(f"  max θ′           = {state.theta_prime.max():.4f} K")
    print(f"  velocity is divergent: ∇·(ρ̄u) ≠ 0 by construction")

    # ------------------------------------------------------------------
    # Apply V8 projection
    # ------------------------------------------------------------------
    proj = AnelasticProjection(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz)
    proj.apply_projection(state, equation_set=None, staggering=None,
                          base=base_cc, dt=dt)

    print(f"\nProjection diagnostics:")
    print(f"  compatibility residual: {proj.last_solve.compatibility_residual:.3e}")
    print(f"  discrete operator residual: {proj.last_solve.discrete_operator_residual:.3e}")

    # ------------------------------------------------------------------
    # Examine the resulting projection potential
    # ------------------------------------------------------------------
    phi = state.projection_potential

    # Bubble center column indices
    xc_idx = nx // 2
    yc_idx = ny // 2

    # Gauge-correct: subtract mean
    phi_centered = phi - np.mean(phi)

    phi_at_surface_below_bubble = phi_centered[xc_idx, yc_idx, 0]
    phi_at_top_above_bubble = phi_centered[xc_idx, yc_idx, -1]
    phi_at_bubble_center = phi_centered[xc_idx, yc_idx, nz // 5]  # ~2 km

    print(f"\nProjection potential φ at bubble center column "
          f"(gauge-corrected):")
    print(f"  φ at surface (z≈0):       {phi_at_surface_below_bubble:9.3f} Pa")
    print(f"  φ at bubble center (~2 km): {phi_at_bubble_center:9.3f} Pa")
    print(f"  φ at top (z≈10 km):       {phi_at_top_above_bubble:9.3f} Pa")

    # Find the location of the most negative φ
    surface_phi = phi_centered[:, :, 0]
    min_idx = np.unravel_index(np.argmin(surface_phi), surface_phi.shape)
    surface_min = surface_phi[min_idx]
    horizontal_offset_km = (
        np.sqrt((min_idx[0] * dx - Lx / 2) ** 2
                + (min_idx[1] * dy - Ly / 2) ** 2)
        / 1000.0
    )
    print(f"\nDeepest surface φ (most negative): {surface_min:.3f} Pa")
    print(f"  located {horizontal_offset_km:.1f} km from bubble center")
    print(f"  bubble radius: {bubble.xr_m / 1000:.1f} km")

    # ------------------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("INTERPRETATION")
    print("=" * 70)

    is_negative = phi_at_surface_below_bubble < 0
    is_collocated = horizontal_offset_km < bubble.xr_m / 1000.0
    is_substantial = abs(phi_at_surface_below_bubble) > 0.1  # > 0.1 Pa

    if is_negative and is_substantial and is_collocated:
        print("✓ V8 projection produced a NEGATIVE pressure deficit,")
        print("  collocated with the warm bubble, with substantial magnitude.")
        print("  This is the core hydrostatic-adjustment response that")
        print("  V7's Boussinesq core was structurally incapable of producing.")
    else:
        print("⚠ Result needs interpretation:")
        if not is_negative:
            print(f"  - φ at surface beneath bubble is not negative "
                  f"({phi_at_surface_below_bubble:.3f} Pa)")
        if not is_substantial:
            print(f"  - magnitude is small ({abs(phi_at_surface_below_bubble):.3f} Pa)")
        if not is_collocated:
            print(f"  - deficit is {horizontal_offset_km:.1f} km from bubble center")

    print()
    print("Caveats:")
    print(" - This is ONE projection step on a manually-divergent velocity,")
    print("   not the full time-evolved hydrostatic adjustment.")
    print(" - The full adjustment requires buoyancy + pressure-gradient +")
    print("   projection coupled over many timesteps; those components")
    print("   are still stubs.")
    print(" - The expected deficit ({:.1f} Pa) is the equilibrium value;".format(
        expected_deficit))
    print("   one step of buoyancy with one projection won't reach equilibrium.")
    print(" - This test demonstrates V8's mechanism, not V8's full response.")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
