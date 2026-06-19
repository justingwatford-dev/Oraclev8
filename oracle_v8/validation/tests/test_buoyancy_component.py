"""
Validation: BuoyancyComponent reproduces the hand-coded buoyancy
=================================================================

This test confirms that the newly-wired BuoyancyComponent +
LH82.compute_buoyancy_tendency + LorenzStaggering.interpolate_full_to_half
chain produces identical buoyancy values to the by-hand computation in
test_v8_projection_on_warm_bubble.py.

If the by-hand b_half and the component-produced b_half match bit-for-bit,
then plugging the component into a warm-bubble step must reproduce the
-64.843 Pa surface deficit we documented as the V8 headline result.

This validates Five's recommended call chain:
    BuoyancyComponent.compute_tendency(...)
      → equation_set.compute_buoyancy_tendency(theta_prime, base, staggering)
          → b_full = g·θ'/θ̄  on full levels
          → b_half = staggering.interpolate_full_to_half(b_full)
          → b_half[..., 0] = b_half[..., -1] = 0  (rigid-w constraint)
      → wrap in Tendency(dw_dt=b_half)
"""

from __future__ import annotations

import sys

import numpy as np

from oracle_v8.solver.tendency import (
    BuoyancyComponent, AnelasticProjection, State,
)
from oracle_v8.solver.equation_set import LH82AnelasticEquationSet
from oracle_v8.grid.staggering import LorenzStaggering
from oracle_v8.validation import base_states as bs
from oracle_v8.validation.tests.test_hydrostatic_adjustment import (
    WarmBubbleParams,
    warm_bubble_theta_perturbation,
    expected_pressure_deficit_at_center,
)


def main() -> int:
    print("=" * 70)
    print("VALIDATION: BuoyancyComponent reproduces by-hand buoyancy")
    print("=" * 70)

    # Same setup as test_v8_projection_on_warm_bubble.py
    Lx = Ly = 100_000.0
    Lz = 10_000.0
    nx = ny = 64
    nz = 32
    dz = Lz / nz
    z_centers = (np.arange(nz) + 0.5) * dz

    # Cell-centered base state (matches the working warm-bubble test exactly)
    theta0_cell = 300.0 * np.exp(0.01**2 * z_centers / 9.81)
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

    class CellCenteredBase:
        z = z_centers
        rho0 = rho0_cell
        theta0 = theta0_cell

    base = CellCenteredBase()

    # Warm bubble
    bubble = WarmBubbleParams()
    x = np.linspace(0, Lx, nx, endpoint=False)
    y = np.linspace(0, Ly, ny, endpoint=False)
    X, Y, Z = np.meshgrid(x, y, z_centers, indexing="ij")
    theta_prime = warm_bubble_theta_perturbation(X, Y, Z, Lx, Ly, bubble)

    print(f"\nSetup: 64x64x32 grid, B&F warm bubble (2 K, 10 km radius, 2 km altitude)")
    print(f"  max θ′ = {theta_prime.max():.4f} K")

    # ------------------------------------------------------------------
    # By-hand buoyancy (copied verbatim from test_v8_projection_on_warm_bubble.py)
    # ------------------------------------------------------------------
    g = 9.81
    buoyancy_full_byhand = g * theta_prime / theta0_cell[None, None, :]
    buoyancy_half_byhand = np.zeros((nx, ny, nz + 1))
    buoyancy_half_byhand[:, :, 1:-1] = 0.5 * (
        buoyancy_full_byhand[:, :, :-1] + buoyancy_full_byhand[:, :, 1:]
    )

    print(f"\n[By-hand buoyancy]")
    print(f"  max |b_full| = {np.max(np.abs(buoyancy_full_byhand)):.6e} m/s²")
    print(f"  max |b_half| (interior) = "
          f"{np.max(np.abs(buoyancy_half_byhand[:, :, 1:-1])):.6e} m/s²")
    print(f"  b_half[..., 0] (surface) = {buoyancy_half_byhand[0,0,0]:.6e}")
    print(f"  b_half[..., -1] (lid) = {buoyancy_half_byhand[0,0,-1]:.6e}")

    # ------------------------------------------------------------------
    # Component-produced buoyancy (the new chain)
    # ------------------------------------------------------------------
    state = State(
        u=np.zeros((nx, ny, nz)),
        v=np.zeros((nx, ny, nz)),
        w=np.zeros((nx, ny, nz + 1)),
        theta_prime=theta_prime,
        projection_potential=np.zeros((nx, ny, nz)),
        t=0.0,
    )
    eq_set = LH82AnelasticEquationSet()
    staggering = LorenzStaggering()
    buoyancy_comp = BuoyancyComponent()

    tendency = buoyancy_comp.compute_tendency(
        state, eq_set, staggering, base, dt=1.0,
    )
    buoyancy_half_component = tendency.dw_dt

    print(f"\n[Component-produced buoyancy]")
    print(f"  shape: {buoyancy_half_component.shape}")
    print(f"  max |b_half| (interior) = "
          f"{np.max(np.abs(buoyancy_half_component[:, :, 1:-1])):.6e} m/s²")
    print(f"  b_half[..., 0] (surface) = {buoyancy_half_component[0,0,0]:.6e}")
    print(f"  b_half[..., -1] (lid) = {buoyancy_half_component[0,0,-1]:.6e}")

    # ------------------------------------------------------------------
    # Bit-for-bit comparison
    # ------------------------------------------------------------------
    diff = buoyancy_half_byhand - buoyancy_half_component
    max_abs_diff = float(np.max(np.abs(diff)))

    print(f"\n[Bit-for-bit comparison]")
    print(f"  max |by_hand - component| = {max_abs_diff:.3e}")

    if max_abs_diff == 0.0:
        print(f"  ✓ EXACT MATCH — component reproduces by-hand buoyancy")
        comparison_passed = True
    elif max_abs_diff < 1e-15:
        print(f"  ✓ machine-precision match — component reproduces by-hand "
              f"buoyancy modulo float rounding")
        comparison_passed = True
    else:
        print(f"  ✗ MISMATCH — investigate")
        comparison_passed = False

    # ------------------------------------------------------------------
    # End-to-end: BuoyancyComponent + AnelasticProjection on warm bubble
    # ------------------------------------------------------------------
    print(f"\n[End-to-end: applying buoyancy step + V8 projection]")
    dt = 1.0
    # Run a manual buoyancy step using the component
    state.w = state.w + dt * buoyancy_half_component
    proj = AnelasticProjection(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz)
    proj.apply_projection(state, eq_set, staggering, base, dt)

    phi_centered = state.projection_potential - np.mean(state.projection_potential)
    xc_idx, yc_idx = nx // 2, ny // 2
    phi_surface = float(phi_centered[xc_idx, yc_idx, 0])
    surface_min = float(np.min(phi_centered[:, :, 0]))

    print(f"  φ at surface beneath bubble: {phi_surface:.3f} Pa")
    print(f"  deepest surface φ:           {surface_min:.3f} Pa")
    print(f"  expected (from baseline):    -64.843 Pa")

    target = -64.843
    end_to_end_match = abs(surface_min - target) < 0.001

    if end_to_end_match:
        print(f"  ✓ MATCHES baseline -64.843 Pa to 3 decimal places")
    else:
        print(f"  ✗ DIFFERS from baseline ({abs(surface_min - target):.3f} Pa off)")

    print("\n" + "=" * 70)
    if comparison_passed and end_to_end_match:
        print("PASSED: BuoyancyComponent integration verified end-to-end")
        print()
        print("The call chain works:")
        print("  BuoyancyComponent → LH82.compute_buoyancy_tendency")
        print("    → LorenzStaggering.interpolate_full_to_half")
        print("    → rigid-w boundary zeroing")
        print("  → Tendency(dw_dt=b_half)")
        print()
        print("Bit-for-bit reproduction of the by-hand result confirms")
        print("Five's recommended architecture is in place.")
        print("=" * 70)
        return 0
    else:
        print("FAILED — investigate before proceeding to PressureGradientComponent")
        print("=" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
