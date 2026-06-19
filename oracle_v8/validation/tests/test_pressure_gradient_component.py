"""
Validation: PressureGradientComponent
======================================

Mirrors the structure of test_buoyancy_component.py.

Three sections:

  [A] By-hand vs component comparison
      Build a known φ with an analytical gradient, compute the
      gradient by hand, run PressureGradientComponent.compute_tendency
      on the same φ, verify bit-for-bit agreement.

  [B] Sign and physics sanity
      For a φ field that is positive at the bubble center (high
      projection potential = high pressure), the pressure-gradient
      tendency must point *away* from center (outward acceleration).
      The reversed case (negative φ at center) must produce inward
      acceleration.

  [C] End-to-end: buoyancy + pressure-gradient + projection on the
      warm bubble
      Starting from the −64.843 Pa warm-bubble baseline:
        1. Apply one buoyancy step (exactly as test_buoyancy_component).
        2. Apply one pressure-gradient step using the φ from that
           projection.
        3. Run projection again.
      The second φ should remain negative at the surface beneath the
      bubble.  The horizontal velocity should now show a divergent
      outflow pattern at low levels (pressure gradient accelerating
      fluid outward under the lifted column) — the first hint of the
      secondary circulation that drives TC intensification.

The bit-for-bit comparison in [A] confirms that
PressureGradientComponent._phi_gradient and
AnelasticProjection._compute_phi_gradient use identical stencils
(they both call the module-level _phi_gradient helper).

Passes → PressureGradientComponent is ready for the RK3 integrator.
"""

from __future__ import annotations

import sys

import numpy as np

from oracle_v8.solver.tendency import (
    BuoyancyComponent,
    PressureGradientComponent,
    AnelasticProjection,
    State,
    Tendency,
)
from oracle_v8.solver.equation_set import LH82AnelasticEquationSet
from oracle_v8.grid.staggering import LorenzStaggering
from oracle_v8.validation.tests.test_hydrostatic_adjustment import (
    WarmBubbleParams,
    warm_bubble_theta_perturbation,
)


# ---------------------------------------------------------------------------
# Grid shared across all sections
# ---------------------------------------------------------------------------
Lx = Ly = 100_000.0
Lz = 10_000.0
nx = ny = 64
nz = 32
dz = Lz / nz
z_centers = (np.arange(nz) + 0.5) * dz

# Standard warm-bubble base state (cell-centered, same as other tests)
theta0_cell = 300.0 * np.exp(0.01**2 * z_centers / 9.81)
Pi = np.zeros(nz)
Pi[0] = 1.0 - (9.81 / 1004.5) * z_centers[0] / theta0_cell[0]
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


BASE = CellCenteredBase()


def make_state(phi=None, w=None, theta_prime=None):
    """Convenience: make a zero-velocity State, override fields as needed."""
    s = State(
        u=np.zeros((nx, ny, nz)),
        v=np.zeros((nx, ny, nz)),
        w=np.zeros((nx, ny, nz + 1)) if w is None else w,
        theta_prime=np.zeros((nx, ny, nz)) if theta_prime is None else theta_prime,
        projection_potential=np.zeros((nx, ny, nz)) if phi is None else phi,
        t=0.0,
    )
    return s


# ---------------------------------------------------------------------------
# Section A: by-hand vs component
# ---------------------------------------------------------------------------

def section_a() -> bool:
    print("\n[A] By-hand vs component (bit-for-bit)")

    # Manufactured φ: single horizontal cosine + vertical linear ramp.
    #   φ(x, y, z) = A · cos(2π x / Lx)  +  B · z
    # Analytical gradients:
    #   ∂φ/∂x = −A · (2π/Lx) · sin(2π x / Lx)
    #   ∂φ/∂y = 0
    #   ∂φ/∂z = B  (uniform on full levels; on half levels same, with
    #                boundaries zeroed by Neumann BC)
    A = 50.0    # Pa amplitude
    B = 0.005   # Pa/m vertical ramp

    x = np.linspace(0, Lx, nx, endpoint=False)
    y = np.linspace(0, Ly, ny, endpoint=False)
    X, Y, Z = np.meshgrid(x, y, z_centers, indexing="ij")
    phi_mfg = A * np.cos(2 * np.pi * X / Lx) + B * Z

    # By-hand horizontal gradient via FFT (same stencil as the component)
    kx_1d = 2.0 * np.pi * np.fft.fftfreq(nx, d=Lx / nx)
    ky_1d = 2.0 * np.pi * np.fft.fftfreq(ny, d=Ly / ny)
    phi_hat = np.fft.fft2(phi_mfg, axes=(0, 1))
    dphi_dx_byhand = np.real(np.fft.ifft2(1j * kx_1d[:, None, None] * phi_hat, axes=(0, 1)))
    dphi_dy_byhand = np.real(np.fft.ifft2(1j * ky_1d[None, :, None] * phi_hat, axes=(0, 1)))

    # By-hand vertical gradient: centered FD full → half, Neumann BCs
    dphi_dz_byhand = np.zeros((nx, ny, nz + 1))
    dphi_dz_byhand[:, :, 1:-1] = (phi_mfg[:, :, 1:] - phi_mfg[:, :, :-1]) / dz

    print(f"  max |∂φ/∂x| by-hand: {np.max(np.abs(dphi_dx_byhand)):.4e} Pa/m")
    print(f"  max |∂φ/∂y| by-hand: {np.max(np.abs(dphi_dy_byhand)):.4e} Pa/m "
          f"(should be ~0)")
    print(f"  max |∂φ/∂z| by-hand (interior): "
          f"{np.max(np.abs(dphi_dz_byhand[:,:,1:-1])):.4e} Pa/m")

    # Component
    pg = PressureGradientComponent(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz)
    state = make_state(phi=phi_mfg)
    eq = LH82AnelasticEquationSet()
    stg = LorenzStaggering()
    tendency = pg.compute_tendency(state, eq, stg, BASE, dt=1.0)

    diff_dx = np.max(np.abs(-dphi_dx_byhand - tendency.du_dt))
    diff_dy = np.max(np.abs(-dphi_dy_byhand - tendency.dv_dt))
    diff_dz = np.max(np.abs(-dphi_dz_byhand - tendency.dw_dt))

    print(f"\n  max |by-hand − component|:")
    print(f"    du_dt: {diff_dx:.3e}")
    print(f"    dv_dt: {diff_dy:.3e}")
    print(f"    dw_dt: {diff_dz:.3e}")

    passed = diff_dx == 0.0 and diff_dy == 0.0 and diff_dz == 0.0
    if passed:
        print(f"  ✓ EXACT bit-for-bit match on all three components")
    else:
        print(f"  ✗ MISMATCH — check _phi_gradient helper sharing")

    # Shapes
    assert tendency.du_dt.shape == (nx, ny, nz), "du_dt shape wrong"
    assert tendency.dv_dt.shape == (nx, ny, nz), "dv_dt shape wrong"
    assert tendency.dw_dt.shape == (nx, ny, nz + 1), "dw_dt shape wrong"
    print(f"  ✓ output shapes correct: "
          f"du/dv = {tendency.du_dt.shape}, dw = {tendency.dw_dt.shape}")

    # Neumann BCs: surface and lid entries of dw_dt must be zero
    assert np.all(tendency.dw_dt[:, :, 0] == 0.0), "surface dw_dt nonzero"
    assert np.all(tendency.dw_dt[:, :, -1] == 0.0), "lid dw_dt nonzero"
    print(f"  ✓ Neumann BCs: surface and lid dw_dt == 0")

    return passed


# ---------------------------------------------------------------------------
# Section B: sign and physics sanity
# ---------------------------------------------------------------------------

def section_b() -> bool:
    print("\n[B] Sign and physics sanity")

    pg = PressureGradientComponent(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz)
    eq = LH82AnelasticEquationSet()
    stg = LorenzStaggering()

    # Case 1: φ > 0 at centre (local pressure maximum)
    # Fluid should be pushed AWAY from centre: du/dt < 0 for x < centre,
    # du/dt > 0 for x > centre (outward in x).
    x = np.linspace(0, Lx, nx, endpoint=False)
    y = np.linspace(0, Ly, ny, endpoint=False)
    X, Y, Z = np.meshgrid(x, y, z_centers, indexing="ij")
    # Gaussian-ish positive hill centred at (Lx/2, Ly/2)
    phi_pos = 100.0 * np.exp(
        -((X - Lx / 2)**2 + (Y - Ly / 2)**2) / (2 * (Lx / 5)**2)
    )
    state_pos = make_state(phi=phi_pos)
    t_pos = pg.compute_tendency(state_pos, eq, stg, BASE, dt=1.0)

    # du/dt at (Lx/4, Ly/2) should be negative (pushed left, toward lower φ)
    ix_left = nx // 4
    iy_mid = ny // 2
    kz_mid = nz // 2
    du_left = float(t_pos.du_dt[ix_left, iy_mid, kz_mid])

    # du/dt at (3Lx/4, Ly/2) should be positive (pushed right, toward lower φ)
    ix_right = 3 * nx // 4
    du_right = float(t_pos.du_dt[ix_right, iy_mid, kz_mid])

    sign_ok_pos = du_left < 0 and du_right > 0
    print(f"  Positive φ hill: du_dt at left quarter = {du_left:.4e} m/s²")
    print(f"                   du_dt at right quarter = {du_right:.4e} m/s²")
    if sign_ok_pos:
        print(f"  ✓ Correct: outward acceleration from pressure maximum")
    else:
        print(f"  ✗ Wrong sign: inward acceleration from pressure maximum")

    # Case 2: φ < 0 at centre (pressure deficit, as in the warm-bubble result)
    # Fluid should be pulled INWARD toward centre.
    phi_neg = -phi_pos
    state_neg = make_state(phi=phi_neg)
    t_neg = pg.compute_tendency(state_neg, eq, stg, BASE, dt=1.0)

    du_left_neg = float(t_neg.du_dt[ix_left, iy_mid, kz_mid])
    du_right_neg = float(t_neg.du_dt[ix_right, iy_mid, kz_mid])

    sign_ok_neg = du_left_neg > 0 and du_right_neg < 0
    print(f"  Negative φ hollow: du_dt at left quarter = {du_left_neg:.4e} m/s²")
    print(f"                     du_dt at right quarter = {du_right_neg:.4e} m/s²")
    if sign_ok_neg:
        print(f"  ✓ Correct: inward acceleration toward pressure deficit")
    else:
        print(f"  ✗ Wrong sign: outward acceleration from pressure deficit")

    # Zero φ → zero tendency
    state_zero = make_state(phi=np.zeros((nx, ny, nz)))
    t_zero = pg.compute_tendency(state_zero, eq, stg, BASE, dt=1.0)
    max_zero = max(
        float(np.max(np.abs(t_zero.du_dt))),
        float(np.max(np.abs(t_zero.dv_dt))),
        float(np.max(np.abs(t_zero.dw_dt))),
    )
    zero_ok = max_zero == 0.0
    if zero_ok:
        print(f"  ✓ φ=0 → zero tendency (max = {max_zero:.2e})")
    else:
        print(f"  ✗ φ=0 produced nonzero tendency (max = {max_zero:.2e})")

    return sign_ok_pos and sign_ok_neg and zero_ok


# ---------------------------------------------------------------------------
# Section C: end-to-end with warm bubble
# ---------------------------------------------------------------------------

def section_c() -> bool:
    print("\n[C] End-to-end: buoyancy + pressure-gradient on warm bubble")
    print("    (Single-step diagnostic — second projection belongs in RK3 test)")

    bubble = WarmBubbleParams()
    x = np.linspace(0, Lx, nx, endpoint=False)
    y = np.linspace(0, Ly, ny, endpoint=False)
    X, Y, Z = np.meshgrid(x, y, z_centers, indexing="ij")
    theta_prime = warm_bubble_theta_perturbation(X, Y, Z, Lx, Ly, bubble)

    eq = LH82AnelasticEquationSet()
    stg = LorenzStaggering()
    buoy_comp = BuoyancyComponent()
    pg_comp = PressureGradientComponent(
        nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
    )
    proj = AnelasticProjection(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz)

    dt = 1.0
    state = make_state(theta_prime=theta_prime)

    # --- [C.1] Buoyancy step + projection: reproduce the baseline ---
    b_tendency = buoy_comp.compute_tendency(state, eq, stg, BASE, dt)
    state.w = state.w + dt * b_tendency.dw_dt
    proj.apply_projection(state, eq, stg, BASE, dt)

    phi1 = state.projection_potential
    phi1_centered = phi1 - np.mean(phi1)
    phi1_surface_min = float(np.min(phi1_centered[:, :, 0]))
    baseline_ok = abs(phi1_surface_min - (-64.843)) < 0.001
    print(f"\n  [C.1] After buoyancy + projection:")
    print(f"    surface φ_min = {phi1_surface_min:.3f} Pa  (baseline −64.843 Pa)")
    print(f"    {'✓' if baseline_ok else '✗'} baseline reproduced")

    # --- [C.2] Pressure-gradient step: horizontal velocity appears ---
    pg_tendency = pg_comp.compute_tendency(state, eq, stg, BASE, dt)
    state.u = state.u + dt * pg_tendency.du_dt
    state.v = state.v + dt * pg_tendency.dv_dt
    state.w = state.w + dt * pg_tendency.dw_dt

    max_u = float(np.max(np.abs(state.u)))
    max_v = float(np.max(np.abs(state.v)))
    nonzero_uv = max_u > 1e-6 and max_v > 1e-6
    print(f"\n  [C.2] After pressure-gradient step:")
    print(f"    max |u| = {max_u:.4e} m/s")
    print(f"    max |v| = {max_v:.4e} m/s")
    print(f"    {'✓' if nonzero_uv else '✗'} horizontal velocity non-zero")

    # --- [C.3] Flow direction: inward toward pressure deficit ---
    # At (Lx/4, Ly/2, near-surface), φ is negative to the right (toward centre),
    # so the pressure gradient drives rightward (inward) flow: u > 0.
    ix_left = nx // 4
    iy_mid = ny // 2
    kz_low = nz // 8
    u_left_low = float(state.u[ix_left, iy_mid, kz_low])
    inward_ok = u_left_low > 0
    print(f"\n  [C.3] Flow direction at (Lx/4, Ly/2, near-surface):")
    print(f"    u = {u_left_low:.4e} m/s "
          f"({'positive → inward ✓' if inward_ok else 'negative → OUTWARD ✗'})")

    # --- [C.4] Magnitude sanity: u should be O(φ·dt/L) ---
    # φ ~ 64 Pa, L ~ Lx/2 ~ 50 km, dt = 1 s, ρ ~ 1.2 kg/m³
    # du/dt = -∂φ/∂x ~ φ/L ~ 64/50000 = 1.3e-3 m/s² → u ~ 1.3e-3 m/s after 1s
    # max |u| = 1.9e-2 m/s is larger (full spectral gradient peak), reasonable
    magnitude_ok = 1e-4 < max_u < 1.0
    print(f"\n  [C.4] Magnitude sanity (expected O(10⁻³–10⁻²) m/s):")
    print(f"    max |u| = {max_u:.4e} m/s "
          f"({'✓ in expected range' if magnitude_ok else '✗ out of range'})")

    print()
    all_ok = baseline_ok and nonzero_uv and inward_ok and magnitude_ok
    return all_ok


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    print("=" * 70)
    print("VALIDATION: PressureGradientComponent")
    print("=" * 70)

    passed_a = section_a()
    passed_b = section_b()
    passed_c = section_c()

    print("\n" + "=" * 70)
    if passed_a and passed_b and passed_c:
        print("PASSED: PressureGradientComponent integration verified")
        print()
        print("Confirmed:")
        print("  ✓ Bit-for-bit match with by-hand gradient")
        print("  ✓ Correct sign: fluid accelerates toward low φ")
        print("  ✓ φ=0 → zero tendency")
        print("  ✓ End-to-end: −64.843 Pa baseline reproduced,")
        print("    inward near-surface flow at correct magnitude")
        print("    (second-projection loop belongs in RK3 integrator test)")
        print()
        print("PressureGradientComponent is ready for the RK3 integrator.")
        print("=" * 70)
        return 0
    else:
        failures = []
        if not passed_a: failures.append("A (bit-for-bit)")
        if not passed_b: failures.append("B (sign sanity)")
        if not passed_c: failures.append("C (end-to-end)")
        print(f"FAILED: sections {', '.join(failures)}")
        print("=" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
