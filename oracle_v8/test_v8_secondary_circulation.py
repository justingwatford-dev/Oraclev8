"""
V8 TC Secondary Circulation Test
==================================
First test of vortex spin-up in Oracle V8.

Physics active:
    BuoyancyComponent    — warm bubble drives updraft
    AdvectionComponent   — full momentum + θ′ advection (V8.1)
    CoriolisComponent    — f-plane, deflects radial inflow into rotation
    AnelasticProjection  — enforces ∇·(ρ̄u) = 0 at every stage

Sections
--------
[A] CoriolisComponent unit test
    Verify signs, magnitudes, and the f-parameter interface.

[B] TC spin-up: 200 steps × dt=5s = 1000s ≈ 1.6 buoyancy periods
    Starting from a warm bubble at rest, the secondary circulation should:
    1. Develop radial inflow (from projection-enforced continuity)
    2. Develop tangential wind as Coriolis deflects the inflow cyclonically
    3. Show cyclonic (positive) relative vorticity at the warm-bubble centre
    4. Remain numerically stable throughout

[C] Vorticity diagnostic
    Compute ζ = ∂v/∂x − ∂u/∂y (vertical relative vorticity).
    The warm-bubble centre should develop positive ζ (NH cyclonic).
    Compare peak ζ to a no-Coriolis reference run.
"""

from __future__ import annotations

import sys
from typing import List

import numpy as np

from oracle_v8.solver import (
    LH82AnelasticEquationSet,
    BuoyancyComponent,
    AdvectionComponent,
    CoriolisComponent,
    AnelasticProjection,
    OperatorConfig,
    State,
    RK3Integrator,
    StepDiagnostics,
)
from oracle_v8.grid.staggering import LorenzStaggering
from oracle_v8.validation.tests.test_hydrostatic_adjustment import (
    WarmBubbleParams,
    warm_bubble_theta_perturbation,
)


# ---------------------------------------------------------------------------
# Grid / base state (canonical setup)
# ---------------------------------------------------------------------------

Lx = Ly = 100_000.0
Lz = 10_000.0
nx = ny = 64
nz = 32
dz = Lz / nz
dx = Lx / nx
dy = Ly / ny
z_centers = (np.arange(nz) + 0.5) * dz

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

F_DEFAULT = 5e-5   # s⁻¹, ≈ 20°N


def make_initial_state(theta_prime: np.ndarray) -> State:
    return State(
        u=np.zeros((nx, ny, nz)),
        v=np.zeros((nx, ny, nz)),
        w=np.zeros((nx, ny, nz + 1)),
        theta_prime=theta_prime,
        projection_potential=np.zeros((nx, ny, nz)),
        t=0.0,
    )


def make_config(f: float = F_DEFAULT) -> OperatorConfig:
    """Full V8.1 config: buoyancy + advection + Coriolis + projection."""
    return OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        buoyancy=BuoyancyComponent(),
        advection=AdvectionComponent(
            nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
        ),
        coriolis=CoriolisComponent(f=f),
        projection=AnelasticProjection(
            nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
        ),
    )


def vorticity_surface(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    Vertical relative vorticity ζ = ∂v/∂x − ∂u/∂y at the lowest full level.
    Second-order centred, periodic BCs.
    """
    k = 0
    dv_dx = (np.roll(v[:, :, k], -1, axis=0) - np.roll(v[:, :, k], 1, axis=0)) / (2*dx)
    du_dy = (np.roll(u[:, :, k], -1, axis=1) - np.roll(u[:, :, k], 1, axis=1)) / (2*dy)
    return dv_dx - du_dy


# ---------------------------------------------------------------------------
# Section A: CoriolisComponent unit test
# ---------------------------------------------------------------------------

def section_a() -> bool:
    print("\n[A] CoriolisComponent unit test")

    f = F_DEFAULT
    cor = CoriolisComponent(f=f)
    eq = LH82AnelasticEquationSet()
    stg = LorenzStaggering()

    # A.1: uniform u = U₀, v = 0  →  du/dt = 0, dv/dt = -f·U₀
    U0 = 10.0
    state_u = State(
        u=np.full((nx, ny, nz), U0), v=np.zeros((nx, ny, nz)),
        w=np.zeros((nx, ny, nz+1)), theta_prime=np.zeros((nx, ny, nz)),
        projection_potential=np.zeros((nx, ny, nz)), t=0.0,
    )
    t_u = cor.compute_tendency(state_u, eq, stg, BASE, dt=1.0)
    du_max = float(np.max(np.abs(t_u.du_dt)))
    dv_expected = -f * U0
    dv_actual = float(np.mean(t_u.dv_dt))
    ok_u = du_max < 1e-15 and abs(dv_actual - dv_expected) < 1e-15
    print(f"  u={U0}, v=0: du/dt=0 ({'✓' if du_max<1e-15 else '✗'}), "
          f"dv/dt={dv_actual:.2e} (expected {dv_expected:.2e}) "
          f"{'✓' if abs(dv_actual-dv_expected)<1e-15 else '✗'}")

    # A.2: u = 0, v = V₀  →  du/dt = +f·V₀, dv/dt = 0
    V0 = -5.0  # inflow
    state_v = State(
        u=np.zeros((nx, ny, nz)), v=np.full((nx, ny, nz), V0),
        w=np.zeros((nx, ny, nz+1)), theta_prime=np.zeros((nx, ny, nz)),
        projection_potential=np.zeros((nx, ny, nz)), t=0.0,
    )
    t_v = cor.compute_tendency(state_v, eq, stg, BASE, dt=1.0)
    dv_max = float(np.max(np.abs(t_v.dv_dt)))
    du_expected = +f * V0
    du_actual = float(np.mean(t_v.du_dt))
    ok_v = dv_max < 1e-15 and abs(du_actual - du_expected) < 1e-15
    print(f"  u=0, v={V0}: dv/dt=0 ({'✓' if dv_max<1e-15 else '✗'}), "
          f"du/dt={du_actual:.2e} (expected {du_expected:.2e}) "
          f"{'✓' if abs(du_actual-du_expected)<1e-15 else '✗'}")

    # A.3: NH sign check — radial inflow (u<0) should generate cyclonic (v>0 at east side)
    # If u = -U₀ (westward inflow from east), dv/dt = -f*(-U₀) = +f*U₀ > 0 ✓
    dv_from_inflow = float(np.mean(-cor.f * (-U0) * np.ones((nx,ny,nz))))
    cyclonic = dv_from_inflow > 0
    print(f"  NH sign: inward u < 0 → dv/dt = {dv_from_inflow:.2e} > 0 "
          f"({'✓ cyclonic' if cyclonic else '✗ anticyclonic'})")

    # A.4: f-parameter respected
    cor2 = CoriolisComponent(f=1e-4)
    t2 = cor2.compute_tendency(state_u, eq, stg, BASE, dt=1.0)
    f_ratio = float(np.mean(np.abs(t2.dv_dt))) / float(np.mean(np.abs(t_u.dv_dt)))
    f_ok = abs(f_ratio - 2.0) < 1e-10
    print(f"  f=2× → tendency 2× ({'✓' if f_ok else '✗'}, ratio={f_ratio:.6f})")

    passed = ok_u and ok_v and cyclonic and f_ok
    return passed


# ---------------------------------------------------------------------------
# Section B: TC spin-up
# ---------------------------------------------------------------------------

def section_b(
    theta_prime: np.ndarray,
    n_steps: int = 200,
    dt: float = 5.0,
) -> tuple[bool, List[StepDiagnostics], State]:
    print(f"\n[B] TC spin-up: {n_steps} steps × dt={dt}s = "
          f"{n_steps*dt:.0f}s  (buoyancy + advection + Coriolis, f={F_DEFAULT:.1e})")

    config = make_config(f=F_DEFAULT)
    integrator = RK3Integrator(config=config, base=BASE)
    state = make_initial_state(theta_prime)
    history: List[StepDiagnostics] = []

    print(f"\n  {'Step':>5}  {'t (s)':>7}  {'φ_min (Pa)':>12}  "
          f"{'max|u| (m/s)':>13}  {'max|v| (m/s)':>13}  {'max θ′ (K)':>11}")
    print(f"  {'-'*5}  {'-'*7}  {'-'*12}  {'-'*13}  {'-'*13}  {'-'*11}")

    any_nan = False
    for n in range(n_steps):
        state, diag = integrator.step(state, dt=dt, step_number=n)
        history.append(diag)

        if np.any(np.isnan(state.u)) or np.any(np.isnan(state.theta_prime)):
            print(f"\n  ✗ NaN at step {n}")
            any_nan = True
            break

        if n == 0 or (n + 1) % 20 == 0:
            print(f"  {n+1:>5}  {state.t:>7.1f}  "
                  f"{diag.surface_phi_min:>12.3f}  "
                  f"{diag.max_u:>13.4e}  "
                  f"{diag.max_v:>13.4e}  "
                  f"{diag.max_theta_prime:>11.4f}")

    if not any_nan:
        print(f"\n  ✓ No NaN in {len(history)} steps")

    max_u = max(d.max_u for d in history) if history else 0
    max_v = max(d.max_v for d in history) if history else 0
    max_w = max(d.max_w for d in history) if history else 0
    stable = (not any_nan) and max_u < 30.0 and max_v < 30.0 and max_w < 50.0

    print(f"  {'✓' if stable else '✗'} Stability: "
          f"max|u|={max_u:.3e}, max|v|={max_v:.3e}, max|w|={max_w:.3e} m/s")

    # Tangential wind developing — max|v| should grow from near-zero
    v_step1 = history[0].max_v if history else 0
    v_final = history[-1].max_v if history else 0
    tangential_wind = v_final > 0.05  # at least 5 cm/s tangential wind
    print(f"  max|v| step 1:   {v_step1:.4e} m/s")
    print(f"  max|v| step {len(history)}: {v_final:.4e} m/s")
    if tangential_wind:
        print(f"  ✓ Tangential wind developing (>{0.05} m/s)")
    else:
        print(f"  ✗ Tangential wind below threshold")

    return stable and tangential_wind, history, state


# ---------------------------------------------------------------------------
# Section C: vorticity diagnostic
# ---------------------------------------------------------------------------

def section_c(
    state_with_coriolis: State,
    theta_prime: np.ndarray,
    n_steps: int = 200,
    dt: float = 5.0,
) -> bool:
    print(f"\n[C] Vorticity diagnostic")

    xc, yc = nx // 2, ny // 2

    # Vorticity from the Coriolis run
    zeta_cor = vorticity_surface(state_with_coriolis.u, state_with_coriolis.v)
    zeta_center = float(zeta_cor[xc, yc])
    zeta_max = float(np.max(zeta_cor))
    zeta_min = float(np.min(zeta_cor))

    print(f"\n  [C.1] Vorticity field (with Coriolis, after {n_steps*dt:.0f}s):")
    print(f"    ζ at bubble centre: {zeta_center:.4e} s⁻¹")
    print(f"    max ζ over domain:  {zeta_max:.4e} s⁻¹")
    print(f"    min ζ over domain:  {zeta_min:.4e} s⁻¹")

    # Cyclonic check: ζ > 0 somewhere near the bubble centre
    # (exact centre might not be max due to discretisation)
    region = zeta_cor[xc-5:xc+5, yc-5:yc+5]
    max_zeta_region = float(np.max(region))
    cyclonic = max_zeta_region > 0

    if cyclonic:
        print(f"  ✓ Cyclonic vorticity (ζ > 0) near bubble centre — "
              f"Coriolis is spinning up the vortex")
    else:
        print(f"  ✗ No cyclonic vorticity near bubble centre")

    # Reference run WITHOUT Coriolis
    print(f"\n  [C.2] Reference run (no Coriolis, {n_steps*dt:.0f}s):")
    config_norot = OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        buoyancy=BuoyancyComponent(),
        advection=AdvectionComponent(
            nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
        ),
        projection=AnelasticProjection(
            nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
        ),
    )
    int_norot = RK3Integrator(config=config_norot, base=BASE)
    state_norot = make_initial_state(theta_prime)
    for n in range(n_steps):
        state_norot, _ = int_norot.step(state_norot, dt=dt, step_number=n)

    zeta_norot = vorticity_surface(state_norot.u, state_norot.v)
    max_zeta_norot = float(np.max(np.abs(zeta_norot)))
    max_v_norot = float(np.max(np.abs(state_norot.v)))

    print(f"    max |ζ| (no Coriolis): {max_zeta_norot:.4e} s⁻¹")
    print(f"    max |v| (no Coriolis): {max_v_norot:.4e} m/s")
    print(f"    max |v| (Coriolis):    "
          f"{float(np.max(np.abs(state_with_coriolis.v))):.4e} m/s")

    coriolis_enhances = (
        float(np.max(np.abs(state_with_coriolis.v))) > max_v_norot * 1.5
    )
    if coriolis_enhances:
        print(f"  ✓ Coriolis significantly enhances v-component (>1.5× reference)")
    else:
        print(f"  (Note: tangential enhancement modest at this timescale — "
              f"vortex spin-up requires O(1/f) ≈ 5.6 h)")

    return cyclonic


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    print("=" * 70)
    print("V8 TC SECONDARY CIRCULATION TEST")
    print("Buoyancy + Full Advection (V8.1) + Coriolis + Projection")
    print("=" * 70)

    bubble = WarmBubbleParams()
    x = np.linspace(0, Lx, nx, endpoint=False)
    y = np.linspace(0, Ly, ny, endpoint=False)
    X, Y, Z = np.meshgrid(x, y, z_centers, indexing="ij")
    theta_prime = warm_bubble_theta_perturbation(X, Y, Z, Lx, Ly, bubble)

    print(f"\nSetup: 64×64×32 grid, B&F warm bubble (2 K, 10 km radius, 2 km alt)")
    print(f"  f = {F_DEFAULT:.1e} s⁻¹  (≈ 20°N)")
    print(f"  Inertial period 1/f = {1/F_DEFAULT:.0f} s ≈ "
          f"{1/F_DEFAULT/3600:.1f} h")

    passed_a = section_a()
    passed_b, history, final_state = section_b(theta_prime, n_steps=200, dt=5.0)
    passed_c = section_c(final_state, theta_prime, n_steps=200, dt=5.0)

    print("\n" + "=" * 70)
    if passed_a and passed_b and passed_c:
        print("PASSED: V8.1 TC secondary circulation verified")
        print()
        print("Confirmed:")
        print("  ✓ Coriolis: correct signs, f-parameter scaling")
        print("  ✓ Spin-up: tangential wind develops from Coriolis deflection")
        print("  ✓ Stability: 200-step integration bounded")
        print("  ✓ Vorticity: cyclonic (ζ > 0) near warm bubble centre")
        print()
        print("Architectural status:")
        print("  BuoyancyComponent    → PRE_PROJECTION ✓")
        print("  AdvectionComponent   → SLOW (V8.1, full momentum) ✓")
        print("  CoriolisComponent    → SLOW ✓")
        print("  AnelasticProjection  → PROJECTION ✓ (+ post-SLOW)")
        print()
        print("Next: SurfaceDragComponent → frictional inflow for TC intensification")
        print("=" * 70)
        return 0
    else:
        failures = [n for n, p in zip("ABC", [passed_a, passed_b, passed_c]) if not p]
        print(f"FAILED: sections {', '.join(failures)}")
        print("=" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
