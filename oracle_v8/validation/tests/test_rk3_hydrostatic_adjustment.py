"""
RK3 Hydrostatic Adjustment Integration Test
=============================================
V8.0 — post ensemble review (Claude + Gemini + Five)

ARCHITECTURAL NOTE (see rk3_diagnostic_brief.md for full derivation):
----------------------------------------------------------------------
PressureGradientComponent is NOT included in the V8.0 integrator config.

In the current full-potential projection architecture, AnelasticProjection
solves ∇·(ρ̄∇φ) = ∇·(ρ̄u*) and applies u_final = u* − ∇φ. The resulting
φ is a timestep-scaled pressure impulse (dimensions of kinematic pressure × Δt),
NOT a physical pressure acceleration. Adding PressureGradientComponent as an
explicit tendency would apply the pressure gradient twice per step:
  1. Explicit: u_prov += Δt × (−∇φ_old)    [PGC tendency]
  2. Implicit: u_final  = u_prov − ∇φ_new   [projection correction]

This double-counting gives the recurrence aₙ₊₁ = Δt(1 − aₙ), with
amplification factor −Δt per step → instability for Δt > 1 s.

PressureGradientComponent is valid as an isolated gradient operator
(tests A, B, C in test_pressure_gradient_component.py all pass) and is
retained in the codebase for:
  (a) V8 paper architectural documentation
  (b) Future V8.x incremental-pressure formulation, where the Poisson
      solver will solve for δφ = φⁿ⁺¹ − φⁿ and the explicit PG uses
      the accumulated φⁿ from the previous step.

PRE-ADVECTION PHYSICS:
----------------------
With fixed θ′ (AdvectionComponent not yet implemented), the buoyancy source
is constant in space and time. The correct pre-advection behavior is:
  - φ stays approximately constant each step (the per-step pressure impulse
    is set by the fixed buoyancy divergence, which does not grow)
  - φ/Δt ≈ −64.843 Pa/[Δt_units] is the Δt-invariant diagnostic
  - Horizontal circulation grows linearly (projection-induced inflow)
  - The surface pressure deficit does NOT deepen toward the −123.57 Pa
    equilibrium without θ′ redistribution — that is a post-advection test

Sections
--------
[A] First-step consistency: φ = −64.843 Pa (pure-buoyancy baseline)
[B] 100-step trajectory: stable, φ approximately constant, u grows
[C] Physics checks: φ constant, φ/Δt invariant, circulation grows
[D] Discriminating test (Five/Gemini recommended):
    Without PGC → φ/Δt invariant across Δt=1,2,5 s
    With PGC    → alternating sign recurrence diagnostic
"""

from __future__ import annotations

import sys
from typing import List

import numpy as np

from oracle_v8.solver import (
    LH82AnelasticEquationSet,
    BuoyancyComponent,
    PressureGradientComponent,
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
    expected_pressure_deficit_at_center,
)


# --------------------------------------------------------------------------
# Shared grid and base state
# --------------------------------------------------------------------------

Lx = Ly = 100_000.0
Lz = 10_000.0
nx = ny = 64
nz = 32
dz = Lz / nz
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


def make_config_v8() -> OperatorConfig:
    """
    Correct V8.0 integrator config: buoyancy + projection only.

    PressureGradientComponent is intentionally excluded.
    See module docstring and rk3_diagnostic_brief.md for the full
    explanation.  Short version: the projection already applies the
    complete pressure impulse; adding PGC would double-count it and
    produce recurrence instability for Δt > 1 s.
    """
    return OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        buoyancy=BuoyancyComponent(),
        # pressure_gradient=None  ← excluded; see module docstring
        projection=AnelasticProjection(
            nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
        ),
    )


def make_config_with_pgc() -> OperatorConfig:
    """
    Config WITH PressureGradientComponent — for discriminating test only.
    Expected to exhibit aₙ = Δt(1 − aₙ₋₁) instability per the brief.
    """
    return OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        buoyancy=BuoyancyComponent(),
        pressure_gradient=PressureGradientComponent(
            nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
        ),
        projection=AnelasticProjection(
            nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
        ),
    )


def make_initial_state(theta_prime: np.ndarray) -> State:
    return State(
        u=np.zeros((nx, ny, nz)),
        v=np.zeros((nx, ny, nz)),
        w=np.zeros((nx, ny, nz + 1)),
        theta_prime=theta_prime,
        projection_potential=np.zeros((nx, ny, nz)),
        t=0.0,
    )


# --------------------------------------------------------------------------
# Section A: first-step consistency
# --------------------------------------------------------------------------

def section_a(theta_prime: np.ndarray) -> tuple[bool, float]:
    print("\n[A] First-step consistency (dt=1s, buoyancy+projection only)")
    integrator = RK3Integrator(config=make_config_v8(), base=BASE)
    state = make_initial_state(theta_prime)
    state, diag = integrator.step(state, dt=1.0, step_number=0)

    phi_c = state.projection_potential - np.mean(state.projection_potential)
    phi_min = float(np.min(phi_c[:, :, 0]))
    phi_over_dt = phi_min / 1.0

    print(f"  φ_min     = {phi_min:.4f} Pa")
    print(f"  φ/Δt      = {phi_over_dt:.4f} (invariant diagnostic)")
    print(f"  max |u|   = {diag.max_u:.4e} m/s")
    print(f"  max |w|   = {diag.max_w:.4e} m/s")
    print(f"  compat    = {diag.stage_compat_residuals}")

    baseline = -64.843
    ok = abs(phi_min - baseline) < 0.01
    print(f"  {'✓' if ok else '✗'} matches pure-buoyancy baseline "
          f"({baseline} Pa)")
    return ok, phi_min


# --------------------------------------------------------------------------
# Section B: 100-step trajectory
# --------------------------------------------------------------------------

def section_b(
    theta_prime: np.ndarray,
    n_steps: int = 100,
    dt: float = 5.0,
) -> tuple[bool, List[StepDiagnostics]]:
    print(f"\n[B] Trajectory: {n_steps} steps × dt={dt}s = "
          f"{n_steps*dt:.0f}s "
          f"(buoyancy+projection only, PGC excluded)")

    integrator = RK3Integrator(config=make_config_v8(), base=BASE)
    state = make_initial_state(theta_prime)
    history: List[StepDiagnostics] = []

    print(f"\n  {'Step':>5}  {'t (s)':>7}  {'φ_min (Pa)':>12}  "
          f"{'φ/Δt':>10}  {'max|u| (m/s)':>13}")
    print(f"  {'-'*5}  {'-'*7}  {'-'*12}  {'-'*10}  {'-'*13}")

    any_nan = False
    for n in range(n_steps):
        state, diag = integrator.step(state, dt=dt, step_number=n)
        history.append(diag)

        if (np.any(np.isnan(state.u))
                or np.any(np.isnan(state.w))
                or np.any(np.isnan(state.projection_potential))):
            print(f"\n  ✗ NaN at step {n} — instability")
            any_nan = True
            break

        if n == 0 or (n + 1) % 10 == 0:
            phi_over_dt = diag.surface_phi_min / dt
            print(f"  {n+1:>5}  {state.t:>7.1f}  "
                  f"{diag.surface_phi_min:>12.3f}  "
                  f"{phi_over_dt:>10.3f}  "
                  f"{diag.max_u:>13.4e}")

    if not any_nan:
        print(f"\n  ✓ No NaN in {len(history)} steps")

    max_u_ever = max(d.max_u for d in history) if history else 0
    max_w_ever = max(d.max_w for d in history) if history else 0
    stable = (not any_nan) and max_u_ever < 10.0 and max_w_ever < 50.0
    if stable:
        print(f"  ✓ Stable: max|u|={max_u_ever:.4e}, max|w|={max_w_ever:.4e} m/s")
    else:
        print(f"  ✗ Unstable: max|u|={max_u_ever:.4e}, max|w|={max_w_ever:.4e}")

    return stable, history


# --------------------------------------------------------------------------
# Section C: physics checks
# --------------------------------------------------------------------------

def section_c(
    history: List[StepDiagnostics],
    dt: float,
) -> bool:
    print(f"\n[C] Physics checks (pre-advection behavior, dt={dt}s)")

    if not history:
        print("  ✗ Empty history")
        return False

    phi_step1 = history[0].surface_phi_min
    phi_final = history[-1].surface_phi_min

    # C.1: φ approximately constant (fixed buoyancy → constant pressure impulse)
    # Theory: without PGC, ∇·(ρ̄u_prov) = ∇·(ρ̄u^n) + dt·D_b = 0 + dt·D_b = const
    # → φ^n = dt·Φ_b at every step.
    max_drift = max(abs(d.surface_phi_min - phi_step1) for d in history)
    drift_frac = max_drift / abs(phi_step1) if phi_step1 != 0 else float('inf')
    drift_ok = drift_frac < 0.05   # within 5%

    print(f"\n  [C.1] φ approximately constant (fixed-θ′ theory predicts exactly constant):")
    print(f"    φ at step 1:   {phi_step1:.3f} Pa")
    print(f"    φ at step {len(history)}: {phi_final:.3f} Pa")
    print(f"    max drift:     {max_drift:.3f} Pa  ({100*drift_frac:.2f}% of step-1)")
    if drift_ok:
        print(f"    ✓ φ constant to within 5% — correct pre-advection physics")
    else:
        print(f"    ✗ φ drifted more than 5% — unexpected")

    # C.2: φ/Δt ≈ −64.843 (Δt-invariant diagnostic; per Five and Gemini)
    # "Raw φ scales with Δt; φ/Δt should be timestep-independent."
    phi_over_dt_values = [d.surface_phi_min / dt for d in history]
    mean_phi_over_dt = float(np.mean(phi_over_dt_values))
    std_phi_over_dt = float(np.std(phi_over_dt_values))
    target = -64.843
    invariant_ok = abs(mean_phi_over_dt - target) < 1.0

    print(f"\n  [C.2] φ/Δt invariant diagnostic:")
    print(f"    mean φ/Δt = {mean_phi_over_dt:.3f}  (target ≈ {target})")
    print(f"    std  φ/Δt = {std_phi_over_dt:.4f}  (should be small)")
    if invariant_ok:
        print(f"    ✓ φ/Δt ≈ {target} — Δt-independent diagnostic confirmed")
    else:
        print(f"    ✗ φ/Δt deviates from {target}")

    # C.3: Horizontal circulation grows (projection-induced inflow developing)
    max_u_step1 = history[0].max_u
    max_u_final = history[-1].max_u
    circulation_grows = max_u_final > 10 * max_u_step1
    print(f"\n  [C.3] Secondary circulation (projection-induced inflow):")
    print(f"    max|u| step 1:   {max_u_step1:.4e} m/s")
    print(f"    max|u| step {len(history)}: {max_u_final:.4e} m/s")
    print(f"    growth factor: {max_u_final/max_u_step1:.1f}×")
    if circulation_grows:
        print(f"    ✓ Horizontal inflow developing (>10× growth)")
    else:
        print(f"    ✗ Horizontal inflow not growing as expected")

    # C.4: Projection residuals at machine precision throughout
    max_compat = max(r for d in history for r in d.stage_compat_residuals)
    max_disc_op = max(r for d in history for r in d.stage_disc_op_residuals)
    residuals_ok = max_compat < 1e-6 and max_disc_op < 1e-6
    print(f"\n  [C.4] Projection residuals (all {len(history)*3} sub-stages):")
    print(f"    max compat   = {max_compat:.3e}  "
          f"({'✓' if max_compat < 1e-6 else '✗'})")
    print(f"    max disc_op  = {max_disc_op:.3e}  "
          f"({'✓' if max_disc_op < 1e-6 else '✗'})")

    return drift_ok and invariant_ok and circulation_grows and residuals_ok


# --------------------------------------------------------------------------
# Section D: discriminating test (Five + Gemini)
# --------------------------------------------------------------------------

def section_d(theta_prime: np.ndarray) -> bool:
    """
    Discriminating test demonstrating that:
      (a) Without PGC: φ/Δt is Δt-invariant (architectural correctness)
      (b) With    PGC: alternating-sign recurrence appears immediately
                       (confirms the double-counting diagnosis)

    Proving (a) and (b) together makes the architectural choice
    empirically verifiable, not just theoretically motivated.
    """
    print(f"\n[D] Discriminating dt-scaling test (Five + Gemini recommendation)")

    # Part D.1: φ/Δt invariant WITHOUT PGC
    print(f"\n  [D.1] φ/Δt invariance across Δt = 1, 2, 5 s (no PGC):")
    print(f"  {'Δt':>5}  {'φ_min (Pa)':>12}  {'φ/Δt':>10}  {'Δt-invariant?':>14}")

    phi_over_dt_vals = []
    for dt_test in [1.0, 2.0, 5.0]:
        integrator = RK3Integrator(config=make_config_v8(), base=BASE)
        state = make_initial_state(theta_prime)
        state, diag = integrator.step(state, dt=dt_test, step_number=0)
        phi_c = state.projection_potential - np.mean(state.projection_potential)
        phi_min = float(np.min(phi_c[:, :, 0]))
        phi_dt = phi_min / dt_test
        phi_over_dt_vals.append(phi_dt)
        target = -64.843
        ok = abs(phi_dt - target) < 0.5
        print(f"  {dt_test:>5.1f}  {phi_min:>12.3f}  {phi_dt:>10.3f}  "
              f"{'✓' if ok else '✗'} (target {target})")

    invariant_ok = all(
        abs(v - phi_over_dt_vals[0]) < 1.0 for v in phi_over_dt_vals
    )
    if invariant_ok:
        print(f"  ✓ φ/Δt invariant — projection is the correct pressure operator")
    else:
        print(f"  ✗ φ/Δt varies — unexpected")

    # Part D.2: alternating-sign recurrence WITH PGC (dt=2s)
    # Theory: a₁ = 2, a₂ = 2(1−2) = −2, a₃ = 2(1−(−2)) = 6, ...
    # The sign alternation and growing magnitude is the diagnostic signature.
    print(f"\n  [D.2] Recurrence signature WITH PGC (Δt=2s, 3 steps):")
    print(f"  Theory: a₁=2, a₂=−2, a₃=6 → φ alternates sign, grows each step")
    print(f"  {'Step':>5}  {'φ_min (Pa)':>12}  {'sign flipped?':>14}")

    integrator_pgc = RK3Integrator(config=make_config_with_pgc(), base=BASE)
    state_pgc = make_initial_state(theta_prime)
    prev_phi = None
    sign_flips = 0
    pgc_phis = []

    for n in range(3):
        state_pgc, diag_pgc = integrator_pgc.step(
            state_pgc, dt=2.0, step_number=n
        )
        phi_c = (state_pgc.projection_potential
                 - np.mean(state_pgc.projection_potential))
        phi_min = float(np.min(phi_c[:, :, 0]))
        pgc_phis.append(phi_min)
        flipped = (prev_phi is not None
                   and np.sign(phi_min) != np.sign(prev_phi))
        if flipped:
            sign_flips += 1
        print(f"  {n+1:>5}  {phi_min:>12.3f}  {'✓ YES' if flipped else 'no'}")
        prev_phi = phi_min

    recurrence_confirmed = sign_flips >= 1 and abs(pgc_phis[-1]) > abs(pgc_phis[0])
    if recurrence_confirmed:
        print(f"  ✓ Recurrence confirmed: alternating sign, growing magnitude")
        print(f"    This validates the double-counting diagnosis.")
    else:
        print(f"  ✗ Recurrence not evident — review diagnosis")

    return invariant_ok and recurrence_confirmed


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main() -> int:
    print("=" * 70)
    print("RK3 HYDROSTATIC ADJUSTMENT: TIME-INTEGRATION TEST")
    print("(V8.0 — post-ensemble review, PGC excluded from integrator)")
    print("=" * 70)

    bubble = WarmBubbleParams()
    x = np.linspace(0, Lx, nx, endpoint=False)
    y = np.linspace(0, Ly, ny, endpoint=False)
    X, Y, Z = np.meshgrid(x, y, z_centers, indexing="ij")
    theta_prime = warm_bubble_theta_perturbation(X, Y, Z, Lx, Ly, bubble)
    equilibrium_deficit = expected_pressure_deficit_at_center(BASE, bubble)

    print(f"\nSetup: 64×64×32 grid, B&F warm bubble (2 K, 10 km radius, 2 km alt)")
    print(f"  max θ′ = {theta_prime.max():.4f} K")
    print(f"  Analytical equilibrium deficit: {equilibrium_deficit:.2f} Pa")
    print(f"  (Deepening toward this requires AdvectionComponent — deferred)")

    passed_a, _  = section_a(theta_prime)
    passed_b, history = section_b(theta_prime, n_steps=100, dt=5.0)
    passed_c = section_c(history, dt=5.0)
    passed_d = section_d(theta_prime)

    print("\n" + "=" * 70)
    if passed_a and passed_b and passed_c and passed_d:
        print("PASSED: RK3 integrator verified (V8.0 architecture)")
        print()
        print("Confirmed:")
        print("  ✓ First step: exact −64.843 Pa baseline (pure buoyancy)")
        print("  ✓ 100-step trajectory: stable, φ constant, u grows linearly")
        print("  ✓ φ/Δt invariant diagnostic ≈ −64.843 (Δt-independent)")
        print("  ✓ Secondary circulation develops (projection-induced)")
        print("  ✓ Discriminating test: φ/Δt invariant without PGC,")
        print("    alternating-sign recurrence with PGC (diagnosis confirmed)")
        print()
        print("Architectural status:")
        print("  BuoyancyComponent    → PRE_PROJECTION  ✓ active")
        print("  PressureGradient     → excluded V8.0   (reserved for V8.x")
        print("                          incremental-pressure formulation)")
        print("  AnelasticProjection  → PROJECTION      ✓ active (= PG force)")
        print()
        print("Next: AdvectionComponent — enables θ′ redistribution, which")
        print("drives the deficit deepening toward −123.57 Pa equilibrium.")
        print("=" * 70)
        return 0
    else:
        failures = []
        if not passed_a: failures.append("A")
        if not passed_b: failures.append("B")
        if not passed_c: failures.append("C")
        if not passed_d: failures.append("D")
        print(f"FAILED: sections {', '.join(failures)}")
        print("=" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
