"""
Dry Hydrostatic Balance Validation Test
=========================================

Verifies V8's full timestep preserves a stratified atmosphere at rest.
A well-implemented dynamical core should produce zero motion when
initialized with a thermodynamically closed, hydrostatically balanced
base state. Any drift indicates an implementation bug.

This is a DRY test: q = 0 throughout. Moist hydrostatic balance is a
separate test (forward-looking).

Five staged sub-tests, in order of increasing operator coverage. Each
stage isolates a different operator group so failures localize:

    Stage 1 — Buoyancy + pressure-gradient only:
        Other operators zeroed. Verifies that the buoyancy-pressure
        coupling preserves the rest state. Catches sign errors and
        staggering inconsistencies in this coupling specifically.

    Stage 2 — Stage 1 + projection:
        Adds the variable-coefficient Poisson projection to enforce the
        anelastic constraint. Verifies that projection does not
        introduce spurious motion when the constraint is already
        satisfied.

    Stage 3 — Stage 2 + Coriolis:
        Adds Coriolis forcing. Since u = v = 0 in the rest state,
        Coriolis should produce no acceleration. Catches Coriolis sign
        bugs that would produce spurious motion from rest.

    Stage 4 — Stage 3 + advection:
        Adds advection. Zero state advected by zero velocity should
        produce zero. Catches advection bugs that produce noise from
        quiescent fields (sign errors, FFT normalization slips,
        boundary aliasing).

    Stage 5 — Full dry timestep:
        All operators active. Final integration test. Passing here is
        a green light to attempt hydrostatic adjustment, not a green
        light to attempt a hurricane.

Three timing tiers:

    SHORT  (~100 timesteps)   — every commit, a few seconds each stage
    MEDIUM (~10,000 timesteps) — pre-merge gate, ~1 minute each stage
    LONG   (~600,000 timesteps) — release validation only, hours

Pass criteria (per stage, after integration):
    - max |u|, |v|, |w| < 1e-4 m/s
    - max |θ′| < 1e-3 K
    - max |∇·(ρ₀**u**)| at machine precision (< 1e-10 typical)
    - max |φ - mean(φ)| stays bounded (no growing pressure gradients)

Until V8 exists, this file's `__main__` block runs harness self-tests:
a "perfect solver" that returns rest state unchanged (must yield zero
drift), and a "broken solver" that injects small noise (must produce
detectable drift). When V8 lands, the same test framework receives V8's
solver step function with no other changes.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Callable

import numpy as np

from oracle_v8.validation import base_states as bs
from oracle_v8.validation.test_harness import compute_anelastic_residual, get_git_commit, utc_timestamp


# -------------------------------------------------------------------------
# Tier and stage definitions
# -------------------------------------------------------------------------


class TestTier(Enum):
    SHORT = ("short", 100)
    MEDIUM = ("medium", 10_000)
    LONG = ("long", 600_000)

    def __init__(self, label: str, n_steps: int):
        self.label = label
        self.n_steps = n_steps


class TestStage(Enum):
    BUOYANCY_PGRAD = "buoyancy_and_pressure_gradient"
    PLUS_PROJECTION = "plus_projection"
    PLUS_CORIOLIS = "plus_coriolis"
    PLUS_ADVECTION = "plus_advection"
    FULL_DRY_STEP = "full_dry_step"


# -------------------------------------------------------------------------
# Pass thresholds (per stage, applied after integration)
# -------------------------------------------------------------------------


@dataclass
class HydrostaticPassThresholds:
    """
    Pass criteria for the rest-state preservation test. Default values
    are conservative (well below physical relevance) but tight enough
    to catch implementation bugs.
    """
    max_velocity_m_s: float = 1e-4
    max_theta_perturbation_K: float = 1e-3
    max_anelastic_residual: float = 1e-10
    max_phi_deviation_normalized: float = 1e-6  # |φ - mean(φ)| / |p₀_max|


# -------------------------------------------------------------------------
# Test state container
# -------------------------------------------------------------------------


@dataclass
class HydrostaticState:
    """3D state for the hydrostatic balance test, broadcast from 1D base."""
    u: np.ndarray
    v: np.ndarray
    w: np.ndarray
    theta_prime: np.ndarray   # perturbation potential temperature
    phi: np.ndarray            # projection potential
    base: bs.DryBaseState      # reference base state, immutable


def initialize_rest_state(
    base: bs.DryBaseState, nx: int, ny: int
) -> HydrostaticState:
    """
    Build a 3D rest state on an (nx, ny, nz) grid by broadcasting the
    1D base state horizontally. All velocities and perturbations are
    zero by construction.
    """
    nz = len(base.z)
    shape = (nx, ny, nz)
    return HydrostaticState(
        u=np.zeros(shape),
        v=np.zeros(shape),
        w=np.zeros(shape),
        theta_prime=np.zeros(shape),
        phi=np.zeros(shape),
        base=base,
    )


# -------------------------------------------------------------------------
# Self-test solvers (used until V8 exists)
# -------------------------------------------------------------------------


def perfect_rest_solver(state: HydrostaticState, dt: float) -> HydrostaticState:
    """
    'Solver' that does nothing — returns state unchanged. Used to verify
    the harness measures zero drift on a truly stationary state.
    """
    return state


def broken_noise_solver(state: HydrostaticState, dt: float,
                          noise_amplitude: float = 1e-3) -> HydrostaticState:
    """
    'Solver' that injects small noise. Used to verify the harness
    detects drift when it occurs. The noise amplitude is well above the
    pass thresholds so the test should fail.
    """
    rng = np.random.default_rng(seed=42)
    new_state = HydrostaticState(
        u=state.u + noise_amplitude * rng.standard_normal(state.u.shape),
        v=state.v + noise_amplitude * rng.standard_normal(state.v.shape),
        w=state.w + noise_amplitude * rng.standard_normal(state.w.shape),
        theta_prime=state.theta_prime
            + noise_amplitude * rng.standard_normal(state.theta_prime.shape),
        phi=state.phi.copy(),
        base=state.base,
    )
    return new_state


# -------------------------------------------------------------------------
# Pass/fail evaluation
# -------------------------------------------------------------------------


@dataclass
class HydrostaticTestResult:
    test_name: str
    stage: str
    tier: str
    timestamp_utc: str
    git_commit: str | None
    grid_shape: tuple
    n_steps: int
    dt: float
    max_velocity: float
    max_theta_perturbation: float
    max_anelastic_residual: float
    max_phi_deviation: float
    passed: bool
    failure_reasons: list[str]
    notes: str = ""


def evaluate_rest_state(
    state: HydrostaticState,
    dx: float, dy: float, dz: float,
    thresholds: HydrostaticPassThresholds,
    stage: TestStage,
    tier: TestTier,
    n_steps: int,
    dt: float,
) -> HydrostaticTestResult:
    """
    Evaluate whether a state has remained at rest within tolerances.
    Returns a result object detailing each criterion.
    """
    failure_reasons = []

    max_velocity = float(max(
        np.max(np.abs(state.u)),
        np.max(np.abs(state.v)),
        np.max(np.abs(state.w)),
    ))
    if max_velocity > thresholds.max_velocity_m_s:
        failure_reasons.append(
            f"max |velocity| = {max_velocity:.3e} m/s exceeds threshold "
            f"{thresholds.max_velocity_m_s:.3e}"
        )

    max_theta = float(np.max(np.abs(state.theta_prime)))
    if max_theta > thresholds.max_theta_perturbation_K:
        failure_reasons.append(
            f"max |θ′| = {max_theta:.3e} K exceeds threshold "
            f"{thresholds.max_theta_perturbation_K:.3e}"
        )

    residual_diag = compute_anelastic_residual(
        state.u, state.v, state.w, state.base.rho0,
        dx, dy, dz, periodic=(True, True, False),
    )
    max_residual = residual_diag["max_abs_residual"]
    if max_residual > thresholds.max_anelastic_residual:
        failure_reasons.append(
            f"max |∇·(ρ₀u)| = {max_residual:.3e} exceeds threshold "
            f"{thresholds.max_anelastic_residual:.3e}"
        )

    # Gauge-corrected phi check
    phi_centered = state.phi - np.mean(state.phi)
    p0_max = float(np.max(state.base.p0))
    max_phi_deviation = float(np.max(np.abs(phi_centered)) / p0_max) if p0_max > 0 else 0.0
    if max_phi_deviation > thresholds.max_phi_deviation_normalized:
        failure_reasons.append(
            f"max |φ - mean(φ)| / p₀_max = {max_phi_deviation:.3e} exceeds "
            f"threshold {thresholds.max_phi_deviation_normalized:.3e}"
        )

    passed = len(failure_reasons) == 0

    return HydrostaticTestResult(
        test_name="dry_hydrostatic_balance",
        stage=stage.value,
        tier=tier.label,
        timestamp_utc=utc_timestamp(),
        git_commit=get_git_commit(),
        grid_shape=tuple(state.u.shape),
        n_steps=n_steps,
        dt=dt,
        max_velocity=max_velocity,
        max_theta_perturbation=max_theta,
        max_anelastic_residual=max_residual,
        max_phi_deviation=max_phi_deviation,
        passed=passed,
        failure_reasons=failure_reasons,
    )


# -------------------------------------------------------------------------
# Stage runners (each takes a configurable solver step function)
# -------------------------------------------------------------------------


def run_stage(
    stage: TestStage,
    tier: TestTier,
    solver_step_fn: Callable[[HydrostaticState, float], HydrostaticState],
    base: bs.DryBaseState,
    nx: int = 32,
    ny: int = 32,
    Lx: float = 2_000_000.0,
    Ly: float = 2_000_000.0,
    dt: float = 1.0,
    thresholds: HydrostaticPassThresholds | None = None,
) -> HydrostaticTestResult:
    """
    Run a single stage at a single tier. The solver_step_fn must
    advance the state by dt seconds. For pre-V8 testing, use
    perfect_rest_solver or broken_noise_solver.

    NOTE: When V8 lands, this function will need a way to tell V8 which
    operators are active for each stage. Currently the stage parameter is
    only recorded in the result; V8 will need to expose a flag-based
    operator-toggling interface, e.g. solver_step_fn(state, dt, *,
    enable_advection=False, enable_coriolis=False) for stages 1-3 to be
    meaningful. Until then, all stages run identically and the test
    primarily verifies the harness rather than V8's per-operator
    correctness.
    """
    if thresholds is None:
        thresholds = HydrostaticPassThresholds()

    state = initialize_rest_state(base, nx, ny)
    dx = Lx / nx
    dy = Ly / ny
    dz = base.z[1] - base.z[0] if len(base.z) > 1 else 1000.0

    for step in range(tier.n_steps):
        state = solver_step_fn(state, dt)

    return evaluate_rest_state(
        state, dx, dy, dz, thresholds, stage, tier, tier.n_steps, dt,
    )


def run_all_stages(
    tier: TestTier,
    solver_step_fn: Callable[[HydrostaticState, float], HydrostaticState],
    base: bs.DryBaseState,
    **kwargs,
) -> list[HydrostaticTestResult]:
    """Run all five stages at the given tier."""
    return [
        run_stage(stage, tier, solver_step_fn, base, **kwargs)
        for stage in TestStage
    ]


# -------------------------------------------------------------------------
# Self-test entry point
# -------------------------------------------------------------------------


if __name__ == "__main__":
    print("=" * 64)
    print("HARNESS SELF-TEST: dry hydrostatic balance")
    print("=" * 64)

    # Build a base state we can reuse
    nz = 32
    Lz = 20_000.0
    z = np.linspace(0, Lz, nz)
    base = bs.constant_N_dry_base_state(
        z, N=0.01, theta_surface=300.0,
        staggering=bs.GridStaggering.UNSTAGGERED_PLACEHOLDER,
    )

    print(f"\nBase state (constant-N dry):")
    print(f"  N = {base.N} 1/s, theta_surface = {base.params['theta_surface']} K")
    print(f"  z[0] to z[-1]: {base.z[0]:.0f} → {base.z[-1]:.0f} m")
    print(f"  θ₀ at top: {base.theta0[-1]:.2f} K (vs surface {base.theta0[0]:.2f})")
    print(f"  p₀ at top: {base.p0[-1]:.2f} Pa (vs surface {base.p0[0]:.2f})")
    print(f"  ρ₀ at top: {base.rho0[-1]:.4f} kg/m³ (vs surface {base.rho0[0]:.4f})")
    print(f"  load-bearing ready: {base.is_load_bearing_ready()}")
    print(f"  integration scheme: {base.integration_scheme}")

    # Diagnose discretization error of the placeholder integration
    err = bs.discrete_vs_continuous_Pi_error(z, N=0.01, theta_surface=300.0)
    print(f"\nPlaceholder Π integration error vs continuous analytical:")
    print(f"  max abs error: {err['max_abs_error']:.3e}")
    print(f"  RMS error:     {err['rms_error']:.3e}")
    print(f"  max rel error: {err['max_rel_error']:.3e}")
    print(f"  (this quantifies how far the placeholder is from continuous;")
    print(f"   replace with V8's discrete operator before declaring load-bearing)")

    # Self-test path A: perfect rest solver — should pass all stages, all tiers
    print("\n[A] Perfect rest solver (does nothing): expecting all PASS")
    print(f"{'Stage':<35} {'Tier':<8} {'max|u|':<11} {'max|θ′|':<11} {'residual':<11} {'pass'}")
    for tier in [TestTier.SHORT, TestTier.MEDIUM]:
        for stage in TestStage:
            result = run_stage(stage, tier, perfect_rest_solver, base)
            tag = "✓" if result.passed else "✗"
            print(f"  {stage.value:<33} {tier.label:<8} "
                  f"{result.max_velocity:<11.2e} "
                  f"{result.max_theta_perturbation:<11.2e} "
                  f"{result.max_anelastic_residual:<11.2e} "
                  f"{tag}")

    # Self-test path B: broken noise solver — should fail all stages
    print("\n[B] Broken noise solver: expecting all FAIL")
    print(f"{'Stage':<35} {'Tier':<8} {'max|u|':<11} {'max|θ′|':<11} {'pass'}")
    for stage in TestStage:
        result = run_stage(stage, TestTier.SHORT, broken_noise_solver, base)
        tag = "✓" if result.passed else "✗"
        print(f"  {stage.value:<33} {TestTier.SHORT.label:<8} "
              f"{result.max_velocity:<11.2e} "
              f"{result.max_theta_perturbation:<11.2e} "
              f"{tag}")

    print("\n" + "=" * 64)
    print("Self-test interpretation:")
    print("  Path A all ✓: harness correctly identifies preserved rest states")
    print("  Path B all ✗: harness correctly identifies broken solvers")
    print("  Discretization error finite & small: placeholder integration is")
    print("    consistent with continuous theory at this resolution")
    print()
    print("This file is READY to receive V8's solver step function. Replace")
    print("perfect_rest_solver with V8's actual stepper to begin testing.")
    print("=" * 64)
