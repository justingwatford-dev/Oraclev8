"""
Hydrostatic Adjustment Validation Test
========================================

The V7-to-V8 distinguishing test. V7's Boussinesq core was structurally
incapable of passing this; V8's anelastic core must pass it for V8 to
have any meaning.

Setup: a stably stratified dry atmosphere at rest is initialized, then
a warm potential-temperature perturbation (a "bubble") is added. The
atmosphere responds: the warm bubble is positively buoyant, the column
beneath it loses mass via the anelastic continuity constraint, and a
pressure deficit develops collocated with the warm anomaly. The system
radiates gravity waves at the Brunt-Väisälä frequency and reaches a
new perturbed hydrostatic balance after ~5-10 buoyancy periods.

V7 BEHAVIOR (documented empirically): max p_range across 15,000 frames
of Hugo, Katrina, Ivan was 0.0031 (dimensional ~1.9 Pa). The Boussinesq
incompressibility constraint forbids hydrostatic mass evacuation; the
pressure field is purely diagnostic and never deepens.

V8 EXPECTATION: pressure deficit of order tens to hundreds of Pa for
a Bryan & Fritsch test bubble (2 K, 10 km radius, 2 km altitude). This
is a 100-1000x signal-to-noise advantage over V7's noise floor.

Reference setup: Bryan & Fritsch (2002), "A Benchmark Simulation for
Moist Nonhydrostatic Numerical Models," Mon. Wea. Rev. 130, 2917-2928,
§3 (the dry warm bubble case is the cited reference standard).

This test does NOT have a closed-form time-dependent analytical solution.
The pass criteria are:

    1. Pressure deficit develops (not zero, not positive).
    2. Magnitude is within a factor of 2 of the analytical expected value
       computed from hydrostatic theory.
    3. Deficit is collocated with the warm anomaly (within ~1 bubble
       radius horizontally, within ~2 bubble radii vertically — gravity
       wave radiation tilts the response).
    4. Anelastic constraint preserved at machine precision throughout.

Tier coverage: SHORT (a few buoyancy periods, ~600 steps) and MEDIUM
(adjustment plus relaxation, ~6000 steps). LONG is not meaningful here
— hydrostatic adjustment completes in minutes of simulated time.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Callable

import numpy as np

from oracle_v8.validation import base_states as bs
from oracle_v8.validation.test_harness import (
    compute_anelastic_residual,
    get_git_commit,
    utc_timestamp,
)
from oracle_v8.validation.tests.test_dry_hydrostatic_balance import (
    HydrostaticState,
    initialize_rest_state,
    TestTier as RestTier,
)


# -------------------------------------------------------------------------
# Bryan & Fritsch (2002) warm bubble shape
# -------------------------------------------------------------------------
#
# B&F prescribe a cosine-cubed bubble that is smooth (continuous to second
# derivatives at the boundary) and finite-support (zero outside the
# specified radius). The functional form is:
#
#     θ′(x, y, z) = θ_max · cos²(π L / 2)  if L ≤ 1
#                  0                        if L > 1
#
# where L is the normalized distance from the bubble center:
#
#     L = sqrt( ((x-xc)/xr)² + ((y-yc)/yr)² + ((z-zc)/zr)² )
#
# For our test, default parameters follow B&F's dry warm bubble case:
#
#     θ_max = 2.0 K            (perturbation amplitude)
#     xr = yr = 10 km          (horizontal bubble radius)
#     zr = 2 km                (vertical bubble radius — vertically smaller
#                                so the bubble is realistic for adjustment;
#                                a tall narrow bubble adjusts cleanly while
#                                a spherical bubble of equal radii would
#                                be vertically smeared by stratification)
#     (xc, yc, zc) = domain center horizontally, 2 km altitude
#
# Reference: Bryan & Fritsch (2002), Mon. Wea. Rev. 130, 2917-2928, §3a.
# -------------------------------------------------------------------------


@dataclass
class WarmBubbleParams:
    """Bryan & Fritsch warm bubble parameters."""
    theta_max_K: float = 2.0
    xr_m: float = 10_000.0
    yr_m: float = 10_000.0
    zr_m: float = 2_000.0
    xc_m: float | None = None  # None → domain center
    yc_m: float | None = None  # None → domain center
    zc_m: float = 2_000.0


def warm_bubble_theta_perturbation(
    X: np.ndarray, Y: np.ndarray, Z: np.ndarray,
    Lx: float, Ly: float, params: WarmBubbleParams,
) -> np.ndarray:
    """
    Compute the Bryan & Fritsch warm bubble θ′ field on the (X, Y, Z) grid.
    """
    xc = params.xc_m if params.xc_m is not None else Lx / 2
    yc = params.yc_m if params.yc_m is not None else Ly / 2
    zc = params.zc_m

    L_sq = (
        ((X - xc) / params.xr_m) ** 2
        + ((Y - yc) / params.yr_m) ** 2
        + ((Z - zc) / params.zr_m) ** 2
    )
    L = np.sqrt(L_sq)

    theta_prime = np.where(
        L <= 1.0,
        params.theta_max_K * np.cos(np.pi * L / 2.0) ** 2,
        0.0,
    )
    return theta_prime


# -------------------------------------------------------------------------
# Analytical expected pressure deficit
# -------------------------------------------------------------------------
#
# Derivation: for a warm potential-temperature perturbation θ′ embedded
# in a hydrostatically balanced base state, the perturbation buoyancy
# is b ≈ g · θ′/θ₀. After hydrostatic adjustment, the column-integrated
# pressure perturbation beneath the bubble is approximately:
#
#     Δp(z) ≈ -∫_z^∞ ρ₀(z') · g · (θ′(z')/θ₀(z')) dz'
#
# This is the *hydrostatic* pressure response — what V8 should produce
# in the limit of complete adjustment. The signal is negative (a deficit)
# and order ρ₀ · g · θ′/θ₀ · zr ≈ 1.0 · 9.81 · 2/300 · 2000 ≈ 130 Pa for
# the default B&F bubble. This is the order-of-magnitude target.
#
# CAVEATS to be honest about:
#   - This is the equilibrium answer, not the time-dependent answer.
#   - The actual V8 response will include gravity wave radiation, so
#     instantaneous deficit will oscillate around the equilibrium value
#     before settling. We require V8 to be within a factor of 2 of this
#     value at any time after the first buoyancy period (~100 s for N=0.01).
#   - At adjustment timescales the integral is approximated; for long-time
#     equilibrium see the steady-state Boussinesq calculation in Durran
#     (2010) §2.4 which is the reference for the order-of-magnitude check.
# -------------------------------------------------------------------------


def expected_pressure_deficit_at_center(
    base: bs.DryBaseState,
    params: WarmBubbleParams,
) -> float:
    """
    Compute the order-of-magnitude expected pressure deficit at the
    bubble center column, after hydrostatic adjustment.

    Approximation: integrate ρ₀ · g · (θ′/θ₀) over the bubble's vertical
    extent at the centerline (where θ′ is maximum).

    Returns the deficit as a NEGATIVE number (pressure perturbation is
    negative in a warm column).
    """
    z = base.z
    rho0 = base.rho0
    theta0 = base.theta0

    # θ′ along the centerline (x=xc, y=yc): bubble shape reduces to
    # θ_max · cos²(π/2 · |z-zc|/zr) for |z-zc| ≤ zr
    z_offset = (z - params.zc_m) / params.zr_m
    theta_prime_centerline = np.where(
        np.abs(z_offset) <= 1.0,
        params.theta_max_K * np.cos(np.pi * z_offset / 2.0) ** 2,
        0.0,
    )

    # Hydrostatic pressure perturbation: dp/dz = -ρ₀ g (θ′/θ₀)
    # Integrate from top down to find p′(z=0) = -∫₀^∞ ρ₀ g θ′/θ₀ dz
    integrand = rho0 * bs.GRAVITY * theta_prime_centerline / theta0
    # Trapezoidal integration over the vertical extent
    nz = len(z)
    dz_vec = np.diff(z)
    p_deficit = 0.0
    for k in range(nz - 1):
        p_deficit += -0.5 * (integrand[k] + integrand[k + 1]) * dz_vec[k]

    return float(p_deficit)


# -------------------------------------------------------------------------
# Test setup and pass evaluation
# -------------------------------------------------------------------------


@dataclass
class AdjustmentPassThresholds:
    """Pass criteria for the hydrostatic adjustment test."""
    deficit_must_be_negative: bool = True
    deficit_magnitude_factor_lower: float = 0.5  # at least 0.5x expected
    deficit_magnitude_factor_upper: float = 2.0  # at most 2x expected
    max_horizontal_offset_radii: float = 1.0     # deficit center within 1 xr
    max_vertical_offset_radii: float = 2.0       # deficit center within 2 zr
    max_anelastic_residual: float = 1e-8         # constraint preservation
                                                  # slightly looser than rest
                                                  # state because real motion
                                                  # is happening


@dataclass
class AdjustmentTestResult:
    test_name: str
    tier: str
    timestamp_utc: str
    git_commit: str | None
    grid_shape: tuple
    n_steps: int
    dt: float
    expected_deficit_Pa: float
    measured_deficit_Pa: float
    deficit_ratio: float                # measured / expected
    horizontal_offset_radii: float
    vertical_offset_radii: float
    max_anelastic_residual: float
    passed: bool
    failure_reasons: list[str]
    notes: str = ""


def initialize_adjustment_state(
    base: bs.DryBaseState,
    nx: int, ny: int,
    Lx: float, Ly: float,
    bubble_params: WarmBubbleParams,
) -> HydrostaticState:
    """
    Build the initial state: rest state + warm bubble θ′ perturbation.
    All velocities and φ remain zero.
    """
    state = initialize_rest_state(base, nx, ny)

    nz = len(base.z)
    x = np.linspace(0, Lx, nx, endpoint=False)
    y = np.linspace(0, Ly, ny, endpoint=False)
    z = base.z
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    state.theta_prime = warm_bubble_theta_perturbation(
        X, Y, Z, Lx, Ly, bubble_params,
    )
    return state


def evaluate_adjustment(
    state: HydrostaticState,
    base: bs.DryBaseState,
    bubble_params: WarmBubbleParams,
    dx: float, dy: float, dz: float,
    Lx: float, Ly: float,
    thresholds: AdjustmentPassThresholds,
    tier: RestTier,
    n_steps: int,
    dt: float,
) -> AdjustmentTestResult:
    """
    Evaluate whether hydrostatic adjustment produced the expected
    pressure deficit.

    The pressure deficit from a warm bubble is column-integrated: most
    negative at the surface beneath the bubble (where the entire warm
    column contributes), monotonically less negative as you ascend
    through and above the bubble. So we sample φ at the surface beneath
    the bubble center, which is where the deficit is deepest. This
    matches the expected_pressure_deficit_at_center calculation, which
    integrates the full bubble vertical extent.
    """
    failure_reasons = []

    expected_deficit = expected_pressure_deficit_at_center(base, bubble_params)

    # Sample at the surface beneath the bubble center, where the
    # column-integrated deficit is deepest.
    nx, ny, nz = state.phi.shape
    xc_idx = int(round((bubble_params.xc_m or Lx / 2) / dx))
    yc_idx = int(round((bubble_params.yc_m or Ly / 2) / dy))
    surface_idx = 0

    # Sample φ at bubble center column, surface (gauge-corrected by
    # subtracting domain mean to remove additive constants from the
    # nullspace of the elliptic operator)
    phi_centered = state.phi - np.mean(state.phi)
    measured_deficit = float(phi_centered[xc_idx, yc_idx, surface_idx])

    # Find the surface-level location of the deepest deficit (where
    # the column response is strongest). For equilibrium hydrostatic
    # adjustment, this should be directly beneath the bubble center.
    surface_phi = phi_centered[:, :, surface_idx]
    if np.all(surface_phi == 0):
        # Degenerate case (V7-style no-response): no minimum to find,
        # report bubble center as the "minimum" so offset is zero
        # (the failure will be flagged via the deficit-not-negative check)
        actual_xc = bubble_params.xc_m or Lx / 2
        actual_yc = bubble_params.yc_m or Ly / 2
    else:
        min_idx = np.unravel_index(np.argmin(surface_phi), surface_phi.shape)
        actual_xc = min_idx[0] * dx
        actual_yc = min_idx[1] * dy

    bubble_xc = bubble_params.xc_m or Lx / 2
    bubble_yc = bubble_params.yc_m or Ly / 2

    horizontal_offset = (
        np.sqrt((actual_xc - bubble_xc) ** 2 + (actual_yc - bubble_yc) ** 2)
        / bubble_params.xr_m
    )
    # Vertical structure check: verify the deficit decreases in magnitude
    # with height (phi at surface < phi at bubble top, since phi is
    # negative below the bubble and ~0 above it). This replaces the
    # "vertical offset" check, which doesn't make sense for hydrostatic
    # equilibrium where the deficit is always at the surface.
    bubble_top_idx = int(round((bubble_params.zc_m + bubble_params.zr_m) / dz))
    bubble_top_idx = min(bubble_top_idx, nz - 1)
    phi_at_surface = phi_centered[xc_idx, yc_idx, surface_idx]
    phi_above_bubble = phi_centered[xc_idx, yc_idx, bubble_top_idx]
    # Set vertical_offset to a structural diagnostic: 0 if the structure
    # is correct (more negative at surface than above bubble), positive
    # otherwise. We track this as informational; the hard pass criterion
    # is the magnitude check below.
    if phi_at_surface < phi_above_bubble:
        vertical_structure_correct = True
        vertical_offset = 0.0
    else:
        vertical_structure_correct = False
        vertical_offset = float("inf")  # signals broken vertical structure

    # Pass checks
    if thresholds.deficit_must_be_negative and measured_deficit >= 0:
        failure_reasons.append(
            f"pressure deficit is not negative: measured = {measured_deficit:.3e} Pa "
            f"(expected ~{expected_deficit:.3e} Pa). This is the V7 failure mode: "
            f"no hydrostatic mass evacuation occurred."
        )

    if expected_deficit != 0:
        deficit_ratio = measured_deficit / expected_deficit
    else:
        deficit_ratio = float("nan")

    if (
        not np.isnan(deficit_ratio)
        and (deficit_ratio < thresholds.deficit_magnitude_factor_lower
             or deficit_ratio > thresholds.deficit_magnitude_factor_upper)
    ):
        failure_reasons.append(
            f"pressure deficit magnitude out of range: ratio "
            f"measured/expected = {deficit_ratio:.2f}, "
            f"acceptable range "
            f"[{thresholds.deficit_magnitude_factor_lower}, "
            f"{thresholds.deficit_magnitude_factor_upper}]"
        )

    if horizontal_offset > thresholds.max_horizontal_offset_radii:
        failure_reasons.append(
            f"deficit center is {horizontal_offset:.2f} bubble radii "
            f"horizontally from bubble center; threshold is "
            f"{thresholds.max_horizontal_offset_radii}"
        )

    if vertical_offset > thresholds.max_vertical_offset_radii:
        failure_reasons.append(
            f"vertical structure incorrect: φ at surface beneath bubble "
            f"({phi_at_surface:.3e} Pa) is not more negative than φ above "
            f"the bubble ({phi_above_bubble:.3e} Pa); column response is "
            f"inverted or absent"
        )

    residual_diag = compute_anelastic_residual(
        state.u, state.v, state.w, base.rho0,
        dx, dy, dz, periodic=(True, True, False),
    )
    max_residual = residual_diag["max_abs_residual"]
    if max_residual > thresholds.max_anelastic_residual:
        failure_reasons.append(
            f"max |∇·(ρ₀u)| = {max_residual:.3e} exceeds threshold "
            f"{thresholds.max_anelastic_residual:.3e}"
        )

    passed = len(failure_reasons) == 0

    return AdjustmentTestResult(
        test_name="dry_hydrostatic_adjustment",
        tier=tier.label,
        timestamp_utc=utc_timestamp(),
        git_commit=get_git_commit(),
        grid_shape=tuple(state.phi.shape),
        n_steps=n_steps,
        dt=dt,
        expected_deficit_Pa=expected_deficit,
        measured_deficit_Pa=measured_deficit,
        deficit_ratio=deficit_ratio,
        horizontal_offset_radii=horizontal_offset,
        vertical_offset_radii=vertical_offset,
        max_anelastic_residual=max_residual,
        passed=passed,
        failure_reasons=failure_reasons,
    )


# -------------------------------------------------------------------------
# Self-test solvers (used until V8 exists)
# -------------------------------------------------------------------------


def perfect_adjustment_solver(
    state: HydrostaticState, dt: float,
    *,
    bubble_params: WarmBubbleParams,
    Lx: float, Ly: float,
) -> HydrostaticState:
    """
    'Solver' that produces the expected hydrostatic equilibrium response:
    sets φ to the analytical column-integrated pressure deficit at every
    point. Used to verify the harness recognizes a correct V8 response.

    This is NOT a real solver — it's a known-answer injection for
    self-testing the harness.
    """
    base = state.base
    nx, ny, nz = state.phi.shape
    z = base.z
    rho0 = base.rho0
    theta0 = base.theta0
    theta_prime = state.theta_prime

    # Compute column-integrated hydrostatic perturbation pressure at every
    # (x, y) column from the local θ′ profile:
    #     p′(x,y,z) = -∫_z^∞ ρ₀ g θ′/θ₀ dz'
    # Integrating top-down for numerical stability.
    integrand = rho0[None, None, :] * bs.GRAVITY * theta_prime / theta0[None, None, :]

    # Trapezoidal integration from top of domain down
    p_perturbation = np.zeros_like(theta_prime)
    for k in range(nz - 2, -1, -1):
        dz_local = z[k + 1] - z[k]
        p_perturbation[:, :, k] = (
            p_perturbation[:, :, k + 1]
            - 0.5 * (integrand[:, :, k] + integrand[:, :, k + 1]) * dz_local
        )

    new_state = HydrostaticState(
        u=state.u.copy(),
        v=state.v.copy(),
        w=state.w.copy(),
        theta_prime=state.theta_prime.copy(),
        phi=p_perturbation,
        base=state.base,
    )
    return new_state


def v7_style_no_deficit_solver(
    state: HydrostaticState, dt: float, **kwargs,
) -> HydrostaticState:
    """
    'Solver' that produces no pressure deficit (mimics V7's Boussinesq
    behavior). φ stays at zero throughout. Used to verify the harness
    correctly fails this case — i.e. detects the V7 failure mode.
    """
    return state


# -------------------------------------------------------------------------
# Test runner
# -------------------------------------------------------------------------


def run_adjustment_test(
    tier: RestTier,
    solver_step_fn: Callable[[HydrostaticState, float], HydrostaticState],
    base: bs.DryBaseState,
    bubble_params: WarmBubbleParams | None = None,
    nx: int = 64,
    ny: int = 64,
    Lx: float = 100_000.0,    # 100 km — small domain for adjustment
    Ly: float = 100_000.0,    # only the inner-domain response matters
    dt: float = 1.0,
    thresholds: AdjustmentPassThresholds | None = None,
) -> AdjustmentTestResult:
    """
    Run the hydrostatic adjustment test at the given tier.
    """
    if bubble_params is None:
        bubble_params = WarmBubbleParams()
    if thresholds is None:
        thresholds = AdjustmentPassThresholds()

    state = initialize_adjustment_state(base, nx, ny, Lx, Ly, bubble_params)
    dx = Lx / nx
    dy = Ly / ny
    dz = base.z[1] - base.z[0]

    for step in range(tier.n_steps):
        state = solver_step_fn(state, dt)

    return evaluate_adjustment(
        state, base, bubble_params,
        dx, dy, dz, Lx, Ly,
        thresholds, tier, tier.n_steps, dt,
    )


# -------------------------------------------------------------------------
# Self-test entry point
# -------------------------------------------------------------------------


if __name__ == "__main__":
    print("=" * 64)
    print("HARNESS SELF-TEST: dry hydrostatic adjustment")
    print("=" * 64)

    # Build a base state (small Lz domain matching adjustment scale)
    nz = 32
    Lz = 10_000.0  # 10 km vertical, deep enough for adjustment to play out
    z = np.linspace(0, Lz, nz)
    base = bs.constant_N_dry_base_state(
        z, N=0.01, theta_surface=300.0,
        staggering=bs.GridStaggering.UNSTAGGERED_PLACEHOLDER,
    )

    # Default Bryan & Fritsch bubble
    bubble = WarmBubbleParams()
    Lx = Ly = 100_000.0
    nx = ny = 64

    # Compute the analytical expected deficit
    expected = expected_pressure_deficit_at_center(base, bubble)
    print(f"\nBubble: B&F default ({bubble.theta_max_K} K, {bubble.xr_m/1000:.0f} km radius, "
          f"{bubble.zc_m/1000:.0f} km altitude)")
    print(f"Expected pressure deficit at center column: {expected:.2f} Pa")
    print(f"  (V7 actual across 15,000 frames: ~0 to 1.9 Pa)")
    print(f"  (V8 must achieve this magnitude within factor of 2 to pass)")

    # Self-test path A: perfect adjustment solver — should pass
    print("\n[A] Perfect adjustment solver (returns analytical equilibrium):")
    print("    expecting PASS — harness recognizes the right answer")

    # Adapter that injects bubble_params into perfect_adjustment_solver
    def perfect_step(state, dt):
        return perfect_adjustment_solver(
            state, dt, bubble_params=bubble, Lx=Lx, Ly=Ly,
        )

    result_a = run_adjustment_test(
        RestTier.SHORT, perfect_step, base, bubble, nx=nx, ny=ny, Lx=Lx, Ly=Ly,
    )
    tag_a = "✓" if result_a.passed else "✗"
    print(f"    measured deficit:  {result_a.measured_deficit_Pa:.2f} Pa")
    print(f"    deficit ratio:     {result_a.deficit_ratio:.3f}  "
          f"(target: 0.5 - 2.0)")
    print(f"    horizontal offset: {result_a.horizontal_offset_radii:.2f} radii")
    print(f"    vertical offset:   {result_a.vertical_offset_radii:.2f} radii")
    print(f"    pass: {tag_a}")
    if result_a.failure_reasons:
        for r in result_a.failure_reasons:
            print(f"      - {r}")

    # Self-test path B: V7-style no-deficit solver — should fail with the
    # V7 failure mode message
    print("\n[B] V7-style no-deficit solver (φ stays zero):")
    print("    expecting FAIL — harness detects the V7 failure mode")

    result_b = run_adjustment_test(
        RestTier.SHORT, v7_style_no_deficit_solver, base, bubble,
        nx=nx, ny=ny, Lx=Lx, Ly=Ly,
    )
    tag_b = "✓" if result_b.passed else "✗"
    print(f"    measured deficit:  {result_b.measured_deficit_Pa:.2f} Pa")
    print(f"    pass: {tag_b}")
    if result_b.failure_reasons:
        for r in result_b.failure_reasons:
            print(f"      - {r}")

    print("\n" + "=" * 64)
    if result_a.passed and not result_b.passed:
        print("HARNESS SELF-TEST: PASSED")
        print("  - Perfect-solver path produces correct deficit")
        print("  - V7-style solver correctly diagnosed as failed adjustment")
        print("  - Ready to receive V8's actual solver step function")
    else:
        print("HARNESS SELF-TEST: FAILED — investigate before plugging in V8")
        print(f"  Path A passed: {result_a.passed} (should be True)")
        print(f"  Path B passed: {result_b.passed} (should be False)")
    print("=" * 64)
