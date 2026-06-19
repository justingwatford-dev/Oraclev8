"""
Stratified Poisson Validation Test
====================================

Verifies V8's modified pressure projection in isolation, before any
time integration is involved.

The anelastic system replaces the V7 Poisson equation
        ∇²φ = f
with the variable-coefficient form
        ∇·(ρ₀(z) ∇φ) = f
where φ is the projection potential (NOT meteorological pressure; see
analytical_solutions.py for nomenclature).

This test uses the method of manufactured solutions: we choose a smooth
φ and compute f analytically, then ask the V8 solver to recover φ from
f and the prescribed ρ₀(z). The error between the recovered φ and the
manufactured φ measures the solver's accuracy.

Pass criteria:
  - L2 error < 1e-3 at the (64, 64, 32) reference resolution.
  - Convergence order ≥ 1.6 across the anisotropic grid sequence
    [(32,32,16), (64,64,32), (128,128,64)] (expected: 2.0 for second-order
    finite differences, higher for spectral methods; 0.8× expected used
    as floor).
  - Grids are anisotropic to match V7's nest aspect ratio (Lx:Lz = 100:1)
    so the test exercises the same anisotropic regime V8 will operate in.

Boundary-condition status:
  - x and y: periodic by manufactured construction.
  - z: φ is periodic in z (cos(2π z/Lz)) but the exponential ρ₀(z) is NOT
    physically periodic. This test exercises the inner stencil; once V8's
    rigid-lid (Neumann) vertical boundary scheme is implemented, an
    end-to-end variant must be added that uses the real boundary treatment.
    See also test_periodic_density_poisson (planned), which removes the
    ρ₀ wraparound issue by using a periodic ρ₀ profile.

If V8 fails this test, no other test will be informative until the
projection is fixed.
"""

from __future__ import annotations

import numpy as np

from oracle_v8.validation.test_harness import ValidationTest
from oracle_v8.validation import analytical_solutions as ana


class StratifiedPoissonTest(ValidationTest):

    name = "stratified_poisson"
    expected_convergence_order = 2.0
    pass_threshold_l2 = 1e-3

    # Domain matches V7's ~2000 km × 2000 km × 20 km nest scale,
    # so the test exercises ρ₀(z) variation across the same depth
    # V8 will operate on.
    default_domain = {
        "Lx": 2_000_000.0,   # 2000 km
        "Ly": 2_000_000.0,   # 2000 km
        "Lz": 20_000.0,      # 20 km
    }

    def initial_state(self, grid: dict) -> dict:
        """
        For a Poisson test, "initial state" is the source term f and the
        base-state density ρ₀(z). There is no time evolution.
        """
        X, Y, Z = grid["X"], grid["Y"], grid["Z"]
        Lx = grid["domain"]["Lx"]
        Ly = grid["domain"]["Ly"]
        Lz = grid["domain"]["Lz"]

        phi_true, f = ana.stratified_poisson_solution(X, Y, Z, Lx, Ly, Lz)
        rho0 = ana.base_state_density(Z)

        return {
            "f": f,
            "rho0": rho0,
            "phi": np.zeros_like(f),  # solver will fill this
            "_phi_true": phi_true,    # stashed for comparison; not used by solver
            "_grid": grid,
        }

    def analytical_solution(self, grid: dict, t: float) -> dict:
        """
        Time-independent. t is ignored. Returns the manufactured φ.
        """
        X, Y, Z = grid["X"], grid["Y"], grid["Z"]
        Lx = grid["domain"]["Lx"]
        Ly = grid["domain"]["Ly"]
        Lz = grid["domain"]["Lz"]
        phi_true = ana.manufactured_phi(X, Y, Z, Lx, Ly, Lz)
        return {"phi": phi_true}

    def field_to_compare(self, state: dict) -> np.ndarray:
        return state["phi"]


# -----------------------------------------------------------------------
# Self-test: verify the harness machinery itself is sound, independent
# of any Poisson solver. We do this by feeding the harness a "perfect
# solver" that just returns the analytical solution. If the harness
# reports zero error, the harness arithmetic is correct.
# -----------------------------------------------------------------------


def _perfect_solver_step(state: dict, dt: float) -> dict:
    """A trivial 'solver' that cheats by returning the analytical answer.

    Used to validate the harness, NOT to validate any real physics.
    If the harness reports zero L2 error against this 'solver', we know
    the harness comparison machinery is sound. Then we can trust the
    error measurements when a real solver is plugged in.
    """
    grid = state["_grid"]
    new_state = dict(state)
    new_state["phi"] = state["_phi_true"].copy()
    return new_state


def _identity_solver_step(state: dict, dt: float) -> dict:
    """A 'solver' that does nothing. Used to confirm the harness reports
    a NON-zero error when the answer is wrong (phi stays zero).
    """
    return state


if __name__ == "__main__":
    from pathlib import Path
    test = StratifiedPoissonTest()

    print("=" * 60)
    print("HARNESS SELF-TEST: stratified Poisson")
    print("=" * 60)

    # Test 1: perfect solver should yield zero error.
    result_perfect = test.run(
        solver_step_fn=_perfect_solver_step,
        grid_shape=(32, 32, 16),
        domain=test.default_domain,
        n_steps=1,
        dt=1.0,
    )
    print(f"\n[A] Perfect solver:  L2 = {result_perfect.error_l2:.3e}  "
          f"(expected ~0, must pass)")
    print(f"    passed = {result_perfect.passed}")

    # Test 2: identity solver should yield error ~1.0 (since true answer
    # is order-1 in magnitude and the solver returns zeros).
    result_identity = test.run(
        solver_step_fn=_identity_solver_step,
        grid_shape=(32, 32, 16),
        domain=test.default_domain,
        n_steps=1,
        dt=1.0,
    )
    print(f"\n[B] Identity solver: L2 = {result_identity.error_l2:.3e}  "
          f"(expected ~1.0, must NOT pass)")
    print(f"    passed = {result_identity.passed}")

    # Sanity check: verify the analytical source f is consistent with the
    # manufactured phi at second-order in dx by reconstructing
    #     f = ρ₀ ∇²φ + (dρ₀/dz)(∂φ/∂z)
    # via clean centered finite differences and comparing to the
    # analytical f stored in state["f"]. We test convergence rather than
    # absolute error: if the FD reconstruction converges to the analytical
    # f at the expected 2nd order, the analytical derivation is correct.
    #
    # (We deliberately do NOT compute "∇·(ρ₀ ∇φ)" by taking divergence-of-flux
    #  numerically: that compounds two FD operations and converges much more
    #  slowly, which masks correctness checks at modest resolution.)
    print("\n[C] Manufactured-source consistency check (FD convergence "
          "to analytical f):")
    print(f"    {'shape':<18} {'rel L2 error':<15} {'ratio (target ~4)'}")

    convergence_ratios = []
    prev_err = None
    for check_shape in [(32, 32, 16), (64, 64, 32), (128, 128, 64)]:
        check_grid = test.make_grid(check_shape, test.default_domain)
        check_state = test.initial_state(check_grid)
        phi = check_state["_phi_true"]
        rho0 = check_state["rho0"]
        f_ana = check_state["f"]
        dx_c = check_grid["dx"]
        dy_c = check_grid["dy"]
        dz_c = check_grid["dz"]

        # FD Laplacian, periodic, second-order accurate
        d2x = (np.roll(phi, -1, 0) - 2 * phi + np.roll(phi, 1, 0)) / dx_c**2
        d2y = (np.roll(phi, -1, 1) - 2 * phi + np.roll(phi, 1, 1)) / dy_c**2
        d2z = (np.roll(phi, -1, 2) - 2 * phi + np.roll(phi, 1, 2)) / dz_c**2
        laplacian_fd = d2x + d2y + d2z

        # FD ∂φ/∂z, periodic
        dphi_dz_fd = (np.roll(phi, -1, 2) - np.roll(phi, 1, 2)) / (2 * dz_c)

        # Reconstruct f using clean form: ρ₀∇²φ + (dρ₀/dz) ∂φ/∂z
        # with dρ₀/dz = -ρ₀ / H
        H = 8500.0
        f_fd = rho0 * laplacian_fd + (-rho0 / H) * dphi_dz_fd

        err = (
            np.sqrt(np.mean((f_fd - f_ana) ** 2))
            / np.sqrt(np.mean(f_ana ** 2))
        )
        ratio_str = f"{prev_err / err:.2f}" if prev_err is not None else "—"
        if prev_err is not None:
            convergence_ratios.append(prev_err / err)
        print(f"    {str(check_shape):<18} {err:<15.3e} {ratio_str}")
        prev_err = err

    # Pass criterion: FD reconstruction converges to analytical f at 2nd
    # order (refinement ratios ~4). Tolerance: any ratio > 3.5 is passing.
    convergence_ok = (
        len(convergence_ratios) >= 1
        and all(r > 3.5 for r in convergence_ratios)
    )

    print("\n" + "=" * 60)
    all_ok = (
        result_perfect.passed
        and not result_identity.passed
        and convergence_ok
    )
    if all_ok:
        print("HARNESS SELF-TEST: PASSED")
        print("  - Perfect-solver path: zero error (harness arithmetic OK)")
        print("  - Identity-solver path: large error (harness detects wrong "
              "answers)")
        print("  - Analytical f converges at 2nd order under refinement "
              "(derivation OK)")
    else:
        print("HARNESS SELF-TEST: FAILED — investigate before plugging in V8")
    print("=" * 60)
