"""
Standalone test of VariableCoefficientPoissonSolver

Uses the manufactured-solution infrastructure already built in
validation/analytical_solutions.py. The solver is given the analytical
source f and ρ̄(z), and we check:

    1. That the recovered φ matches the manufactured φ at second order
       under grid refinement.
    2. That the compatibility residual at (kx=0, ky=0) is at machine
       precision for a properly-constructed source.
    3. That the discrete operator residual is at machine precision
       (the tridiagonal solve is doing what it claims).

For these tests we use the periodic-density manufactured solution
(ρ̄ = ρ_s [1 + ε cos(2π z/Lz)]) since it's fully periodic and matches
the Neumann-vertical convention without wraparound contamination.

The exponential ρ̄ case is harder to set up cleanly here because the
manufactured φ chosen for that test is z-periodic (cos(2π z/Lz)) which
does NOT satisfy Neumann BCs at the rigid boundaries. That's a separate
analytical-solution exercise: derive a manufactured φ with Neumann z BCs
under exponential ρ̄. We'll add that to analytical_solutions.py in a
follow-up; for now, the periodic-density baseline is the load-bearing
test.
"""

from __future__ import annotations

import sys

import numpy as np

from oracle_v8.solver.poisson import VariableCoefficientPoissonSolver
from oracle_v8.validation import analytical_solutions as ana


def make_grid(nx: int, ny: int, nz: int, Lx: float, Ly: float, Lz: float):
    # Horizontal: periodic, vertex-based (np.linspace endpoint=False is fine
    # for periodic FFT directions).
    x = np.linspace(0, Lx, nx, endpoint=False)
    y = np.linspace(0, Ly, ny, endpoint=False)
    # Vertical: cell-centered Lorenz convention. dz = Lz/nz, cells at
    # z[k] = (k + 0.5)*dz, so cell 0 spans [0, dz] (surface at lower face)
    # and cell nz-1 spans [Lz-dz, Lz] (lid at upper face).
    dz = Lz / nz
    z = (np.arange(nz) + 0.5) * dz
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    return X, Y, Z, z


def test_neumann_compatible_manufactured(verbose: bool = True) -> dict:
    """
    Use the Neumann-compatible manufactured solution. Both boundary
    conditions (φ's Neumann BCs and the solver's Neumann BCs) match,
    so the L2 error should converge cleanly at second order.
    """
    Lx, Ly, Lz = 2_000_000.0, 2_000_000.0, 20_000.0
    epsilon = 0.1

    results = []
    for nx, ny, nz in [(16, 16, 16), (32, 32, 32), (64, 64, 64)]:
        X, Y, Z, z_1d = make_grid(nx, ny, nz, Lx, Ly, Lz)

        phi_true, source_f = ana.neumann_poisson_solution(
            X, Y, Z, Lx, Ly, Lz, epsilon=epsilon,
        )

        rho_bar_full = ana.base_state_density_periodic(
            z_1d, epsilon=epsilon, Lz=Lz,
        )
        rho_bar_half = np.zeros(nz + 1)
        rho_bar_half[1:-1] = 0.5 * (rho_bar_full[:-1] + rho_bar_full[1:])
        # For Neumann BC, ghost ρ̄ at the boundary is the cell-center value
        # (mirror reflection). This is consistent with the Neumann condition
        # on φ — the discrete operator collapses cleanly.
        rho_bar_half[0] = rho_bar_full[0]
        rho_bar_half[-1] = rho_bar_full[-1]

        solver = VariableCoefficientPoissonSolver(nx, ny, nz, Lx, Ly, Lz)
        result = solver.solve(source_f, rho_bar_full, rho_bar_half)

        # Mean-subtract for the gauge
        phi_centered = result.phi - np.mean(result.phi)
        phi_true_centered = phi_true - np.mean(phi_true)

        l2_err = (
            np.sqrt(np.mean((phi_centered - phi_true_centered)**2))
            / np.sqrt(np.mean(phi_true_centered**2))
        )

        results.append({
            "nx": nx, "ny": ny, "nz": nz,
            "l2_error": float(l2_err),
            "compat_residual": result.compatibility_residual,
            "disc_op_residual": result.discrete_operator_residual,
        })

        if verbose:
            print(f"  ({nx:3d},{ny:3d},{nz:3d}): "
                  f"l2_err={l2_err:.3e}  "
                  f"compat={result.compatibility_residual:.3e}  "
                  f"disc_op={result.discrete_operator_residual:.3e}")

    # Compute convergence ratios
    if verbose and len(results) >= 2:
        print(f"  Convergence ratios (target ~4 for second-order):")
        for i in range(1, len(results)):
            ratio = results[i-1]["l2_error"] / results[i]["l2_error"]
            print(f"    {results[i-1]['nx']}->{results[i]['nx']}: "
                  f"{ratio:.2f}x error reduction")

    return {"results": results}


def test_zero_source_zero_phi(verbose: bool = True) -> bool:
    """
    Sanity check: if d=0, then φ should be 0 (modulo gauge).

    Tests that the gauge pinning at (0,0) and the Thomas algorithm
    don't introduce spurious nonzero solutions.
    """
    Lx, Ly, Lz = 2_000_000.0, 2_000_000.0, 20_000.0
    nx, ny, nz = 16, 16, 16
    X, Y, Z, z_1d = make_grid(nx, ny, nz, Lx, Ly, Lz)
    rho_bar_full = ana.base_state_density_periodic(z_1d, epsilon=0.1, Lz=Lz)
    rho_bar_half = np.zeros(nz + 1)
    rho_bar_half[1:-1] = 0.5 * (rho_bar_full[:-1] + rho_bar_full[1:])
    rho_bar_half[0] = rho_bar_full[0]
    rho_bar_half[-1] = rho_bar_full[-1]

    solver = VariableCoefficientPoissonSolver(nx, ny, nz, Lx, Ly, Lz)
    d_zero = np.zeros((nx, ny, nz))
    result = solver.solve(d_zero, rho_bar_full, rho_bar_half)

    # phi should be zero (or pure gauge-constant — but gauge is pinned at 0)
    phi_max = float(np.max(np.abs(result.phi)))
    passed = phi_max < 1e-12

    if verbose:
        print(f"  d=0 case: max|phi| = {phi_max:.3e}, "
              f"compat = {result.compatibility_residual:.3e} "
              f"({'PASS' if passed else 'FAIL'})")
    return passed


def main():
    print("=" * 70)
    print("VARIABLE-COEFFICIENT POISSON SOLVER: STANDALONE TEST")
    print("=" * 70)

    print("\n[A] Zero-source sanity check (d=0 → φ=0)")
    zero_passed = test_zero_source_zero_phi()

    print("\n[B] Neumann-compatible manufactured solution")
    print(f"  Both φ's BCs and the solver's BCs are Neumann; expect")
    print(f"  clean second-order convergence (~4x error reduction per doubling).")
    nb_result = test_neumann_compatible_manufactured()

    print()
    print("=" * 70)

    # Summary judgments
    print("\nSUMMARY:")

    # Compatibility residual
    compat_residuals = [r["compat_residual"] for r in nb_result["results"]]
    max_compat = max(compat_residuals)
    if max_compat < 1e-8:
        print(f"  ✓ compatibility residual at machine precision "
              f"(max {max_compat:.3e}) — solver is well-posed")
    else:
        print(f"  ✗ compatibility residual high ({max_compat:.3e}) — investigate")

    # Discrete operator residual
    disc_residuals = [r["disc_op_residual"] for r in nb_result["results"]]
    max_disc = max(disc_residuals)
    if max_disc < 1e-6:
        print(f"  ✓ discrete operator residual at solver precision "
              f"(max {max_disc:.3e}) — Thomas algorithm working")
    else:
        print(f"  ⚠ discrete operator residual elevated ({max_disc:.3e})")

    # Zero-source case
    if zero_passed:
        print(f"  ✓ d=0 → φ=0 sanity check passed")
    else:
        print(f"  ✗ d=0 → φ=0 sanity check FAILED")

    # Convergence
    if len(nb_result["results"]) >= 2:
        ratios = []
        for i in range(1, len(nb_result["results"])):
            ratio = (
                nb_result["results"][i-1]["l2_error"]
                / nb_result["results"][i]["l2_error"]
            )
            ratios.append(ratio)
        min_ratio = min(ratios)
        # Second-order ideal is 4x; 3.0x is the conventional pass threshold
        if min_ratio > 3.0:
            print(f"  ✓ second-order convergence verified "
                  f"(min ratio {min_ratio:.2f}, target ~4.0)")
            convergence_passed = True
        else:
            print(f"  ⚠ convergence sub-optimal "
                  f"(min ratio {min_ratio:.2f}, target ~4.0)")
            convergence_passed = False
    else:
        convergence_passed = False

    print("=" * 70)

    all_passed = (
        zero_passed
        and max_compat < 1e-8
        and max_disc < 1e-6
        and convergence_passed
    )

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
