"""
Oracle V8 Analytical Solutions
===============================

Reference analytical solutions for validation tests. Each function returns
a known answer derived from a primary source. These solutions are the ground
truth against which V8's numerical solver is measured.

Nomenclature note: in the anelastic system, the field obtained by inverting
the variable-coefficient Poisson equation is the *projection potential* φ —
the scalar whose gradient corrects the unprojected velocity to satisfy
∇·(ρ₀**u**) = 0. φ has units of pressure but is not the meteorological
pressure that a barometer measures; that is a separate diagnostic field
in V8. We use "φ" or "projection potential" throughout to avoid the
ambiguity.

CRITICAL: Each solution below must be cross-checked against its cited primary
source before the validation suite is treated as authoritative. The Woodruff
physicist's burden-of-proof challenge requires that the validation does not
rest on LLM-recalled formulations. Citations are provided to make verification
straightforward, not to substitute for it.

Verification status (as of V8.0.0-anelastic project start):
    stratified_poisson_solution (exponential ρ₀, x/y periodic):
        - Math: GREEN. Verified by GPT-5.5 (independent re-derivation),
          Gemini (independent re-derivation), and human author against
          Durran (2010) §7.2 and Roache (1998) §3.12.
        - Boundary conditions: x/y periodic by manufactured construction.
          z is technically NOT periodic in physics: ρ₀(z) = ρ_s exp(-z/H)
          jumps from ρ_s exp(-Lz/H) ≈ 0.1 ρ_s back to ρ_s at the wraparound.
          φ itself is z-periodic (cos(2π z/Lz)), so the harness's np.roll
          z-derivatives don't see a discontinuity in φ — but they do see a
          ρ₀ discontinuity. Use this test for inner-stencil verification only.
          For end-to-end V8 testing once a non-periodic vertical boundary
          scheme is implemented, this test must be replaced or extended.
    periodic_density_poisson_solution (cosine-perturbed ρ₀, fully periodic):
        - Built specifically to remove the ρ₀ wraparound issue above. ρ₀(z)
          is a small periodic perturbation around a constant, so the
          stratification structure of the operator is exercised without
          a wraparound discontinuity. Cleaner mathematical baseline.

Currently implemented:
  - stratified_poisson_solution: exponential ρ₀(z), inner-stencil test only
  - periodic_density_poisson_solution: periodic ρ₀(z), fully-periodic test

Planned:
  - taylor_green_2d
  - hydrostatic_balance_initial_condition
  - hydrostatic_adjustment_warm_bubble
"""

from __future__ import annotations

import numpy as np


# -------------------------------------------------------------------------
# Stratified Poisson — tests the V8 modified pressure projection
# -------------------------------------------------------------------------
#
# Reference: Durran (2010), "Numerical Methods for Fluid Dynamics", §7.2,
# anelastic system. The modified Poisson equation for pressure in an
# anelastic system with base-state density ρ₀(z) is:
#
#     ∇·(ρ₀(z) ∇φ) = f(x, y, z)
#
# where f is the source term (in V8: the divergence of the unprojected
# velocity field, weighted by ρ₀). To validate the solver, we manufacture
# a known φ and compute its f directly, then ask the solver to recover φ
# from f and check the error.
#
# Manufactured solution approach (textbook practice for verification):
# Choose ρ₀(z) = ρ_s exp(-z/H) (atmospheric scale height H ≈ 8500 m)
# Choose φ(x,y,z) = sin(2πx/Lx) sin(2πy/Ly) cos(2πz/Lz)
# Compute f analytically from the LHS.
#
# This is a "method of manufactured solutions" verification — the standard
# technique in CFD code verification. See Roache (1998), "Verification and
# Validation in Computational Science and Engineering", §3.12.
# -------------------------------------------------------------------------


def base_state_density(z: np.ndarray, rho_surface: float = 1.225,
                        scale_height: float = 8500.0) -> np.ndarray:
    """
    Exponential atmospheric base-state density.

    ρ₀(z) = ρ_surface × exp(-z / H)

    Default values:
      ρ_surface = 1.225 kg/m³  (sea level standard atmosphere)
      H         = 8500 m       (atmospheric scale height)

    Returns ρ₀ as an array matching z's shape.
    """
    return rho_surface * np.exp(-z / scale_height)


def manufactured_phi(X: np.ndarray, Y: np.ndarray, Z: np.ndarray,
                      Lx: float, Ly: float, Lz: float) -> np.ndarray:
    """
    Manufactured projection potential φ for the stratified Poisson
    verification.

        φ(x, y, z) = sin(2π x/Lx) · sin(2π y/Ly) · cos(2π z/Lz)

    φ is the scalar whose gradient corrects the unprojected velocity to
    satisfy ∇·(ρ₀**u**) = 0 in the anelastic system. It has units of
    pressure but is not the meteorological pressure field.

    Chosen so that:
      - smooth and infinitely differentiable
      - mean-zero in x, y, and z (compatible with periodic BCs and gauge fix)
      - has nontrivial z-dependence so that the ρ₀(z) coupling is exercised
    """
    return (
        np.sin(2.0 * np.pi * X / Lx)
        * np.sin(2.0 * np.pi * Y / Ly)
        * np.cos(2.0 * np.pi * Z / Lz)
    )


def manufactured_source(X: np.ndarray, Y: np.ndarray, Z: np.ndarray,
                         Lx: float, Ly: float, Lz: float,
                         rho_surface: float = 1.225,
                         scale_height: float = 8500.0) -> np.ndarray:
    """
    Source term f(x,y,z) for the stratified Poisson equation, computed
    analytically from the manufactured φ.

    Equation:        ∇·(ρ₀(z) ∇φ) = f
    Expanded:        ρ₀ ∇²φ + (dρ₀/dz)(∂φ/∂z) = f

    With ρ₀(z) = ρ_s exp(-z/H), we have dρ₀/dz = -ρ₀ / H.

    Expanding ∇²φ for the manufactured φ:
        ∂²φ/∂x² = -(2π/Lx)² φ
        ∂²φ/∂y² = -(2π/Ly)² φ
        ∂²φ/∂z² = -(2π/Lz)² φ

    And ∂φ/∂z = -(2π/Lz) sin(2πx/Lx) sin(2πy/Ly) sin(2πz/Lz)

    So:
        f = ρ₀ [ -(2π/Lx)² - (2π/Ly)² - (2π/Lz)² ] φ
            + (-ρ₀/H) [ -(2π/Lz) sin(2πx/Lx) sin(2πy/Ly) sin(2πz/Lz) ]

    NOTE TO VERIFIER: this derivation should be checked by hand. The signs
    on the second term in particular are easy to flip. Reference: any text
    on variable-coefficient elliptic operators.
    """
    rho0 = base_state_density(Z, rho_surface, scale_height)
    phi = manufactured_phi(X, Y, Z, Lx, Ly, Lz)

    kx = 2.0 * np.pi / Lx
    ky = 2.0 * np.pi / Ly
    kz = 2.0 * np.pi / Lz

    laplacian_phi = -(kx**2 + ky**2 + kz**2) * phi

    # ∂φ/∂z = -kz · sin(kx x) sin(ky y) sin(kz z)
    dphi_dz = (
        -kz
        * np.sin(kx * X)
        * np.sin(ky * Y)
        * np.sin(kz * Z)
    )

    drho0_dz = -rho0 / scale_height

    f = rho0 * laplacian_phi + drho0_dz * dphi_dz
    return f


def stratified_poisson_solution(X: np.ndarray, Y: np.ndarray, Z: np.ndarray,
                                  Lx: float, Ly: float, Lz: float
                                  ) -> tuple[np.ndarray, np.ndarray]:
    """
    Convenience function: returns (phi_analytical, source_f) for use in
    the test driver. The solver is given source_f and ρ₀(Z), and is
    expected to return phi_analytical to within discretization error.

    Returns:
        phi (array): the analytical projection potential
        f (array):   the source term to feed the solver
    """
    phi = manufactured_phi(X, Y, Z, Lx, Ly, Lz)
    f = manufactured_source(X, Y, Z, Lx, Ly, Lz)
    return phi, f


# -------------------------------------------------------------------------
# Periodic-density Poisson — fully periodic in all three directions
# -------------------------------------------------------------------------
#
# Purpose: same operator, ∇·(ρ₀(z) ∇φ) = f, but with a base-state density
# that is itself periodic in z. This gives a fully periodic manufactured
# problem with no wraparound discontinuity in ρ₀, useful as a cleaner
# baseline test of the variable-coefficient elliptic operator.
#
# Density profile: ρ₀(z) = ρ_s [1 + ε cos(2π z / Lz)]  with small ε (default 0.1)
# This preserves the variable-coefficient structure of the operator
# (dρ₀/dz ≠ 0 except at the extrema of the cosine) without introducing a
# discontinuity at the z-boundary.
#
# Manufactured projection potential: same form as the stratified case,
#     φ(x,y,z) = sin(kx x) sin(ky y) cos(kz z)
# with kx=2π/Lx, ky=2π/Ly, kz=2π/Lz.
#
# Source derivation: same expansion,
#     f = ρ₀ ∇²φ + (dρ₀/dz)(∂φ/∂z)
# with
#     dρ₀/dz = -ρ_s ε (2π/Lz) sin(2π z / Lz)
# and the same ∇²φ and ∂φ/∂z as in the exponential case.
# -------------------------------------------------------------------------


def base_state_density_periodic(z: np.ndarray, rho_surface: float = 1.225,
                                 epsilon: float = 0.1,
                                 Lz: float = 20_000.0) -> np.ndarray:
    """
    Periodic base-state density for fully-periodic Poisson verification.

    ρ₀(z) = ρ_surface × [1 + ε cos(2π z / Lz)]

    epsilon controls how strongly the variable-coefficient structure is
    exercised: ε=0 reduces to the constant-density Laplacian; ε=0.1 gives
    a 10% periodic density variation; ε too large risks ρ₀ approaching zero
    and ill-conditioning the operator. Default 0.1 is a healthy mid-range.
    """
    return rho_surface * (1.0 + epsilon * np.cos(2.0 * np.pi * z / Lz))


def manufactured_source_periodic(X: np.ndarray, Y: np.ndarray, Z: np.ndarray,
                                   Lx: float, Ly: float, Lz: float,
                                   rho_surface: float = 1.225,
                                   epsilon: float = 0.1) -> np.ndarray:
    """
    Source term f for the periodic-density variable-coefficient Poisson.

    With ρ₀(z) = ρ_s [1 + ε cos(2π z/Lz)] and
         φ = sin(kx x) sin(ky y) cos(kz z), kx=2π/Lx, ky=2π/Ly, kz=2π/Lz:

        f = ρ₀ ∇²φ + (dρ₀/dz)(∂φ/∂z)

    where:
        ∇²φ = -(kx² + ky² + kz²) φ
        dρ₀/dz = -ρ_s ε (2π/Lz) sin(2π z/Lz) = -ρ_s ε kz sin(kz z)
        ∂φ/∂z = -kz sin(kx x) sin(ky y) sin(kz z)

    Therefore:
        f = -(kx²+ky²+kz²) ρ₀ φ
            + (-ρ_s ε kz sin(kz z)) · (-kz sin(kx x) sin(ky y) sin(kz z))
        f = -(kx²+ky²+kz²) ρ₀ φ
            + ρ_s ε kz² sin(kx x) sin(ky y) sin²(kz z)

    NOTE TO VERIFIER: as with the exponential case, this derivation should
    be hand-checked. The signs in the second term in particular are easy to
    flip. Reference: Durran (2010) §7.2 for the operator form; Roache (1998)
    §3.12 for the manufactured-solutions methodology.
    """
    rho0 = base_state_density_periodic(Z, rho_surface, epsilon, Lz)
    phi = manufactured_phi(X, Y, Z, Lx, Ly, Lz)

    kx = 2.0 * np.pi / Lx
    ky = 2.0 * np.pi / Ly
    kz = 2.0 * np.pi / Lz

    laplacian_phi = -(kx**2 + ky**2 + kz**2) * phi

    # ∂φ/∂z = -kz · sin(kx x) sin(ky y) sin(kz z)
    dphi_dz = (
        -kz
        * np.sin(kx * X)
        * np.sin(ky * Y)
        * np.sin(kz * Z)
    )

    # dρ₀/dz = -ρ_s ε kz sin(kz z)
    drho0_dz = -rho_surface * epsilon * kz * np.sin(kz * Z)

    f = rho0 * laplacian_phi + drho0_dz * dphi_dz
    return f


def periodic_density_poisson_solution(X: np.ndarray, Y: np.ndarray, Z: np.ndarray,
                                        Lx: float, Ly: float, Lz: float,
                                        epsilon: float = 0.1
                                        ) -> tuple[np.ndarray, np.ndarray]:
    """
    Convenience function: returns (phi_analytical, source_f) for the
    periodic-density verification. Fully periodic in x, y, and z — no
    boundary discontinuities. Use this as the cleaner mathematical baseline
    for V8 elliptic operator verification.
    """
    phi = manufactured_phi(X, Y, Z, Lx, Ly, Lz)
    f = manufactured_source_periodic(X, Y, Z, Lx, Ly, Lz, epsilon=epsilon)
    return phi, f


# -------------------------------------------------------------------------
# Neumann-compatible manufactured solution
# -------------------------------------------------------------------------
#
# A manufactured solution that satisfies dφ/dz = 0 at z=0 AND z=Lz, so it
# matches the V8 Poisson solver's rigid-vertical-boundary configuration.
#
# Choice:
#     φ(x, y, z) = sin(2π x/Lx) · sin(2π y/Ly) · cos(π z/Lz)
#
# Boundary check: ∂φ/∂z = -(π/Lz) sin(2πx/Lx) sin(2πy/Ly) sin(π z/Lz)
#     At z=0:  sin(0) = 0 ✓
#     At z=Lz: sin(π) = 0 ✓
#
# Density profile: cosine-perturbed periodic ρ̄ (same as before, since the
# operator structure is the same; only φ's z-dependence has changed).
#
# Source derivation: same expansion,
#     f = ρ̄ ∇²φ + (dρ̄/dz)(∂φ/∂z)
# with
#     ∇²φ = -[(2π/Lx)² + (2π/Ly)² + (π/Lz)²] φ
#     ∂φ/∂z = -(π/Lz) sin(kx x) sin(ky y) sin(π z/Lz)
#     dρ̄/dz = -ρ_s ε (2π/Lz) sin(2π z/Lz)
# -------------------------------------------------------------------------


def manufactured_phi_neumann(X: np.ndarray, Y: np.ndarray, Z: np.ndarray,
                              Lx: float, Ly: float, Lz: float) -> np.ndarray:
    """
    Manufactured projection potential with Neumann BCs at z=0 and z=Lz.

        φ(x, y, z) = sin(2π x/Lx) · sin(2π y/Ly) · cos(π z/Lz)

    Note kz = π/Lz (one half-wavelength across the domain), not 2π/Lz.
    This is the simplest non-trivial choice satisfying Neumann at both
    z-boundaries.
    """
    return (
        np.sin(2.0 * np.pi * X / Lx)
        * np.sin(2.0 * np.pi * Y / Ly)
        * np.cos(np.pi * Z / Lz)
    )


def manufactured_source_neumann(X: np.ndarray, Y: np.ndarray, Z: np.ndarray,
                                  Lx: float, Ly: float, Lz: float,
                                  rho_surface: float = 1.225,
                                  epsilon: float = 0.1) -> np.ndarray:
    """
    Source term f for the Neumann-compatible manufactured solution
    paired with the periodic-density base state ρ̄(z) = ρ_s [1 + ε cos(2π z/Lz)].

    Derivation (verified manually before becoming load-bearing):
        ∇²φ = -[(2π/Lx)² + (2π/Ly)² + (π/Lz)²] φ
        ∂φ/∂z = -(π/Lz) sin(kx x) sin(ky y) sin(π z/Lz)
        dρ̄/dz = -ρ_s ε (2π/Lz) sin(2π z/Lz)
        f = ρ̄ · ∇²φ + (dρ̄/dz)(∂φ/∂z)
    """
    rho0 = base_state_density_periodic(Z, rho_surface, epsilon, Lz)
    phi = manufactured_phi_neumann(X, Y, Z, Lx, Ly, Lz)

    kx = 2.0 * np.pi / Lx
    ky = 2.0 * np.pi / Ly
    kz = np.pi / Lz

    laplacian_phi = -(kx**2 + ky**2 + kz**2) * phi

    dphi_dz = (
        -kz
        * np.sin(kx * X)
        * np.sin(ky * Y)
        * np.sin(kz * Z)
    )

    drho0_dz = -rho_surface * epsilon * (2.0 * np.pi / Lz) * np.sin(
        2.0 * np.pi * Z / Lz
    )

    f = rho0 * laplacian_phi + drho0_dz * dphi_dz
    return f


def neumann_poisson_solution(X: np.ndarray, Y: np.ndarray, Z: np.ndarray,
                              Lx: float, Ly: float, Lz: float,
                              epsilon: float = 0.1
                              ) -> tuple[np.ndarray, np.ndarray]:
    """
    Convenience function: returns (phi_analytical, source_f) for the
    Neumann-compatible verification. Use this for measuring the V8
    Poisson solver's convergence under refinement, since both BCs match.
    """
    phi = manufactured_phi_neumann(X, Y, Z, Lx, Ly, Lz)
    f = manufactured_source_neumann(X, Y, Z, Lx, Ly, Lz, epsilon=epsilon)
    return phi, f
