"""
Oracle V8 Base States
======================

Thermodynamically closed dry-atmosphere base states for validation tests.

A "base state" is a horizontally homogeneous, time-independent reference
atmosphere (θ₀(z), Π₀(z), p₀(z), T₀(z), ρ₀(z), u₀=v₀=w₀=0) that satisfies
hydrostatic balance exactly within the discrete operators V8 uses. It is
the rest-state reference about which V8 evolves perturbations.

KEY DESIGN PRINCIPLE: do not derive base-state thermodynamic variables by
analytical integration and then sample them onto V8's grid. Construct them
by integrating V8's own discrete vertical operator. This guarantees that
the base state is in V8's notion of discrete hydrostatic balance to machine
precision, regardless of what integration scheme V8 uses.

Until V8's vertical integration operator exists, this module uses a
trapezoidal-rule placeholder. ALL CALLERS MUST CHECK the
`integration_scheme` attribute of the returned base state and refuse to
treat the test as load-bearing if it reads "trapezoidal_placeholder".

Currently implemented:
  - constant_N_dry_base_state: stable dry atmosphere with constant
    Brunt-Väisälä frequency N. Standard textbook setup for anelastic
    hydrostatic-adjustment tests (Ogura & Phillips 1962; Durran 1989).

Planned:
  - isothermal_dry_base_state (sanity-check fallback)
  - us_standard_atmosphere_dry (realistic profile, piecewise)
  - moist_constant_N (when moist tests are added)

Verification:
  - The math for the constant-N profile must be cross-checked against
    Durran (2010), "Numerical Methods for Fluid Dynamics", §2.2 and §7.
    See also Ogura & Phillips (1962), "Scale Analysis of Deep and Shallow
    Convection in the Atmosphere", J. Atmos. Sci. 19, 173-179.
  - Independent re-derivation by GPT-5.5 and Gemini, plus human author
    cross-check, is the verification standard before a base state becomes
    load-bearing for V8 testing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np


# Standard dry-air thermodynamic constants. Sources:
#   - Cp, Rd: AMS Glossary of Meteorology, "specific heat" and "gas constant"
#   - g, p_ref: WMO standard atmosphere
# All of these now come from V8's canonical constants module (single source
# of truth); aliased to the local names this module already uses.
from oracle_v8.constants import (
    GRAVITY,
    C_P as DRY_AIR_CP,
    R_D as DRY_AIR_RD,
    P_REF as P_REFERENCE,
)


class GridStaggering(Enum):
    """
    Vertical staggering options for the base-state grid.

    LORENZ:   θ on full (mass) levels k=0..nz-1; w on half levels k+1/2.
              θ and ρ co-located on full levels.
    CHARNEY_PHILLIPS: θ on half levels; w on half levels co-located with θ;
                       ρ on full levels.
    UNSTAGGERED_PLACEHOLDER: all variables at the same nz cell centers.
                              Use only until V8 commits to a real staggering.

    V8's choice will determine which discrete hydrostatic relation the base
    state must satisfy. Until V8 commits, this module supports
    UNSTAGGERED_PLACEHOLDER and constructs the base state on a single grid.
    Real staggering support is a forward-looking dependency.
    """
    LORENZ = "lorenz"
    CHARNEY_PHILLIPS = "charney_phillips"
    UNSTAGGERED_PLACEHOLDER = "unstaggered_placeholder"


@dataclass
class DryBaseState:
    """
    Container for a thermodynamically closed dry-atmosphere base state.

    All fields are 1D arrays in z (the base state is horizontally
    homogeneous). Callers tile to 3D as needed for tests.

    Attributes:
        z:       vertical coordinate, shape (nz,) or (nz+1,) depending on
                 staggering; meters above surface.
        theta0:  potential temperature, K
        Pi0:     Exner function, dimensionless: Π = (p/p_ref)^(Rd/cp)
        p0:      pressure, Pa
        T0:      temperature, K
        rho0:    density, kg/m³
        N:       Brunt-Väisälä frequency, 1/s (scalar; constant for
                 constant-N profile)
        staggering: GridStaggering enum
        integration_scheme: string identifier; "trapezoidal_placeholder"
                            until V8's discrete operator is available
        params: dictionary of profile parameters (theta_surface, p_surface,
                etc.) for documentation
    """
    z: np.ndarray
    theta0: np.ndarray
    Pi0: np.ndarray
    p0: np.ndarray
    T0: np.ndarray
    rho0: np.ndarray
    N: float
    staggering: GridStaggering
    integration_scheme: str
    params: dict = field(default_factory=dict)

    def is_load_bearing_ready(self) -> bool:
        """
        Returns True only when this base state can be trusted as a V8
        validation reference. Currently always False because we use a
        trapezoidal placeholder; will return True once V8's discrete
        vertical integration operator is plugged in.
        """
        return self.integration_scheme != "trapezoidal_placeholder"


# -------------------------------------------------------------------------
# Constant-N dry base state
# -------------------------------------------------------------------------
#
# Derivation (dry ideal gas, constant Brunt-Väisälä frequency N):
#
# 1. Brunt-Väisälä frequency for dry air:
#       N² = (g / θ₀) (dθ₀/dz)
#    For constant N, this integrates to:
#       θ₀(z) = θ_surface · exp(N² z / g)                              (1)
#
# 2. Hydrostatic balance in Exner-form for dry air:
#       dΠ/dz = -g / (cp · θ₀(z))                                      (2)
#    Integrate (2) from surface (Π = 1, since p = p_ref at z=0 by choice
#    of p_surface = p_ref):
#       Π(z) = 1 - (g/cp) ∫₀ᶻ dz' / θ₀(z')
#
#    With θ₀(z') = θ_s exp(N² z'/g), the integral has a closed analytical
#    form:
#       ∫₀ᶻ exp(-N² z'/g) / θ_s dz' = (g / (θ_s N²)) [1 - exp(-N² z / g)]
#    so:
#       Π(z) = 1 - (g²/(cp θ_s N²)) [1 - exp(-N² z / g)]               (3)
#
#    HOWEVER: per Five's Point 1 and our extension of it, we do NOT use
#    (3) directly. We use it only to verify the discrete answer is close.
#    The actual base state is built by discrete integration of (2) using
#    the placeholder trapezoidal rule, which will be replaced by V8's
#    operator when available. This is the correct posture: the test base
#    state must be in V8's notion of hydrostatic balance, not in the
#    continuous notion.
#
# 3. Diagnose remaining variables from Π and θ₀:
#       p₀(z) = p_ref · Π(z)^(cp/Rd)                                   (4)
#       T₀(z) = θ₀(z) · Π(z)                                            (5)
#       ρ₀(z) = p₀(z) / (Rd · T₀(z))                                   (6)
#
# Reference: Durran (2010), "Numerical Methods for Fluid Dynamics",
# §2.2 (potential temperature and Exner function) and §7.2 (anelastic
# system). Ogura & Phillips (1962) for the original anelastic derivation.
# -------------------------------------------------------------------------


def _trapezoidal_integrate_dPi_dz(
    z: np.ndarray, theta0: np.ndarray
) -> np.ndarray:
    """
    Placeholder discrete integration of dΠ/dz = -g/(cp θ₀) from surface.

    Trapezoidal rule:
        Π[k+1] = Π[k] - (g / cp) · (dz/2) · (1/θ₀[k] + 1/θ₀[k+1])

    THIS IS A PLACEHOLDER. When V8's discrete vertical integration
    operator becomes available, replace this function with one that
    invokes V8's operator directly. Failure to do this swap before
    treating the test as load-bearing will produce O(dz²) spurious motion
    when V8's integration scheme differs from trapezoidal.
    """
    nz = len(z)
    Pi = np.zeros(nz)
    Pi[0] = 1.0  # Π at surface (corresponds to p = p_ref)
    for k in range(nz - 1):
        dz_local = z[k + 1] - z[k]
        Pi[k + 1] = Pi[k] - (GRAVITY / DRY_AIR_CP) * (dz_local / 2.0) * (
            1.0 / theta0[k] + 1.0 / theta0[k + 1]
        )
    return Pi


def constant_N_dry_base_state(
    z: np.ndarray,
    N: float = 0.01,
    theta_surface: float = 300.0,
    staggering: GridStaggering = GridStaggering.UNSTAGGERED_PLACEHOLDER,
) -> DryBaseState:
    """
    Construct a dry, hydrostatically balanced, constant-N base state on
    a 1D vertical grid z.

    Parameters
    ----------
    z : np.ndarray
        Vertical grid points, shape (nz,), meters above surface,
        monotonically increasing, z[0] should be 0 (surface).
    N : float
        Brunt-Väisälä frequency in 1/s. Default 0.01 corresponds to a
        moderately stable troposphere; typical values are 0.005 (weakly
        stable, near tropics) to 0.02 (strongly stable stratosphere).
    theta_surface : float
        Surface potential temperature, K. Default 300 (warm tropical
        ocean surface).
    staggering : GridStaggering
        Staggering choice. UNSTAGGERED_PLACEHOLDER is the only option
        currently fully supported; LORENZ and CHARNEY_PHILLIPS are stubs
        for future V8 integration.

    Returns
    -------
    DryBaseState with all fields populated. The returned object's
    `is_load_bearing_ready()` method returns False until V8's discrete
    vertical integration operator replaces the trapezoidal placeholder.
    """
    if staggering != GridStaggering.UNSTAGGERED_PLACEHOLDER:
        raise NotImplementedError(
            f"Staggering {staggering.value} not yet supported. Only "
            f"UNSTAGGERED_PLACEHOLDER is implemented until V8 commits to "
            f"a vertical staggering."
        )

    if z[0] != 0.0:
        raise ValueError(
            f"Base state requires z[0] = 0 (surface); got z[0] = {z[0]}"
        )

    # Step 1: θ₀(z) from constant-N relation
    theta0 = theta_surface * np.exp(N**2 * z / GRAVITY)

    # Step 2: Π(z) by discrete integration (placeholder!)
    Pi0 = _trapezoidal_integrate_dPi_dz(z, theta0)

    # Sanity check: Π must remain positive over the domain. If it crosses
    # zero, the chosen N or domain depth are nonphysical.
    if np.any(Pi0 <= 0):
        raise ValueError(
            f"Π became non-positive somewhere in the domain. Check that "
            f"N={N}, theta_surface={theta_surface}, and Lz={z[-1]} are "
            f"consistent with a physical atmosphere. Min Π = {np.min(Pi0)}"
        )

    # Steps 3-5: diagnose remaining variables
    p0 = P_REFERENCE * Pi0 ** (DRY_AIR_CP / DRY_AIR_RD)
    T0 = theta0 * Pi0
    rho0 = p0 / (DRY_AIR_RD * T0)

    return DryBaseState(
        z=z,
        theta0=theta0,
        Pi0=Pi0,
        p0=p0,
        T0=T0,
        rho0=rho0,
        N=N,
        staggering=staggering,
        integration_scheme="trapezoidal_placeholder",
        params={
            "theta_surface": theta_surface,
            "p_surface": P_REFERENCE,
            "N": N,
            "domain_top_z": float(z[-1]),
            "n_levels": len(z),
        },
    )


def constant_N_continuous_Pi(z: np.ndarray, N: float, theta_surface: float) -> np.ndarray:
    """
    Continuous analytical Π(z) for the constant-N profile.

        Π(z) = 1 - (g²/(cp θ_s N²)) [1 - exp(-N² z / g)]

    Used for measuring the discretization error of the placeholder
    integration. NOT used to build the base state directly.
    """
    return 1.0 - (
        GRAVITY**2 / (DRY_AIR_CP * theta_surface * N**2)
    ) * (1.0 - np.exp(-N**2 * z / GRAVITY))


def discrete_vs_continuous_Pi_error(
    z: np.ndarray, N: float, theta_surface: float
) -> dict:
    """
    Diagnostic: compare the placeholder discrete Π to the continuous
    analytical Π. Returns max and L2 errors.

    This quantifies how far our trapezoidal placeholder is from the
    continuous solution at a given grid resolution. When V8's discrete
    operator replaces the placeholder, this same function can compare
    V8's discrete answer to the continuous analytical answer to measure
    V8's vertical integration accuracy.
    """
    theta0 = theta_surface * np.exp(N**2 * z / GRAVITY)
    Pi_discrete = _trapezoidal_integrate_dPi_dz(z, theta0)
    Pi_continuous = constant_N_continuous_Pi(z, N, theta_surface)
    diff = Pi_discrete - Pi_continuous
    return {
        "max_abs_error": float(np.max(np.abs(diff))),
        "rms_error": float(np.sqrt(np.mean(diff**2))),
        "max_rel_error": float(np.max(np.abs(diff / Pi_continuous))),
        "Pi_discrete": Pi_discrete,
        "Pi_continuous": Pi_continuous,
    }
