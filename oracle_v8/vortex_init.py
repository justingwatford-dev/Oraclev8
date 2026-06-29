"""
Oracle V8 — TC Vortex Initialization
======================================
Holland (1980) gradient-wind-balanced vortex with hydrostatic θ′ derivation.

Algorithm (Five's recommended approach)
-----------------------------------------
1.  Holland tangential wind profile on a 1-D radial grid:

        V_t(r) = Vmax · (Rmax/r)^B · exp(1 − (Rmax/r)^B)

    Maximum wind at r = Rmax; decays to zero at r = 0 and r → ∞.

2.  Vertical structure — Gaussian peaked at upper troposphere:

        S(z) = exp(−(z − z_peak)² / (2·σ_z²))

    z_peak = 8 km by default (~350 hPa).  S(0) ≈ 0.16 at the surface.
    Full 3-D wind: V_t(r, z) = V_t(r) · S(z)

3.  Gradient-wind balance → P′(r, z) at each level:

        dP′/dr = ρ̄(z) · (V_t²/r + f · V_t)      (cyclostrophic + Coriolis)
        P′(r, z) = −∫_r^R_env  ρ̄(z) · (V_t²/r′ + f·V_t(r′,z))  dr′

    P′(r) < 0 everywhere (pressure deficit); P′(R_env) = 0 by definition.

4.  Hydrostatic balance → θ′(r, z):

        θ′(r, z) = [θ̄(z) / (ρ̄(z) · g)] · ∂P′(r, z)/∂z

    This is the physically correct derivation: θ′ gets a real job —
    supporting the surface pressure deficit hydrostatically.  Do NOT map
    the pressure deficit to θ′ pointwise via the ideal-gas equation of
    state; that double-counts the thermodynamic coupling (Five, ensemble
    review May 2026).

    Sign check: ∂P′/∂z > 0 for a warm-core vortex (deficit weakens with
    height) → θ′ > 0 (warm core) ✓.

5.  Wind decomposition to Cartesian:

        u(x,y,z) = −V_t(r,z) · sin(φ) + u_env
        v(x,y,z) = +V_t(r,z) · cos(φ) + v_env

    where φ = atan2(y−y_c, x−x_c).  NH cyclonic (counterclockwise).

6.  w = 0 at initialisation (secondary circulation spins up from the
    balanced primary vortex during the integration).

Parameters
----------
Vmax, Rmax : from HURDAT2 Best Track at the chosen t=0.
B          : Holland shape parameter — typically 1.5 for Cat 5.
             This is the SINGLE tunable knob.  Freeze before seeing
             landfall errors (Five, ensemble review May 2026).
u_env, v_env : deep-layer mean (DLM) steering wind (m/s).
               Also passed to SpongeDampingComponent so the sponge
               damps perturbations from the steering flow, not from zero.

References
----------
Holland, G. J. (1980). An analytic model of the wind and pressure
profiles in hurricanes. Mon. Wea. Rev., 108, 1212–1218.
"""

from __future__ import annotations

import numpy as np

# Gravitational acceleration (m/s²)
from oracle_v8.constants import GRAVITY as _G


class HollandVortexInit:
    """
    Initialise a TC vortex using the Holland (1980) wind profile.

    Usage
    -----
    init = HollandVortexInit(Vmax=60.0, Rmax=40_000.0, B=1.5, f=5e-5,
                             u_env=-5.0, v_env=1.0)
    state = init.build_state(nx, ny, nz, Lx, Ly, base)

    The returned State has:
      - u, v  : Holland tangential wind + steering flow
      - w     : zero (secondary circulation spins up during integration)
      - theta_prime : warm-core θ′ from hydrostatic balance
      - projection_potential : zero (first step will project)
    """

    def __init__(
        self,
        Vmax: float,
        Rmax: float,
        B: float = 1.5,
        f: float = 5e-5,
        R_env: float = None,
        u_env: float = 0.0,
        v_env: float = 0.0,
        n_radial: int = 4000,
        z_peak: float = 8_000.0,
        sigma_z: float = 4_000.0,
        H_wind: float = 15_000.0,
        wind_taper: bool = False,
        taper_start_frac: float = 0.5,
    ) -> None:
        """
        Parameters
        ----------
        Vmax : float
            Maximum sustained wind speed (m/s).  From HURDAT2.
        Rmax : float
            Radius of maximum winds (m).  From HURDAT2.
        B : float
            Holland shape parameter (dimensionless).  Default 1.5.
            Tunable but must be frozen before seeing landfall errors.
        f : float
            Coriolis parameter (s⁻¹).  Match storm latitude at t=0.
        R_env : float, optional
            Outer radius where P′ vanishes (m).  Default 10×Rmax.
        u_env : float
            Zonal steering wind (m/s).  Added to u; also pass to
            SpongeDampingComponent(u_env=u_env).
        v_env : float
            Meridional steering wind (m/s).
        n_radial : int
            Resolution of 1-D radial integration grid.  Default 4000.
        z_peak : float
            Altitude of maximum warm-core buoyancy (m) for θ′.
            Default 8 000 m (~350 hPa, upper troposphere).
        sigma_z : float
            Gaussian half-width of the warm-core θ′ structure (m).
            Default 4 000 m.
        H_wind : float
            e-folding height (m) for the WIND vertical structure.
            Separate from θ′: the Holland profile represents surface
            winds (S_wind(0) = 1.0), decaying upward.  Default 15 000 m.
        """
        self.Vmax     = float(Vmax)
        self.Rmax     = float(Rmax)
        self.B        = float(B)
        self.f        = float(f)
        self.R_env    = float(R_env) if R_env is not None else 10.0 * Rmax
        self.u_env    = float(u_env)
        self.v_env    = float(v_env)
        self.n_radial = int(n_radial)
        self.z_peak   = float(z_peak)
        self.sigma_z  = float(sigma_z)
        self.H_wind   = float(H_wind)
        # Outer wind taper (V8.5.1): the bare Holland profile decays only as
        # r^(-B/2) (~r^-0.75 for B=1.5) — still 10-27 m/s at 500-2000 km, so it
        # fills the domain and gives an unrealistically large β-drift + taper
        # interaction.  R_env was meant to bound the vortex but only bounds the
        # pressure integral (θ′, which is passive).  With wind_taper=True the
        # tangential wind is smoothly brought to 0 between taper_start_frac·R_env
        # and R_env, so R_env becomes the real vortex-size knob.  Default OFF →
        # bit-identical to the validated runs (Hugo) until explicitly enabled.
        self.wind_taper       = bool(wind_taper)
        self.taper_start_frac = float(taper_start_frac)

        # Pre-build 1-D radial grid (CPU, done once)
        self._r1d = np.linspace(0.0, self.R_env, self.n_radial)

    # ------------------------------------------------------------------
    # 1-D profile builders
    # ------------------------------------------------------------------

    def tangential_wind(self, r: np.ndarray) -> np.ndarray:
        """
        Holland (1980) tangential wind profile.

        V_t(r) = Vmax · (Rmax/r)^B · exp(1 − (Rmax/r)^B)

        Handles r = 0 singularity: V_t(0) = 0.
        Returns array with same shape as r.
        """
        r_safe  = np.where(r > 0.0, r, 1.0)
        x       = (self.Rmax / r_safe) ** self.B
        Vt      = self.Vmax * np.sqrt(x * np.exp(1.0 - x))
        Vt      = np.where(r > 0.0, Vt, 0.0)
        if self.wind_taper:
            # cosine taper: 1 inside r0, smoothly → 0 at R_env, 0 beyond.
            # Radial (axisymmetric) ⇒ multiplies tangential flow → ~no divergence
            # (the pre-balance projection mops up the small residual).
            r0   = self.taper_start_frac * self.R_env
            frac = np.clip((r - r0) / max(self.R_env - r0, 1e-9), 0.0, 1.0)
            Vt   = Vt * 0.5 * (1.0 + np.cos(np.pi * frac))
        return Vt

    def vertical_structure(self, z: np.ndarray) -> np.ndarray:
        """
        Upper-troposphere Gaussian warm-core vertical structure.

        S(z) = exp(−(z − z_peak)² / (2·σ_z²))

        Peaks at z_peak (default 8 km, upper troposphere).
        S(z=0) ≈ 0.16 (small surface contribution) — this is crucial:
        placing the peak at the surface causes ∂P′/∂z to be enormous
        over the first dz, giving unrealistically large θ′ (>40 K).
        With z_peak = 8 km, the surface θ′ is ~5 K, consistent with
        observed TC eye temperature anomalies.
        """
        return np.exp(-0.5 * ((z - self.z_peak) / self.sigma_z) ** 2)

    def wind_vertical_structure(self, z: np.ndarray) -> np.ndarray:
        """
        Vertical structure for the WIND field.

        S_wind(z) = exp(−z / H_wind)

        S_wind(0) = 1.0 — the Holland profile represents surface winds,
        so k=0 must have full Vmax.  Decays upward with e-folding height
        H_wind (default 15 km).  Kept separate from the θ′ warm-core
        structure so the two can be tuned independently.
        """
        return np.exp(-z / self.H_wind)

    # ------------------------------------------------------------------
    # 2.  Pressure deficit from gradient-wind balance
    # ------------------------------------------------------------------

    def _pressure_deficit_1d(self, rho_bar_z: float, S_z: float) -> np.ndarray:
        """
        Compute P′(r) at one vertical level.

        Parameters
        ----------
        rho_bar_z : float
            Base-state density ρ̄ at this level (kg/m³).
        S_z : float
            Vertical structure factor S(z) at this level.

        Returns
        -------
        P_prime : (n_radial,) array, P′ ≤ 0, P′(R_env) = 0.
        """
        r   = self._r1d
        Vt  = self.tangential_wind(r) * S_z         # wind at this level

        # Radial pressure gradient: dP′/dr = ρ̄·(Vt²/r + f·Vt)
        r_safe = np.where(r > 0.0, r, 1.0)
        dPdr   = rho_bar_z * (
            np.where(r > 0.0, Vt**2 / r_safe, 0.0)
            + self.f * Vt
        )

        # Integrate from R_env inward (right-to-left trapezoid rule)
        # P′(R_env) = 0; P′(r) = −∫_r^R_env dP/dr dr
        dr      = r[1] - r[0]                        # uniform spacing
        P_prime = np.zeros(self.n_radial)
        for i in range(self.n_radial - 2, -1, -1):
            P_prime[i] = P_prime[i + 1] - 0.5 * (dPdr[i] + dPdr[i + 1]) * dr

        return P_prime

    # ------------------------------------------------------------------
    # 3.  Full 3-D initialization
    # ------------------------------------------------------------------

    def build_state(
        self,
        nx: int,
        ny: int,
        nz: int,
        Lx: float,
        Ly: float,
        base,
        x_c: float = None,
        y_c: float = None,
    ):
        """
        Build a balanced initial State.

        Parameters
        ----------
        nx, ny, nz : grid dimensions
        Lx, Ly     : domain extent (m)
        base       : base state object with z, rho0, theta0 arrays
        x_c, y_c   : vortex centre (m).  Default: domain centre.

        Returns
        -------
        State  (numpy arrays — integrator wraps them to GPU via wrap_base)
        """
        from oracle_v8.solver.tendency import State
        from oracle_v8.backend import asarray

        # Vortex centre
        x_c = x_c if x_c is not None else Lx / 2.0
        y_c = y_c if y_c is not None else Ly / 2.0

        # Cell-centred coordinate arrays
        dx   = Lx / nx
        dy   = Ly / ny
        x    = (np.arange(nx) + 0.5) * dx         # (nx,)
        y    = (np.arange(ny) + 0.5) * dy         # (ny,)
        X, Y = np.meshgrid(x, y, indexing="ij")   # (nx, ny)

        # Radius and azimuth from vortex centre
        dX   = X - x_c
        dY   = Y - y_c
        R2d  = np.sqrt(dX**2 + dY**2)             # (nx, ny)
        PHI  = np.arctan2(dY, dX)                 # azimuthal angle

        # Base-state arrays (ensure numpy for init math)
        from oracle_v8.backend import to_numpy
        z      = to_numpy(base.z)     # (nz,)
        rho0   = to_numpy(base.rho0)  # (nz,)
        theta0 = to_numpy(base.theta0)  # (nz,)

        # --- Warm-core θ′ from gradient-wind + hydrostatic balance --------
        # Implements steps 3-4 of this module's docstring, replacing the
        # earlier PRESCRIBED Gaussian warm core that did NOT balance the wind
        # (red-team round 3: the Gaussian core was ~40% too wide radially and
        # peaked at 8 km vs the balanced ~1 km, so enabling buoyancy with it
        # would shock the vortex at t=0).
        #
        #   1. P′(r,z): radially integrate the gradient-wind equation
        #        dP′/dr = ρ̄(z)·(V_t²/r + f·V_t)
        #      using the SAME tangential wind the u,v field uses
        #      (Holland × S_wind, including any outer taper).
        #   2. θ′(r,z) = [θ̄(z)/(ρ̄(z)·g)]·∂P′/∂z   (discrete vertical derivative)
        #
        # The radial AND vertical structure of θ′ thus emerge from the wind;
        # nothing is prescribed.  Sign: ∂P′/∂z > 0 (deficit weakens with
        # height) → θ′ > 0 (warm core).
        S_wind = self.wind_vertical_structure(z)          # (nz,) exp from surface

        # P′(r, z_k) at every level on the 1-D radial grid
        Pp = np.zeros((self.n_radial, nz))
        for k in range(nz):
            Pp[:, k] = self._pressure_deficit_1d(float(rho0[k]), float(S_wind[k]))

        # ∂P′/∂z (centered interior, one-sided at surface/lid; handles non-uniform z)
        dPp_dz = np.empty_like(Pp)
        dPp_dz[:, 1:-1] = (Pp[:, 2:] - Pp[:, :-2]) / (z[2:] - z[:-2])[None, :]
        dPp_dz[:, 0]    = (Pp[:, 1]  - Pp[:, 0])  / (z[1]  - z[0])
        dPp_dz[:, -1]   = (Pp[:, -1] - Pp[:, -2]) / (z[-1] - z[-2])

        # Hydrostatic θ′ on the radial grid, then map to the 2-D grid per level
        theta_rad = (theta0[None, :] / (rho0[None, :] * _G)) * dPp_dz   # (n_radial, nz)
        r_flat    = R2d.ravel()
        theta_2d  = np.empty((nx, ny, nz))
        for k in range(nz):
            theta_2d[:, :, k] = np.interp(
                r_flat, self._r1d, theta_rad[:, k]
            ).reshape(nx, ny)

        # Tangential wind → u, v  (each level independently)
        # u = -Vt·sin(φ) + u_env,  v = +Vt·cos(φ) + v_env  (NH cyclonic)
        sin_phi = np.sin(PHI)
        cos_phi = np.cos(PHI)
        u_3d  = np.zeros((nx, ny, nz))
        v_3d  = np.zeros((nx, ny, nz))

        for k in range(nz):
            Vt_k = self.tangential_wind(R2d) * float(S_wind[k])  # surface-peaked
            u_3d[:, :, k] = -Vt_k * sin_phi + self.u_env
            v_3d[:, :, k] = +Vt_k * cos_phi + self.v_env

        # --- Assemble State (move to device) ----------------------------
        return State(
            u                   = asarray(u_3d),
            v                   = asarray(v_3d),
            w                   = asarray(np.zeros((nx, ny, nz + 1))),
            theta_prime         = asarray(theta_2d),
            projection_potential= asarray(np.zeros((nx, ny, nz))),
            t                   = 0.0,
        )

    # ------------------------------------------------------------------
    # Diagnostic helpers
    # ------------------------------------------------------------------

    def surface_pressure_deficit(self, rho_bar_surface: float) -> float:
        """
        Return the central surface pressure deficit P′(r=0) [Pa].
        Negative for a valid cyclonic vortex.
        """
        P1d = self._pressure_deficit_1d(rho_bar_surface, S_z=1.0)
        return float(P1d[0])

    def peak_wind_check(self) -> tuple[float, float]:
        """
        Verify that V_t peaks at Rmax with value Vmax.
        Returns (r_at_peak, V_at_peak).
        """
        Vt  = self.tangential_wind(self._r1d)
        idx = int(np.argmax(Vt))
        return float(self._r1d[idx]), float(Vt[idx])

    def summary(self) -> str:
        """One-line summary for logging."""
        r_peak, v_peak = self.peak_wind_check()
        P_min = self.surface_pressure_deficit(rho_bar_surface=1.15)
        taper = (
            f"wind_taper=ON start={self.taper_start_frac * self.R_env / 1000:.0f} km "
            f"(frac={self.taper_start_frac:.2f})  "
            if self.wind_taper else
            "wind_taper=OFF  "
        )
        return (
            f"HollandVortex: Vmax={self.Vmax:.0f} m/s  "
            f"Rmax={self.Rmax/1000:.0f} km  B={self.B:.2f}  "
            f"f={self.f:.2e}  R_env={self.R_env/1000:.0f} km  "
            f"{taper}"
            f"P_min~{P_min:.0f} Pa  "
            f"z_peak={self.z_peak/1000:.0f} km  sz={self.sigma_z/1000:.0f} km  "
            f"u_env={self.u_env:.1f}  v_env={self.v_env:.1f} m/s"
        )
