"""
Oracle V8 — Production configuration (single source of truth)
=============================================================
The ONE storm-agnostic configuration.  Every storm imports build_production_config()
and the constants below; there is no per-storm copy to drift out of sync (the failure
mode that left Hugo without newtonian_cooling / surface_drag, ignoring REQUIRE_ERA5,
and on the old one-shot tracker).

This module is the literal embodiment of the paper's "storm-agnostic configuration /
no per-storm parameter adjustments" claim: nothing here depends on which storm runs.
Per-storm inputs (init lat/lon/Vmax, ERA5 file) are passed in; the physics is fixed.

Lifted verbatim from the validated run_ivan.py config block (Ivan/Katrina path).
BuoyancyComponent is intentionally OMITTED — runs are barotropic (theta'=0 init);
buoyancy is re-enabled only after 30 h stability is confirmed in a separate study.
"""
from __future__ import annotations

import math

# ---------------------------------------------------------------------------
# Physics constants (storm-agnostic — identical for every run)
# ---------------------------------------------------------------------------
DT          = 30.0          # s  (CFL ~0.115 at Vmax 60; 60 s caused advection instability)
DIAG_EVERY  = 60            # steps (= every 30 min)
N_PREBAL    = 5             # projection-only pre-balance iterations

RMAX_RUN_M       = 75_000.0     # 5×dx — numerically representative (obs Rmax dominated by steering+β)
R_ENV_M          = 500_000.0    # vortex-size knob (with WIND_TAPER)
TAPER_START_M    = 200_000.0
TAPER_START_FRAC = TAPER_START_M / R_ENV_M    # = 0.40 at R_env=500km.  taper-start 200 km is
# set from β-drift physics (gate-beta testbed: in-band 1.5–2.5 m/s, max westward component before
# the core limit), NOT from landfall — see run_translation_test.gate_beta_only and legacy/run_katrina.
WIND_TAPER       = True         # winds → 0 by R_env; sets the real vortex size

NU4         = 3.0e11        # m^4/s — ∇⁴ hyperdiffusion (10× below dx^4/(64 dt) limit)
EPSILON     = 0.5           # Helmholtz divergence-damping fraction per SLOW half-step
TAU_COOL    = 1800.0        # s — Newtonian cooling (30 min; bounds θ' ~5K)
CD          = 1.5e-3        # surface-drag coefficient
H_BL        = 1000.0        # m — boundary-layer depth for drag

VMAX_CAP_MS = 70.0          # m/s ceiling on |V'| (loop-applied; None disables)
TAU_CAP     = 300.0         # s — cap relaxation timescale

REQUIRE_ERA5       = True   # abort rather than silently fall back to constant steering
TIME_VARYING_STEER = True   # relax background toward local ERA5 DLM each diag step
TAU_STEER          = 10800.0  # s — steering relaxation timescale (3 h)

# Vertical grid (storm-agnostic)
LZ = 20_000.0
NZ = 32

# Horizontal grid resolution (fixed; domain SIZE varies per storm geometry)
DX = 15_625.0              # m  (Lx/nx held constant; nx chosen per storm)
STANDARD_NX = (256, 320, 384, 448, 512)   # composite FFT sizes, dx=15.625 km
TAPER_FRAC  = 0.8          # reversed-β taper zone starts at y > 0.8·Ly
DOMAIN_MARGIN_KM = 200.0   # clearance between the storm's outer circulation and the taper
N_STEPS_MARGIN_H = 10.0    # run this many hours past landfall


# ---------------------------------------------------------------------------
# Geometry-derived run parameters (computed, NOT hand-set → storm-agnostic)
# ---------------------------------------------------------------------------

def choose_domain(init_lat: float, threshold_lat: float) -> tuple[int, float]:
    """Smallest standard (nx, Ly) such that the storm's R_env outer circulation at
    threshold latitude stays below the reversed-β taper zone (the contamination that
    invalidated Hugo's old 4000 km run).  Generalizes run_ivan.py's IVAN_DOMAIN logic:

        y_thresh = Ly/2 + (threshold-init)·111 km
        require  y_thresh + R_env + margin < TAPER_FRAC·Ly
        →  Ly > (Δlat_km + R_env_km + margin) / (TAPER_FRAC - 0.5)
    """
    dlat_km = (threshold_lat - init_lat) * 111.0
    ly_needed_m = ((dlat_km + R_ENV_M / 1e3 + DOMAIN_MARGIN_KM)
                   / (TAPER_FRAC - 0.5)) * 1e3
    for nx in STANDARD_NX:
        if nx * DX >= ly_needed_m:
            return nx, nx * DX
    nx = STANDARD_NX[-1]
    return nx, nx * DX


def n_steps_for(landfall_time_h: float) -> int:
    """Steps to run past landfall with the standard margin (same cushion Ivan got)."""
    hours = math.ceil(landfall_time_h + N_STEPS_MARGIN_H)
    return int(round(hours * 3600.0 / DT))


# ---------------------------------------------------------------------------
# Base-state stratification (storm-agnostic; depends only on the vertical grid)
# ---------------------------------------------------------------------------

def build_base_state(nz: int = NZ, Lz: float = LZ):
    """Return (z_centers, rho0, theta0) for the standard anelastic base state."""
    import numpy as np
    from oracle_v8.constants import GRAVITY
    dz = Lz / nz
    z_centers = (np.arange(nz) + 0.5) * dz
    theta0 = 300.0 * np.exp(0.01 ** 2 * z_centers / GRAVITY)
    Pi = np.zeros(nz)
    Pi[0] = 1.0 - (GRAVITY / 1004.5) * z_centers[0] / theta0[0]
    for k in range(nz - 1):
        dl = z_centers[k + 1] - z_centers[k]
        Pi[k + 1] = Pi[k] - (GRAVITY / 1004.5) * (dl / 2.0) * (
            1.0 / theta0[k] + 1.0 / theta0[k + 1])
    p0 = 100_000.0 * Pi ** (1004.5 / 287.04)
    rho0 = p0 / (287.04 * theta0 * Pi)
    return z_centers, rho0, theta0


# ---------------------------------------------------------------------------
# The config factories
# ---------------------------------------------------------------------------

def build_prebal_config(nx, ny, nz, Lx, Ly, Lz):
    """Projection-only config for the pre-balance step (removes init divergence)."""
    from oracle_v8.solver import (
        LH82AnelasticEquationSet, AnelasticProjection, OperatorConfig)
    from oracle_v8.grid.staggering import LorenzStaggering
    return OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        projection=AnelasticProjection(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz),
    )


def build_production_config(nx, ny, nz, Lx, Ly, Lz, f, u_env, v_env):
    """The full production physics — identical for every storm.  Per-storm scalars
    (f, u_env, v_env) and grid dims are passed in; the components and their
    coefficients are fixed.  BuoyancyComponent intentionally omitted (barotropic)."""
    from oracle_v8.solver import (
        LH82AnelasticEquationSet, AdvectionComponent, CoriolisComponent,
        SurfaceDragComponent, HyperDiffusionComponent,
        HelmholtzDivergenceDampingComponent, NewtonianCoolingComponent,
        AnelasticProjection, OperatorConfig)
    from oracle_v8.grid.staggering import LorenzStaggering
    return OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        advection=AdvectionComponent(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz),
        coriolis=CoriolisComponent(f=f, mode="beta_plane", Ly=Ly, ny=ny,
                                   u_env=u_env, v_env=v_env, periodic_taper=True),
        horiz_diffusion=HyperDiffusionComponent(nu4=NU4, Lx=Lx, Ly=Ly, nx=nx, ny=ny),
        divergence_damping=HelmholtzDivergenceDampingComponent(
            epsilon=EPSILON, Lx=Lx, Ly=Ly, nx=nx, ny=ny),
        newtonian_cooling=NewtonianCoolingComponent(tau=TAU_COOL),
        surface_drag=SurfaceDragComponent(Cd=CD, H_bl=H_BL, u_env=u_env, v_env=v_env),
        projection=AnelasticProjection(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz),
    )


if __name__ == "__main__":
    # Show the geometry-derived domain / run length for the three calibration storms.
    print("storm-agnostic domain & run-length rules:")
    for nm, ilat, tlat, lf_h in [("Hugo", 27.2, 32.5, 28.0),
                                 ("Katrina", 24.8, 29.1, 35.17),
                                 ("Ivan", 23.0, 30.0, 42.83)]:
        nx, Ly = choose_domain(ilat, tlat)
        ns = n_steps_for(lf_h)
        print(f"  {nm:<8} init {ilat}°N → thr {tlat}°N : nx={nx} Ly={Ly/1e3:.0f}km, "
              f"N_STEPS={ns} ({ns*DT/3600:.0f}h)")
