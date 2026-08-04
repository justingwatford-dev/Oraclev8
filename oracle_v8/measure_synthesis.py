"""
Synthesis experiment — can upwind advection + a REDUCED-ε Helmholtz damper give
BOTH stability AND a physical eyewall secondary circulation?

Background (all measured this session):
  - damper ε=0.5 is load-bearing for stability (barotropic damper-off blows ~step 314)
  - but ε=0.5 throttles the diabatically-driven updraft ~10× (max|w| 0.18 vs 2.0 m/s)
  - upwind5h advection adds dissipation but alone doesn't replace the damper

Hypothesis: upwind dissipation lets us drop ε low enough to stop throttling the
circulation while still controlling the divergent-mode instability.

Config: balanced vortex + buoyancy + eyewall heating (Q=5e-3), advection=upwind5h,
sweeping ε.  Tracks max|u|, max|w| (the circulation), max|θ′|, and stability.

Bump NSTEPS to ~1000 on the GPU for the definitive stability check (must clear the
~314-step barotropic blow-up point with margin).
"""
import os
os.environ.setdefault("ORACLE_GPU", "1")

import numpy as np
from oracle_v8.backend import xp, to_numpy
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import (
    RK3Integrator, BuoyancyComponent, DiabaticHeatingComponent,
    AdvectionComponent, HelmholtzDivergenceDampingComponent)
from oracle_v8.production_config import (
    build_base_state, build_production_config, build_prebal_config, N_PREBAL)

NX = NY = 128
NZ = 32
DX = 15_625.0
LX = LY = NX * DX
LZ = 20_000.0
F = 5.7e-5
DT = 30.0
NSTEPS = 1000          # bump to ~1000 on GPU
Q_MAX = 5.0e-3

zc, rho0_arr, theta0_arr = build_base_state(NZ, LZ)

class Base:
    z = zc; rho0 = rho0_arr; theta0 = theta0_arr


def mx(a):
    return float(to_numpy(xp.max(xp.abs(a))))


def run(label, scheme, eps):
    init = HollandVortexInit(
        Vmax=64.0, Rmax=75_000.0, B=1.5, f=F,
        R_env=500_000.0, wind_taper=True, taper_start_frac=0.40,
        u_env=0.0, v_env=0.0)
    state = init.build_state(NX, NY, NZ, LX, LY, Base())
    pre = RK3Integrator(config=build_prebal_config(NX, NY, NZ, LX, LY, LZ), base=Base())
    for i in range(N_PREBAL):
        state, _ = pre.step(state, dt=1.0, step_number=i)

    cfg = build_production_config(NX, NY, NZ, LX, LY, LZ, F, 0.0, 0.0)
    cfg.advection = AdvectionComponent(nx=NX, ny=NY, nz=NZ, Lx=LX, Ly=LY, Lz=LZ,
                                       scheme=scheme)
    cfg.buoyancy = BuoyancyComponent()
    cfg.diabatic_heating = DiabaticHeatingComponent(
        Q_max=Q_MAX, r_eyewall=75_000.0, width_r=30_000.0,
        z_peak=5_000.0, width_z=3_000.0,
        nx=NX, ny=NY, nz=NZ, Lx=LX, Ly=LY, Lz=LZ)
    if eps <= 0.0:
        cfg.divergence_damping = None
    else:
        cfg.divergence_damping = HelmholtzDivergenceDampingComponent(
            epsilon=eps, Lx=LX, Ly=LY, nx=NX, ny=NY)
    integ = RK3Integrator(config=cfg, base=Base())

    wmax_peak = 0.0
    for n in range(1, NSTEPS + 1):
        state, _ = integ.step(state, dt=DT, step_number=n)
        mu = mx(state.u)
        wmax_peak = max(wmax_peak, mx(state.w))
        if not np.isfinite(mu) or mu > 400.0:
            print(f"  {label:34s} BLEW UP @step {n:3d}")
            return
    print(f"  {label:34s} survived {NSTEPS}  max|u|={mu:5.1f}  "
          f"max|w|_end={mx(state.w):.3f}  max|w|_peak={wmax_peak:.3f}  "
          f"max|θ′|={mx(state.theta_prime):.2f}")


if __name__ == "__main__":
    print(f"Grid {NX}x{NY}x{NZ}  dt={DT}s  heating Q={Q_MAX} steps={NSTEPS}  xp={xp.__name__}\n")
    print("reference (centered advection):")
    run("centered + ε=0.5 (production)", "centered2", 0.5)
    print("upwind advection, ε sweep:")
    run("upwind + ε=0.5", "upwind5h", 0.5)
    run("upwind + ε=0.1", "upwind5h", 0.1)
    run("upwind + ε=0.0 (off)", "upwind5h", 0.0)
