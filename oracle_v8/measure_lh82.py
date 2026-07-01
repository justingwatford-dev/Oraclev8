"""
LH82 small-perturbation study — Phase 1.

Question: does the diabatically-driven eyewall reach the θ′/θ̄ regime where
LH82's small-perturbation assumption (|θ′| ≪ θ̄) is flagged as questionable,
and how big does the small parameter get as a function of forcing?

Config: the validated stable recipe — LH82 + upwind5h advection + ε=0 (no
Helmholtz damper) + buoyancy + prescribed annular eyewall heating. Sweep Q.
At steady state report max|θ′|, the small parameter max(θ′/θ̄), max|w|, max|u|.

Steady-state estimate: θ′_eq ≈ Q·τ_cool (τ_cool=1800 s), so
  Q = 2.5e-3 → ~4.5 K, 5e-3 → ~9 K, 1e-2 → ~18 K, 1.5e-2 → ~27 K
spanning θ′/θ̄ ≈ 1.5% → 9%.
"""
import os
os.environ.setdefault("ORACLE_GPU", "1")

import numpy as np
from oracle_v8.backend import xp, to_numpy
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import (
    RK3Integrator, BuoyancyComponent, DiabaticHeatingComponent, AdvectionComponent)
from oracle_v8.production_config import (
    build_base_state, build_production_config, build_prebal_config, N_PREBAL)

NX = NY = 128
NZ = 32
DX = 15_625.0
LX = LY = NX * DX
LZ = 20_000.0
F = 5.7e-5
DT = 30.0
NSTEPS = int(os.environ.get("LH82_STEPS", "300"))   # bump on GPU for firmer steady state

zc, rho0_arr, theta0_arr = build_base_state(NZ, LZ)

class Base:
    z = zc; rho0 = rho0_arr; theta0 = theta0_arr

theta0_dev = xp.asarray(theta0_arr)   # θ̄(z) for the small-parameter ratio


def mx(a):
    return float(to_numpy(xp.max(xp.abs(a))))


def run(Q_max):
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
                                       scheme="upwind5h")
    cfg.buoyancy = BuoyancyComponent()
    cfg.diabatic_heating = DiabaticHeatingComponent(
        Q_max=Q_max, r_eyewall=75_000.0, width_r=30_000.0,
        z_peak=5_000.0, width_z=3_000.0,
        nx=NX, ny=NY, nz=NZ, Lx=LX, Ly=LY, Lz=LZ)
    cfg.divergence_damping = None            # ε=0 (upwind supplies the dissipation)
    integ = RK3Integrator(config=cfg, base=Base())

    for n in range(1, NSTEPS + 1):
        state, _ = integ.step(state, dt=DT, step_number=n)
        if mx(state.u) > 400.0 or not np.isfinite(mx(state.u)):
            print(f"  Q={Q_max:.1e}  BLEW UP @step {n}")
            return
    # steady-state small parameter: max over the domain of |θ′| / θ̄(z)
    tp = state.theta_prime
    small_param = float(to_numpy(xp.max(xp.abs(tp) / theta0_dev[None, None, :])))
    print(f"  Q={Q_max:.1e}  max|u|={mx(state.u):6.2f}  max|w|={mx(state.w):6.3f}  "
          f"max|θ′|={mx(tp):6.2f} K  max(θ′/θ̄)={100*small_param:5.2f}%")


if __name__ == "__main__":
    _qenv = os.environ.get("LH82_QLIST")
    Q_LIST = ([float(q) for q in _qenv.split(",")] if _qenv
              else [2.5e-3, 5.0e-3, 1.0e-2, 1.5e-2, 2.0e-2])
    print(f"Grid {NX}x{NY}x{NZ}  dt={DT}s  steps={NSTEPS}  LH82+upwind+ε=0+heating  xp={xp.__name__}")
    print("Phase 1 — eyewall θ′/θ̄ vs heating forcing:")
    for Q in Q_LIST:
        run(Q)
