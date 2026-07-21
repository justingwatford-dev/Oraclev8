"""
Which stabilizer is load-bearing? Ablate the production SLOW stack on a Cat-4
vortex (barotropic, no cap) and see which removals cause blow-up.

Directly tests the reviewer's "fix advection and delete half these components"
claim: if a component's removal blows the run up, it is load-bearing given the
current centered advection; if not, it is a candidate for removal.
"""
import os
os.environ["ORACLE_GPU"] = "0"

import numpy as np
from oracle_v8.backend import xp
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import RK3Integrator
from oracle_v8.production_config import build_base_state, build_production_config
from dataclasses import replace as dc_replace

NX = NY = 128
NZ = 32
DX = 15_625.0
LX = LY = NX * DX
LZ = 20_000.0
F = 5.7e-5
DT = 30.0
MAXSTEPS = 300
BLOWUP = 250.0

zc, rho0_arr, theta0_arr = build_base_state(NZ, LZ)

class Base:
    z = zc; rho0 = rho0_arr; theta0 = theta0_arr


def vortex():
    init = HollandVortexInit(
        Vmax=64.0, Rmax=75_000.0, B=1.5, f=F,
        R_env=500_000.0, wind_taper=True, taper_start_frac=0.40,
        u_env=0.0, v_env=0.0)
    s = init.build_state(NX, NY, NZ, LX, LY, Base())
    return dc_replace(s, theta_prime=xp.zeros_like(s.theta_prime))


def make_cfg(drop=()):
    cfg = build_production_config(NX, NY, NZ, LX, LY, LZ, F, 0.0, 0.0)
    for slot in drop:
        setattr(cfg, slot, None)
    return cfg


def run(label, drop=()):
    cfg = make_cfg(drop)
    active = [c.name for c in cfg.slow_components()]
    integ = RK3Integrator(config=cfg, base=Base())
    s = vortex()
    last = MAXSTEPS
    for n in range(1, MAXSTEPS + 1):
        s, _ = integ.step(s, dt=DT, step_number=n)
        mu = float(xp.max(xp.abs(s.u)))          # reduce on-device (CuPy-safe)
        if not np.isfinite(mu) or mu > BLOWUP:
            print(f"  {label:32s} BLEW UP at step {n:3d} (max|u|={mu:.0f})   slow={active}")
            return
    mu = float(xp.max(xp.abs(s.u)))
    print(f"  {label:32s} survived {MAXSTEPS}  (max|u|={mu:6.1f})   slow={active}")


if __name__ == "__main__":
    print(f"Grid {NX}x{NY}x{NZ}  dt={DT}s  barotropic, no intensity cap  xp={xp.__name__}\n")
    run("full production stack")
    run("− hyperdiffusion",        drop=("horiz_diffusion",))
    run("− helmholtz divergence",  drop=("divergence_damping",))
    run("− surface drag",          drop=("surface_drag",))
    run("− newtonian cooling",     drop=("newtonian_cooling",))
    run("− hyperdiff − helmholtz", drop=("horiz_diffusion", "divergence_damping"))
