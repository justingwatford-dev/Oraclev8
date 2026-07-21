"""
Verify the WS2002 5th-order upwind advection and run the damper-off ablation
(the reviewer's acceptance test).

(0) Unit checks on _upwind5:
    - smooth field: derivative matches analytic to high order (5th-order accurate)
    - 2Δx sawtooth under pure advection: energy DECREASES (dissipative) for upwind,
      whereas centered forward-Euler GROWS it (the instability).
(1) Ablation: does upwind advection survive with the Helmholtz damper OFF?
    centered+damper(on/off) vs upwind+damper(on/off).
"""
import os
os.environ.setdefault("ORACLE_GPU", "1")

import numpy as np
from oracle_v8.backend import xp, to_numpy
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import RK3Integrator, AdvectionComponent
from oracle_v8.production_config import build_base_state, build_production_config
from dataclasses import replace as dc_replace

NX = NY = 128
NZ = 32
DX = 15_625.0
LX = LY = NX * DX
LZ = 20_000.0
F = 5.7e-5
DT = 30.0
MAXSTEPS = 2000
BLOWUP = 250.0

zc, rho0_arr, theta0_arr = build_base_state(NZ, LZ)

class Base:
    z = zc; rho0 = rho0_arr; theta0 = theta0_arr


# ---------------------------------------------------------------------------
# (0) unit checks
# ---------------------------------------------------------------------------
def unit_checks():
    print("=== (0) _upwind5 unit checks ===")
    # smooth-field accuracy at two resolutions → measure convergence order
    def err(n):
        L = 1.0
        x = (np.arange(n)) * (L / n)
        q = np.sin(2 * np.pi * x / L)[:, None, None] * np.ones((n, 1, 1))
        a = np.ones_like(q)
        dq = to_numpy(AdvectionComponent._upwind5(xp.asarray(q), xp.asarray(a), L / n, 0))
        exact = (2 * np.pi / L) * np.cos(2 * np.pi * x / L)[:, None, None]
        return np.max(np.abs(dq - exact))
    e1, e2 = err(64), err(128)
    print(f"  smooth sin: err(64)={e1:.3e}  err(128)={e2:.3e}  "
          f"order≈{np.log2(e1/e2):.2f} (expect ~5)")

    # 2Δx sawtooth energy under pure advection (forward Euler, uniform a)
    n = 128
    saw = ((-1.0) ** np.arange(n))[:, None, None] * np.ones((n, 1, 1))
    a = 5.0 * np.ones_like(saw)
    dx = LX / n
    dt = 0.1 * dx / 5.0
    def evolve(use_upwind, steps=200):
        q = xp.asarray(saw.copy()); A = xp.asarray(a)
        E0 = float(to_numpy(xp.sum(q ** 2)))
        for _ in range(steps):
            if use_upwind:
                dq = AdvectionComponent._upwind5(q, A, dx, 0)
            else:
                dq = (xp.roll(q, -1, 0) - xp.roll(q, 1, 0)) / (2 * dx)
            q = q - dt * A * dq
        return float(to_numpy(xp.sum(q ** 2))) / E0
    print(f"  2Δx sawtooth energy ratio after 200 steps:")
    print(f"    centered2 : {evolve(False):.3e}   (>1 ⇒ amplifying, unstable)")
    print(f"    upwind5h  : {evolve(True):.3e}   (<1 ⇒ dissipative, stable)")


# ---------------------------------------------------------------------------
# (1) ablation
# ---------------------------------------------------------------------------
def vortex():
    init = HollandVortexInit(
        Vmax=64.0, Rmax=75_000.0, B=1.5, f=F,
        R_env=500_000.0, wind_taper=True, taper_start_frac=0.40,
        u_env=0.0, v_env=0.0)
    s = init.build_state(NX, NY, NZ, LX, LY, Base())
    return dc_replace(s, theta_prime=xp.zeros_like(s.theta_prime))


def run(label, scheme, damper):
    cfg = build_production_config(NX, NY, NZ, LX, LY, LZ, F, 0.0, 0.0)
    cfg.advection = AdvectionComponent(nx=NX, ny=NY, nz=NZ, Lx=LX, Ly=LY, Lz=LZ,
                                       scheme=scheme)
    if not damper:
        cfg.divergence_damping = None
    integ = RK3Integrator(config=cfg, base=Base())
    s = vortex()
    for n in range(1, MAXSTEPS + 1):
        s, _ = integ.step(s, dt=DT, step_number=n)
        mu = float(to_numpy(xp.max(xp.abs(s.u))))
        if not np.isfinite(mu) or mu > BLOWUP:
            print(f"  {label:34s} BLEW UP at step {n:3d} (max|u|={mu:.0f})")
            return
    mu = float(to_numpy(xp.max(xp.abs(s.u))))
    print(f"  {label:34s} survived {MAXSTEPS}  (max|u|={mu:6.1f})")


if __name__ == "__main__":
    print(f"Grid {NX}x{NY}x{NZ}  dt={DT}s  xp={xp.__name__}\n")
    unit_checks()
    print("\n=== (1) damper-off ablation (barotropic Cat-4 vortex) ===")
    run("centered2 + helmholtz ON",  "centered2", True)
    run("centered2 + helmholtz OFF", "centered2", False)
    run("upwind5h  + helmholtz OFF", "upwind5h",  False)
    run("upwind5h  + helmholtz ON",  "upwind5h",  True)
