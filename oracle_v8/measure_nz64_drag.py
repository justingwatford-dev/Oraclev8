"""
NZ=64 drag-artifact confirmation — Run 3, bird 2 (AM-budget study).

Run 2's lesson: the doubled column drag at nz=64 (S = 0.75 → 1.594; analytic +
component-verified) is invisible to vortex intensity in a DECAYING quiescent
vortex — surface drag couples to intensity only through a driven secondary
circulation.  The regime where the spin-down was OBSERVED is the heated,
buoyancy-on LH82 Phase-3B configuration (Q=1e-2, upwind5h, ε=0, dt=15), where
the heating-driven boundary-layer inflow is exactly the pathway a doubled BL
sink would strangle.  So the kill shot runs THERE.

2×2 design (all Q=1e-2, ε=0, dt=15, same physical time as Phase 3B):

    nz=32  historical    ← Phase-3B reference row (recorded max|u| = 84.3)
    nz=32  NORMALIZED    ← the 0.75→1.00 column-drag effect, in-regime
    nz=64  historical    ← reproduce the spin-down (recorded max|u| = 48.1)
    nz=64  NORMALIZED    ← THE KILL SHOT: drag fix at dz/2

Registered predictions: ENVELOPE_INTENSIFICATION.md "Run 3" (frozen before
any GPU run).  Metric = max|u| (the Phase-3B instrument, for direct
comparability; Run-1/2 instrument lesson respected — plus max|w| and θ′/θ̄).

Usage:  $env:LH82_STEPS = 1000 ; python -m oracle_v8.measure_nz64_drag
"""
import os
os.environ.setdefault("ORACLE_GPU", "1")

import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import time

import numpy as np
from oracle_v8.backend import xp, to_numpy
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import (
    RK3Integrator, BuoyancyComponent, DiabaticHeatingComponent,
    AdvectionComponent, SurfaceDragComponent)
from oracle_v8.production_config import (
    build_production_config, build_prebal_config, N_PREBAL, CD, H_BL)
from oracle_v8.measure_lh82_phase3 import Grid, LX, LY, LZ, F, mx


def run(g, Q_max, dt, nsteps, drag_normalized):
    init = HollandVortexInit(
        Vmax=64.0, Rmax=75_000.0, B=1.5, f=F,
        R_env=500_000.0, wind_taper=True, taper_start_frac=0.40,
        u_env=0.0, v_env=0.0)
    state = init.build_state(g.nx, g.ny, g.nz, LX, LY, g.base)
    pre = RK3Integrator(config=build_prebal_config(g.nx, g.ny, g.nz, LX, LY, LZ),
                        base=g.base)
    for i in range(N_PREBAL):
        state, _ = pre.step(state, dt=1.0, step_number=i)

    cfg = build_production_config(g.nx, g.ny, g.nz, LX, LY, LZ, F, 0.0, 0.0)
    cfg.advection = AdvectionComponent(nx=g.nx, ny=g.ny, nz=g.nz,
                                       Lx=LX, Ly=LY, Lz=LZ, scheme="upwind5h")
    cfg.buoyancy = BuoyancyComponent()
    cfg.diabatic_heating = DiabaticHeatingComponent(
        Q_max=Q_max, r_eyewall=75_000.0, width_r=30_000.0,
        z_peak=5_000.0, width_z=3_000.0,
        nx=g.nx, ny=g.ny, nz=g.nz, Lx=LX, Ly=LY, Lz=LZ)
    cfg.divergence_damping = None                     # ε=0, the Phase-3B recipe
    cfg.surface_drag = SurfaceDragComponent(
        Cd=CD, H_bl=H_BL, u_env=0.0, v_env=0.0,
        column_normalized=drag_normalized)
    integ = RK3Integrator(config=cfg, base=g.base)

    heartbeat = max(nsteps // 4, 1)
    for n in range(1, nsteps + 1):
        state, _ = integ.step(state, dt=dt, step_number=n)
        if mx(state.u) > 400.0 or not np.isfinite(mx(state.u)):
            return dict(blew=n, tphys=n * dt / 3600.0)
        if n % heartbeat == 0 and n < nsteps:
            small_now = float(to_numpy(xp.max(xp.abs(state.theta_prime)
                                              / g.theta0_d)))
            print(f"      … step {n}/{nsteps}  max|u|={mx(state.u):6.1f}  "
                  f"max(θ′/θ̄)={100*small_now:5.2f}%", flush=True)
    small, flux_ratio = g.neglected_diag(state)
    return dict(blew=None, max_u=mx(state.u), max_w=mx(state.w),
                small=small, flux_ratio=flux_ratio)


if __name__ == "__main__":
    N30 = int(os.environ.get("LH82_STEPS", "1000"))
    N15 = 2 * N30
    t0 = time.time()
    print(f"NZ=64 drag confirmation  xp={xp.__name__}  N15={N15}  "
          f"(heated Phase-3B regime, Q=1e-2, ε=0, dt=15)")
    print("  refs (Phase 3B, historical drag): nz=32 max|u|=84.3 ; nz=64 48.1")
    print("  analytic column-drag factor: 0.75 (nz=32) vs 1.594 (nz=64)")

    rows = []
    for lbl, dims, norm in (
            ("nz=32  historical  (ref 84.3)", (128, 128, 32), False),
            ("nz=32  NORMALIZED",             (128, 128, 32), True),
            ("nz=64  historical  (ref 48.1)", (128, 128, 64), False),
            ("nz=64  NORMALIZED  (kill shot)", (128, 128, 64), True)):
        print(f"\n  {lbl}:")
        r = run(Grid(*dims), 1.0e-2, 15.0, N15, norm)
        rows.append((lbl, r))
        if r["blew"]:
            print(f"    BLEW UP @step {r['blew']} (t={r['tphys']:.1f} h)")
        else:
            print(f"    max|u|={r['max_u']:6.1f}  max|w|={r['max_w']:5.2f}  "
                  f"max(θ′/θ̄)={100*r['small']:5.2f}%  "
                  f"continuity={100*r['flux_ratio']:5.2f}%")

    print("\n" + "=" * 70)
    print("SUMMARY:")
    for lbl, r in rows:
        v = "BLEW" if r["blew"] else f"{r['max_u']:6.1f}"
        print(f"  {lbl:>34}  max|u| {v}")
    print("\nREAD (registered rules in ENVELOPE_INTENSIFICATION.md):")
    print("  nz64-normalized ≥65 and ≥15 above nz64-historical")
    print("    → drag discretization CONFIRMED as the spin-down mechanism;")
    print("      LH82 caveat-3 flag closes.")
    print("  nz64-normalized ≈ nz64-historical → drag factor insufficient;")
    print("    next suspects: vertical advection / heating-layer resolution.")
    print(f"\nWall time: {time.time()-t0:.0f}s")
