"""
dz-sensitivity of heated intensification — the last open flag (2026-07).

With the drag artifact removed (column_normalized on both grids), the heated
Phase-3B configuration still reads max|u| 69.6 at nz=32 vs 54.7 at nz=64 — a
~15 m/s genuinely resolution-driven gap, riding on the secondary circulation
(max|w| dropped ~3× at nz=64 in the LH82 study).  Working hypothesis, stated
before running: **the anomaly may be the COARSE grid's** — discrete vertical
operators under-estimate the wavenumbers of marginally-resolved heating
(width_z = 3 km ≈ 4.8 cells at dz=625), so the anelastic response meets less
opposition and the coarse grid OVER-produces w.  If so, nz≥64 is the
convergent regime and the "nz=64 spin-down" was never a spin-down at all.

ROWS (all Q=1e-2, ε=0, dt=15, upwind5h, buoyancy ON, drag NORMALIZED —
the one known dz-artifact held fixed):

    R1  nz=32              (drag-matched ref: 69.6)
    R2  nz=64              (drag-matched ref: 54.7)
    R3  nz=96              ← THE row: convergence direction
    R4  nz=64, heating sampled at nz=32 centers (z_sample_nz=32)
                           ← forcing-representation control

Registered predictions: DZ_HEATED_SENSITIVITY.md (frozen before any GPU run).
A blow-up at nz=96 is a datum (vertical CFL headroom shrinks); fall back to
$env:DZH_DT=10 and note it.

Usage:  $env:LH82_STEPS = 1000 ; python -m oracle_v8.measure_dz_heated
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


def run(g, Q_max, dt, nsteps, z_sample_nz=None):
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
        nx=g.nx, ny=g.ny, nz=g.nz, Lx=LX, Ly=LY, Lz=LZ,
        z_sample_nz=z_sample_nz)
    cfg.divergence_damping = None                     # ε=0, the Phase-3B recipe
    cfg.surface_drag = SurfaceDragComponent(          # drag-matched everywhere
        Cd=CD, H_bl=H_BL, u_env=0.0, v_env=0.0, column_normalized=True)
    integ = RK3Integrator(config=cfg, base=g.base)

    heartbeat = max(nsteps // 4, 1)
    for n in range(1, nsteps + 1):
        state, _ = integ.step(state, dt=dt, step_number=n)
        if mx(state.u) > 400.0 or not np.isfinite(mx(state.u)):
            return dict(blew=n, tphys=n * dt / 3600.0)
        if n % heartbeat == 0 and n < nsteps:
            print(f"      … step {n}/{nsteps}  max|u|={mx(state.u):6.1f}  "
                  f"max|w|={mx(state.w):5.2f}", flush=True)
    small, flux_ratio = g.neglected_diag(state)
    return dict(blew=None, max_u=mx(state.u), max_w=mx(state.w),
                small=small, flux_ratio=flux_ratio)


if __name__ == "__main__":
    N30 = int(os.environ.get("LH82_STEPS", "1000"))
    N15 = 2 * N30
    DT = float(os.environ.get("DZH_DT", "15"))
    nst = int(round(N15 * 15.0 / DT))                 # same physical time
    t0 = time.time()
    print(f"dz-sensitivity of heated intensification  xp={xp.__name__}  "
          f"dt={DT:.0f}  steps={nst}  (drag NORMALIZED everywhere)")
    print("  drag-matched refs: nz=32 → 69.6 ; nz=64 → 54.7 (gap 14.9)")

    rows = []
    for lbl, dims, zs in (
            ("R1 nz=32              (ref 69.6)", (128, 128, 32), None),
            ("R2 nz=64              (ref 54.7)", (128, 128, 64), None),
            ("R3 nz=96  ← convergence direction", (128, 128, 96), None),
            ("R4 nz=64, heating@nz32 (forcing ctl)", (128, 128, 64), 32)):
        print(f"\n  {lbl}:")
        r = run(Grid(*dims), 1.0e-2, DT, nst, z_sample_nz=zs)
        rows.append((lbl, r))
        if r["blew"]:
            print(f"    BLEW UP @step {r['blew']} (t={r['tphys']:.1f} h) — "
                  f"datum; retry with DZH_DT=10")
        else:
            print(f"    max|u|={r['max_u']:6.1f}  max|w|={r['max_w']:5.2f}  "
                  f"max(θ′/θ̄)={100*r['small']:5.2f}%  "
                  f"continuity={100*r['flux_ratio']:5.2f}%")

    print("\n" + "=" * 70)
    print("SUMMARY:")
    for lbl, r in rows:
        v = "BLEW" if r["blew"] else (f"{r['max_u']:6.1f}   max|w| "
                                      f"{r['max_w']:5.2f}")
        print(f"  {lbl:>40}  max|u| {v}")
    print("\nREAD (registered rules in DZ_HEATED_SENSITIVITY.md):")
    print("  P-D1: 32→64→96 monotone with SHRINKING increment (|Δ(96−64)| ≤")
    print("        0.5·|Δ(64−32)|) → coarse over-response; fine side convergent.")
    print("  P-D2: R4 ≈ R2 (within ~2) → forcing representation NOT the driver.")
    print("  P-D3: max|w| ordering tracks max|u| ordering.")
    print(f"\nWall time: {time.time()-t0:.0f}s")
