"""
Oracle V8 — f-plane translation harness (V8.4.1)
================================================

The clean discriminator the ensemble converged on, replacing the contaminated
test_translation.py (which was short, Hugo-like, tracker-dependent, and still
included BuoyancyComponent — Five).

LOGIC
-----
On an f-plane (f = f0 constant, β OFF) a balanced vortex embedded in a uniform
flow U must translate at exactly U — there are no β-gyres, no taper, nothing to
deflect it. So the meridional translation efficiency

        eff_y = (y_center(T) - y_center(0)) / (v_env · T)

is a pure measure of advection/component fidelity:

  eff_y ≈ 1 for all intensities      → advection is faithful; the Run 12 deficit
                                        REQUIRES β → it is the reversed-taper /
                                        β-gyre mechanism (Gemini). Go to the
                                        domain-enlargement run to confirm.
  eff_y < 1, and falls with Vmax      → the intense vortex under-translates even
                                        with no β → numerical / component-
                                        interaction lag (Copilot / Five). The
                                        ε-sweep and resolution rows then localize
                                        which component.

This isolates β (Gemini) from numerics (Copilot/Five) in a single campaign,
before touching beta-plane Katrina again.

NOTE: this imports the real solver, so run it in the repo (woe_env). It does not
run in the NumPy test rig. Center is measured with diagnostics.vorticity_center
(window seeded by dead reckoning) — independent of the θ′ tracker.
"""
from __future__ import annotations
import time
import numpy as np

from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.storm_data import KATRINA, IVAN
from oracle_v8.solver import (
    LH82AnelasticEquationSet,
    AdvectionComponent,
    CoriolisComponent,
    SurfaceDragComponent,
    IntensityCapComponent,
    HyperDiffusionComponent,
    HorizontalDiffusionComponent,
    HelmholtzDivergenceDampingComponent,
    NewtonianCoolingComponent,
    AnelasticProjection,
    OperatorConfig,
    RK3Integrator,
)
from oracle_v8.grid.staggering import LorenzStaggering
from oracle_v8.backend import xp, to_numpy


# --- GPU requirement (V8.6.2) -------------------------------------------------
# Running outside woe_env makes oracle_v8.backend silently fall back to numpy
# (cupy ImportError) — the whole stack then runs 5-10x slower on CPU with no
# warning.  Third silent-fallback found this week (ERA5 constant steering,
# synthetic obs tracks, now this) — production runs must announce their backend
# and abort on the wrong one.  Set REQUIRE_GPU=False for a deliberate CPU run.
REQUIRE_GPU = True
if xp.__name__ != "cupy":
    if REQUIRE_GPU:
        raise SystemExit(
            "\nFATAL: oracle_v8.backend selected NUMPY (cupy import failed)."
            "\nYou are probably outside woe_env — activate it and re-run."
            "\nSet REQUIRE_GPU=False in this script for a deliberate CPU run."
        )
    print("[backend] WARNING: running on NUMPY (CPU) — deliberate CPU run")
else:
    try:
        import cupy as _cp
        _d = _cp.cuda.runtime.getDeviceProperties(0)["name"]
        print(f"[backend] cupy ACTIVE on "
              f"{_d.decode() if isinstance(_d, bytes) else _d}")
    except Exception:
        print("[backend] cupy ACTIVE")
from oracle_v8 import diagnostics as dg
from dataclasses import replace as dc_replace

# ---- fixed domain (vertical); horizontal set per run for the resolution row --
Lx_BASE = Ly = 2_000_000.0
Lz       = 20_000.0
nz       = 32
dz       = Lz / nz
NU4      = 3.0e11
RMAX_M   = 75_000.0
R_ENV_M  = 500_000.0
KATRINA_F_REF = float(KATRINA["f"])
IVAN_F_REF = float(IVAN["f"])
DT       = 30.0
N_PREBAL = 5

z_centers = (np.arange(nz) + 0.5) * dz
theta0_arr = 300.0 * np.exp(0.01**2 * z_centers / 9.81)
Pi = np.zeros(nz)
Pi[0] = 1.0 - (9.81 / 1004.5) * z_centers[0] / theta0_arr[0]
for k in range(nz - 1):
    dl = z_centers[k + 1] - z_centers[k]
    Pi[k + 1] = Pi[k] - (9.81 / 1004.5) * (dl / 2.0) * (
        1.0 / theta0_arr[k] + 1.0 / theta0_arr[k + 1])
p0_arr   = 100_000.0 * Pi ** (1004.5 / 287.04)
rho0_arr = p0_arr / (287.04 * theta0_arr * Pi)


class _Base:
    z = z_centers; rho0 = rho0_arr; theta0 = theta0_arr


def _subcell_refine(state, x_m, y_m, dx, patch=9):
    """Sub-cell refinement of a cell-resolution centre fix (V8.6.2).

    The trace experiment showed dg.vorticity_center returns LATTICE-SNAPPED
    positions with a flicker amplitude of ±1.5 cells (the Galilean control
    occupied exactly two grid values for 14 h) → eff noise ±0.09-0.13 at
    10-14 h horizons.  This refinement recomputes low-level cyclonic ζ with
    the storm_tracker conventions (roll-based central differences, lower
    nz//2 column mean, [1,2,1] smooth) and replaces the cell fix with the
    ζ²-weighted centroid of a (2·patch+1)² patch about it.  ζ² weighting
    locks onto the compact eye peak and resists the asymmetric-skirt bias
    that made plain CoM lag south-west in Katrina Run 9.  Precision ≈ 0.1
    cell ≈ 1.6 km → eff floor < ±0.01.
    """
    nzh  = state.u.shape[2] // 2
    v_lo = state.v[:, :, :nzh]
    u_lo = state.u[:, :, :nzh]
    dvdx = (xp.roll(v_lo, -1, axis=0) - xp.roll(v_lo, 1, axis=0)) / (2.0 * dx)
    dudy = (xp.roll(u_lo, -1, axis=1) - xp.roll(u_lo, 1, axis=1)) / (2.0 * dx)
    z = xp.mean(xp.maximum(dvdx - dudy, 0.0), axis=2)
    z = (xp.roll(z, 1, 0) + 2.0 * z + xp.roll(z, -1, 0)) * 0.25
    z = (xp.roll(z, 1, 1) + 2.0 * z + xp.roll(z, -1, 1)) * 0.25
    nx_, ny_ = z.shape
    # Core-scale ITERATIVE centroid (V8.6.3).  The first trace re-run showed
    # the ζ argmax hops between two near-degenerate peaks ±1 cell from true
    # centre (ring degeneracy / trochoidal core wobble); a small patch refines
    # around the WRONG question.  A wide patch (±9 cells ≈ ±140 km ≈ 2·Rmax)
    # centroids over the whole core — the circulation centre, hop-free (NHC
    # best tracks smooth trochoidal wobble for the same reason).  Two
    # iterations recentre the patch on its own answer, removing the
    # truncation bias of a patch seeded off-centre.
    xc, yc = x_m / dx - 0.5, y_m / dx - 0.5
    if not (np.isfinite(xc) and np.isfinite(yc)):
        return x_m, y_m              # blown-up/NaN field → fall back to seed
    for _ in range(2):
        i0, j0 = int(round(xc)), int(round(yc))
        sx = slice(max(i0 - patch, 0), min(i0 + patch + 1, nx_))
        sy = slice(max(j0 - patch, 0), min(j0 + patch + 1, ny_))
        w = z[sx, sy] ** 2
        tot = float(xp.sum(w))
        if not np.isfinite(tot) or tot < 1e-20:
            return x_m, y_m
        ixs = xp.arange(sx.start, sx.stop, dtype=xp.float64)[:, None]
        iys = xp.arange(sy.start, sy.stop, dtype=xp.float64)[None, :]
        xc = float(xp.sum(ixs * w) / tot)
        yc = float(xp.sum(iys * w) / tot)
    return (xc + 0.5) * dx, (yc + 0.5) * dx


def run_translation(Vmax, u_env=0.0, v_env=5.0, epsilon=0.5,
                    drag_on=True, nx=128, hours=10.0, v_cap=None, dom=None,
                    r_env=None, beta=False, wind_taper=False,
                    taper_start_frac=0.5, keep_theta=False, steer_ramp=None,
                    Rmax=None, B=None, f_ref=None, snapshot_hours=None,
                    diff_form="hyper", nu_H=2.0e5,
                    subcell=True, verbose=False):
    """One translation run. f-plane by default (beta=False, β OFF — eff≈1 ⇒
    advection faithful).  beta=True turns on the β-plane Coriolis (with the
    periodic taper, exactly as run_katrina) so that with u_env=v_env=0 the
    storm's motion IS its β-drift.  Returns eff_x/eff_y, the raw drift vector,
    and the Vmax track."""
    Lx = Ly = dom if dom is not None else Lx_BASE
    ny = nx
    dx = Lx / nx
    f0 = float(f_ref if f_ref is not None else KATRINA_F_REF)
    renv = r_env if r_env is not None else R_ENV_M
    rmax = Rmax if Rmax is not None else RMAX_M
    bb   = B if B is not None else KATRINA["B"]             # <-- NEW

    init = HollandVortexInit(Vmax=Vmax, Rmax=rmax, B=bb, f=f0,
                             R_env=renv, u_env=u_env, v_env=v_env,
                             wind_taper=wind_taper,
                             taper_start_frac=taper_start_frac)
    state = init.build_state(nx, ny, nz, Lx, Ly, _Base())
    if not keep_theta:
        # historical default: kill the (passive) warm core so eff is pure advection
        state = dc_replace(state, theta_prime=xp.zeros_like(state.theta_prime))

    # pre-balance (projection only) — identical to run_katrina
    prebal = RK3Integrator(
        config=OperatorConfig(
            equation_set=LH82AnelasticEquationSet(),
            staggering=LorenzStaggering(),
            projection=AnelasticProjection(nx=nx, ny=ny, nz=nz,
                                           Lx=Lx, Ly=Ly, Lz=Lz)),
        base=_Base())
    for i in range(N_PREBAL):
        state, _ = prebal.step(state, dt=1.0, step_number=i)

    comps = dict(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        advection=AdvectionComponent(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz),
        coriolis=(CoriolisComponent(f=f0, mode="beta_plane", Ly=Ly, ny=ny,
                                    u_env=u_env, v_env=v_env, periodic_taper=True)
                  if beta else
                  CoriolisComponent(f=f0, mode="f_plane",
                                    u_env=u_env, v_env=v_env)),  # ref → no inertial loop
        horiz_diffusion=(HorizontalDiffusionComponent(nu_H=nu_H, Lx=Lx, Ly=Ly,
                                                      nx=nx, ny=ny)
                         if diff_form == "laplacian" else
                         HyperDiffusionComponent(nu4=NU4, Lx=Lx, Ly=Ly,
                                                 nx=nx, ny=ny)),
        newtonian_cooling=NewtonianCoolingComponent(tau=1800.0),
        projection=AnelasticProjection(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz),
    )
    if epsilon > 0:
        comps["divergence_damping"] = HelmholtzDivergenceDampingComponent(
            epsilon=epsilon, Lx=Lx, Ly=Ly, nx=nx, ny=ny)
    if drag_on:
        comps["surface_drag"] = SurfaceDragComponent(
            Cd=1.5e-3, H_bl=1000.0, u_env=u_env, v_env=v_env)
    integrator = RK3Integrator(config=OperatorConfig(**comps), base=_Base())
    cap_comp = (IntensityCapComponent(v_cap=v_cap, tau=300.0,
                                      u_env=u_env, v_env=v_env)
                if v_cap is not None else None)

    cx0, cy0 = Lx / 2.0, Ly / 2.0          # init center (domain center)
    n_steps = int(hours * 3600.0 / DT)
    diag_every = 60
    # live background (constant unless steer_ramp drives the lockstep relaxation)
    u_envc, v_envc = u_env, v_env
    alpha_steer = min(1.0, (diag_every * DT) / 10800.0)   # TAU_STEER = 3 h
    gx, gy = cx0, cy0                       # dead-reckon seed, env-integrated
    env_ix = env_iy = 0.0                   # ∫ applied env dt (eff denominator)
    t_arr, xt, yt, vmt = [], [], [], []
    _snaps = []
    _snap_targets = sorted(float(h) for h in snapshot_hours) if snapshot_hours else []
    def _to_host(a):
        try:
            return a.get()          # cupy -> numpy (copy to host, GPU state untouched)
        except AttributeError:
            return np.asarray(a).copy()
    for n in range(n_steps + 1):
        if n > 0:
            state, _ = integrator.step(state, dt=DT, step_number=n - 1)
            if cap_comp is not None:
                _ct = cap_comp.compute_tendency(state, None, None, _Base(), DT)
                state = dc_replace(state, u=state.u + DT * _ct.du_dt,
                                   v=state.v + DT * _ct.dv_dt)
        if n % diag_every == 0:
            t = n * DT
            if n > 0:
                # advance the dead-reckon seed + env integral with the APPLIED env
                gx += u_envc * (diag_every * DT)
                gy += v_envc * (diag_every * DT)
                env_ix += u_envc * (diag_every * DT)
                env_iy += v_envc * (diag_every * DT)
                # lockstep relaxation toward the ramped target — mirrors
                # run_ivan.py exactly: shift the state background AND every
                # env-relative reference (coriolis/drag/cap) by the same Δ
                if steer_ramp is not None:
                    u0, v0, u1, v1 = steer_ramp
                    frac  = n / n_steps
                    u_tgt = u0 + (u1 - u0) * frac
                    v_tgt = v0 + (v1 - v0) * frac
                    du = (u_tgt - u_envc) * alpha_steer
                    dv = (v_tgt - v_envc) * alpha_steer
                    state = dc_replace(state, u=state.u + du, v=state.v + dv)
                    u_envc += du
                    v_envc += dv
                    comps["coriolis"].set_env(u_envc, v_envc)
                    if "surface_drag" in comps:
                        comps["surface_drag"].set_env(u_envc, v_envc)
                    if cap_comp is not None:
                        cap_comp.set_env(u_envc, v_envc)
            vc = dg.vorticity_center(state, gx, gy, dx, dx, base=_Base(),
                                     window_m=250e3)
            xv, yv = vc["xv_m"], vc["yv_m"]
            if subcell:                       # V8.6.2 — kill the lattice flicker
                xv, yv = _subcell_refine(state, xv, yv, dx)
            iv = dg.low_level_vmax(state, _Base(), xv, yv,
                                   dx, dx, u_env=u_envc, v_env=v_envc)
            t_arr.append(t); xt.append(xv); yt.append(yv)
            vmt.append(iv["vmax_lowlvl"])
            if _snap_targets and abs(t - _snap_targets[0] * 3600.0) < diag_every * DT * 0.5:
                _snaps.append(dict(t=t, cx=xv, cy=yv,
                                   u=_to_host(state.u), v=_to_host(state.v),
                                   w=_to_host(state.w),
                                   theta_prime=_to_host(state.theta_prime),
                                   projection_potential=_to_host(state.projection_potential)))
                _snap_targets.pop(0)
            if verbose:
                print(f"    t={t/3600:5.2f}h  c=({vc['xv_m']/1e3:.0f},"
                      f"{vc['yv_m']/1e3:.0f})km  Vmax'={iv['vmax_lowlvl']:.0f}")

    T = t_arr[-1]
    if steer_ramp is None:                  # bit-identical to all prior modes
        eff_x = (xt[-1] - xt[0]) / (u_env * T) if abs(u_env) > 1e-6 else float("nan")
        eff_y = (yt[-1] - yt[0]) / (v_env * T) if abs(v_env) > 1e-6 else float("nan")
    else:                                   # ramped: eff vs ∫ applied env dt
        eff_x = (xt[-1] - xt[0]) / env_ix if abs(env_ix) > 1e3 else float("nan")
        eff_y = (yt[-1] - yt[0]) / env_iy if abs(env_iy) > 1e3 else float("nan")
    return {"Vmax": Vmax, "u_env": u_env, "v_env": v_env, "eps": epsilon,
            "drag": drag_on, "nx": nx, "eff_x": eff_x, "eff_y": eff_y,
            "vmax_end": vmt[-1], "vmax_max": max(vmt),
            "f_ref": f0,
            "dx_km": dx / 1e3, "dom_km": Lx / 1e3, "r_env_km": renv / 1e3,
            "drift_x": (xt[-1] - xt[0]) / T, "drift_y": (yt[-1] - yt[0]) / T,
            "T_h": T / 3600.0,
            "track": (t_arr, xt, yt, vmt),   # full center/intensity history
            "snapshots": _snaps}             # + full-state host snapshots


def _row(d):
    fx = f"{d['eff_x']:+.2f}" if not np.isnan(d["eff_x"]) else "  -- "
    return (f"  Vmax={d['Vmax']:3.0f}  nx={d['nx']:3d}  dx={d['dx_km']:4.1f}km  "
            f"dom={d['dom_km']:4.0f}km  eps={d['eps']:.2f}  drag={'on ' if d['drag'] else 'off'}"
            f" | eff_y={d['eff_y']:+.2f}  eff_x={fx}  "
            f"| Vmax'_end={d['vmax_end']:5.1f}  peak={d['vmax_max']:5.1f}")


def main():
    t0 = time.time()
    print("=" * 78)
    print("f-PLANE TRANSLATION TEST  (β OFF — eff≈1 ⇒ advection faithful ⇒ "
          "Run 12 deficit needs β)")
    print("=" * 78)

    print("\n[1] Intensity sweep  (v_env=5, ε=0.5, drag on, nx=128):")
    for Vm in (60, 90, 120):
        print(_row(run_translation(Vm, u_env=0.0, v_env=5.0)))

    print("\n[2] ε-sweep  (Vmax=120, v_env=5, drag on, nx=128):")
    for eps in (0.0, 0.25, 0.5):
        print(_row(run_translation(120, u_env=0.0, v_env=5.0, epsilon=eps)))

    print("\n[3] drag off/on  (Vmax=120, v_env=5, ε=0.5, nx=128):")
    for dr in (False, True):
        print(_row(run_translation(120, u_env=0.0, v_env=5.0, drag_on=dr)))

    print("\n[4] resolution  (Vmax=120, v_env=5, ε=0.5, drag on) — Copilot's "
          "under-resolution check:")
    for nxr in (128, 256):
        print(_row(run_translation(120, u_env=0.0, v_env=5.0, nx=nxr)))

    print("\n[5] oblique background  (Vmax=120, U=(-2,5), ε=0.5, drag on):")
    print(_row(run_translation(120, u_env=-2.0, v_env=5.0)))

    print("\n[6] CAPPED high intensity  (Vmax=120 init, |V'|→70 cap, v_env=5) — "
          "closes the Finding-A gap the blowup denied us:")
    for Vm in (120, 150):
        print(_row(run_translation(Vm, u_env=0.0, v_env=5.0, v_cap=70.0)))
    print("   → if eff_y≈1 here, advection is faithful even at Katrina-equivalent")
    print("     intensity once the runaway is contained → deficit is purely β/taper.")

    print("\n[7] CAPPED resolution & domain isolation (|V'|→70, v_env=5, β OFF) — "
          "what drove the Run-14 overshoot?")
    print("    baseline [6] was nx=128, dx=15.6km, dom=2000km → eff_y≈1.13")
    print("    (a) finer dx, same domain (isolates resolution):")
    print(_row(run_translation(120, u_env=0.0, v_env=5.0, v_cap=70.0, nx=256)))
    print("    (b) Run-14 grid: same dx=15.6km, BIG domain (isolates domain size):")
    print(_row(run_translation(120, u_env=0.0, v_env=5.0, v_cap=70.0,
                               nx=256, dom=4_000_000.0)))
    print("   → (b)≈1.1 (matches [6]) ⇒ big domain is numerically faithful with β OFF,")
    print("     so Run-14's overshoot is β-drift physics → chase R_env.")
    print("   → (b)≫1.1 ⇒ the big domain over-translates even without β → a separate,")
    print("     numerical problem to fix first.  (a) flags any pure resolution effect.")

    print("\n" + "=" * 78)
    print("READ: if eff_y≈1 across [1] → advection faithful, deficit is β "
          "(→ enlarge-domain run).")
    print("      if eff_y<1 and falls with Vmax in [1] → numerical; [2]-[4] "
          "say which component.")
    print(f"Wall time: {time.time()-t0:.0f}s")


def discriminator_only():
    """Run ONLY the Run-14 overshoot discriminator (section [7] rows) — skips
    [1]-[6] so you don't re-run the expensive uncapped 256² blowup to reach it.
    Usage:  python run_translation_test.py 7
    """
    t0 = time.time()
    print("=" * 78)
    print("RUN-14 OVERSHOOT DISCRIMINATOR  (capped |V'|→70, v_env=5, β OFF)")
    print("=" * 78)
    print("  baseline [6]: nx=128, dx=15.6km, dom=2000km → eff_y≈1.13")
    print("\n  (a) finer dx, same 2000km domain (isolates resolution → dx=7.8km):")
    print(_row(run_translation(120, u_env=0.0, v_env=5.0, v_cap=70.0, nx=256)))
    print("\n  (b) Run-14 grid: 256² in a 4000km box (same dx=15.6km, isolates domain):")
    print(_row(run_translation(120, u_env=0.0, v_env=5.0, v_cap=70.0,
                               nx=256, dom=4_000_000.0)))
    print("\n  READ:")
    print("   (b) ≈ 1.13  (matches [6]) ⇒ big domain is numerically faithful with β")
    print("               OFF ⇒ Run-14's overshoot is β-drift physics → chase R_env.")
    print("   (b) ≫ 1.13  ⇒ the big grid over-translates even without β → a separate")
    print("               numerical problem to fix before R_env.")
    print("   (a) vs 1.13 ⇒ isolates any pure finer-dx resolution effect.")
    print(f"\nWall time: {time.time()-t0:.0f}s")


def _compass(deg):
    dirs = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
            "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    return dirs[int((deg + 11.25) // 22.5) % 16]


def betadrift_only():
    """β-drift calibration — the R_env lever for the Run-14 overshoot.

    β-plane Coriolis ON (with the periodic taper, as run_katrina), but ZERO
    background flow: with no steering the storm's motion IS its β self-
    propagation.  Intensity is capped to 70 so only R_env (vortex size) varies,
    isolating how size sets the poleward β-drift.  At 128²/10h the storm drifts
    only ~100 km, staying parked near the centre, far from the taper zones.
    Usage:  python run_translation_test.py beta
    """
    import math
    t0 = time.time()
    print("=" * 78)
    print("β-DRIFT CALIBRATION  (β-plane, u_env=v_env=0, |V'|→70 cap, nx=128, dx=15.6km)")
    print("  WIND TAPER ON → R_env now bounds the vortex (bare Holland is R_env-inert)")
    print("  no steering ⇒ motion = pure β self-propagation; sweep R_env")
    print("=" * 78)
    print(f"  {'R_env(km)':>9}  {'drift(m/s)':>10}  {'toward':>7}  "
          f"{'(u,v) m/s':>15}  {'Vmax_end':>8}")
    for renv in (300e3, 400e3, 500e3):
        d = run_translation(120, u_env=0.0, v_env=0.0, v_cap=70.0,
                            beta=True, r_env=renv, wind_taper=True)
        ux, vy = d["drift_x"], d["drift_y"]
        spd = math.hypot(ux, vy)
        ang = math.degrees(math.atan2(ux, vy)) % 360   # 0=N,90=E,180=S,270=W
        print(f"  {renv/1e3:9.0f}  {spd:10.2f}  {_compass(ang):>7}  "
              f"({ux:+5.2f},{vy:+5.2f})  {d['vmax_end']:8.1f}")
    print("\n  READ: with the taper ON the drift should now FALL with R_env (smaller")
    print("        vortex → weaker β-drift).  observed large-TC β-drift ≈ 1–3 m/s NW.")
    print("        Pick the R_env whose drift lands there → run big-domain Katrina at it.")
    print("        (bare-Holland baseline was 2.17 m/s N, identical for all R_env.)")
    print(f"\nWall time: {time.time()-t0:.0f}s")


def effx_only():
    """Zonal vs meridional translation-efficiency check (f-plane, β OFF, capped).

    Tests for a translation ANISOTROPY: does the model follow a westward flow as
    faithfully as a poleward one?  If eff_x < eff_y, NW-moving storms get rotated
    toward vertical (north fine, west short) — exactly Hugo's east-bias signature,
    and a deficit the steering port could NOT fix.  Backgrounds are 5 m/s so the
    storm drifts ~180 km ≈ 11 cells in 10 h — well clear of the quantization floor
    that muddied betadrift.  Usage:  python run_translation_test.py effx
    """
    t0 = time.time()
    print("=" * 78)
    print("ZONAL vs MERIDIONAL TRANSLATION  (f-plane β OFF, |V'|→70 cap, nx=128, dx=15.6km)")
    print("  eff = storm displacement / (background·T);  1.0 = faithful advection")
    print("=" * 78)
    print(f"  {'bg (u,v) m/s':>14}  {'eff_x':>6}  {'eff_y':>6}  {'Vmax_end':>8}  note")
    cases = [((0.0, 5.0),  "pure N  (baseline, cf [1]/[6])"),
             ((-5.0, 0.0), "pure W  ← the test"),
             ((5.0, 0.0),  "pure E  (W/E symmetry check)"),
             ((-5.0, 5.0), "NW 45°"),
             ((-4.0, 8.0), "Hugo-like DLM direction")]
    for (ue, ve), note in cases:
        d = run_translation(120, u_env=ue, v_env=ve, v_cap=70.0)
        fx = f"{d['eff_x']:+.2f}" if not np.isnan(d["eff_x"]) else "  --  "
        fy = f"{d['eff_y']:+.2f}" if not np.isnan(d["eff_y"]) else "  --  "
        print(f"  ({ue:+4.1f},{ve:+4.1f})    {fx:>6}  {fy:>6}  {d['vmax_end']:8.1f}  {note}")
    print("\n  READ:")
    print("   eff_x ≈ eff_y ≈ 1.1  → ISOTROPIC: Hugo's east bias is the frozen steering,")
    print("                          and the time-varying steering port should fix it.")
    print("   eff_x <  eff_y        → ZONAL under-translation: the model rotates NW storms")
    print("                          toward vertical → steering port only partly helps,")
    print("                          there's a real westward-translation deficit to report.")
    print("   pure-W vs pure-E      → any W/E asymmetry flags a directional numerical bias.")
    print(f"\nWall time: {time.time()-t0:.0f}s")


def _grow(label, d):
    fx = f"{d['eff_x']:+.3f}" if not np.isnan(d["eff_x"]) else "  --  "
    return (f"  {label:<26} nx={d['nx']:3d}  dom={d['dom_km']:4.0f}km  "
            f"dx={d['dx_km']:4.1f}km | eff_y={d['eff_y']:+.3f}  eff_x={fx} "
            f"| Vmax'_end={d['vmax_end']:5.1f}")


def gate_only():
    """IVAN-GRID TRANSLATION GATE (V8.6.1) — the experiment that gates all
    further propagation/taper work.

    Ivan run-2 out-translated its own sampled background poleward by a growing
    +1.4…+2.9 m/s on a 5.4–6.0 m/s background → implied eff_y ≈ 1.26–1.50 IF
    the excess is grid-borne.  This mode measures pure-translation fidelity on
    the EXACT Ivan grid vs the validated production grid, with the production
    physics (wind taper ON, cap present, Vmax≈64 — note: all HISTORICAL eff
    baselines were UNTAPERED at 128²/2000 km, so compare gate rows only to
    each other, not to the old 1.13).

    Phase 1 (always):       G0 = 256²/4000 km (Hugo/Katrina grid, tapered
                            baseline) and G1 = 320²/5000 km (Ivan grid).
    Phase 2 (auto if        G2 = 256²/5000 km and G3 = 320²/4000 km —
    |Δ| ≥ 0.08):            triangulates WHICH variable eff_y tracks:
                              tracks nx  → G1≈G3 ≠ G0≈G2  (spectral mode-count
                                           indexing suspect)
                              tracks dom → G1≈G2 ≠ G0≈G3  (domain-scaled term)
                              tracks dx  → monotonic G3(12.5) > G0/G1(15.6)
                                           > G2(19.5)  (known finer-dx trend)
    Wall estimate: phase 1 ≈ 50 min, +50 min if phase 2 triggers (320² row
    ≈ 30 min, 256² ≈ 20 min at 10 h sim).
    Usage:  python run_translation_test.py gate
    """
    t0 = time.time()
    print("=" * 78)
    print("IVAN-GRID TRANSLATION GATE  (f-plane β OFF, v_env=+5 N, taper ON, "
          "r_env=500km,")
    print("  Vmax=64, |V'|→70 cap present — the production Ivan physics)")
    print("  Ivan run-2 implied eff_y ≈ 1.26–1.50 if grid-borne; tracker "
          "floor ≈ ±0.04")
    print("=" * 78)

    common = dict(u_env=0.0, v_env=5.0, v_cap=70.0,
                  f_ref=IVAN_F_REF,
                  wind_taper=True, r_env=500e3)

    print("\n[phase 1] the gate:")
    g0 = run_translation(64, nx=256, dom=4_000_000.0, **common)
    print(_grow("G0 production baseline", g0))
    g1 = run_translation(64, nx=320, dom=5_000_000.0, **common)
    print(_grow("G1 IVAN grid", g1))
    delta = g1["eff_y"] - g0["eff_y"]
    print(f"\n  Δ eff_y (Ivan − baseline) = {delta:+.3f}")

    g2 = g3 = None
    if abs(delta) >= 0.08:
        print("\n[phase 2] Δ exceeds floor — triangulating the scaling "
              "variable:")
        g2 = run_translation(64, nx=256, dom=5_000_000.0, **common)
        print(_grow("G2 256² in 5000km", g2))
        g3 = run_translation(64, nx=320, dom=4_000_000.0, **common)
        print(_grow("G3 320² in 4000km", g3))
        by_nx  = abs(g1["eff_y"] - g3["eff_y"]) + abs(g0["eff_y"] - g2["eff_y"])
        by_dom = abs(g1["eff_y"] - g2["eff_y"]) + abs(g0["eff_y"] - g3["eff_y"])
        print(f"\n  within-nx spread  = {by_nx:.3f}   (small ⇒ eff tracks nx)")
        print(f"  within-dom spread = {by_dom:.3f}   (small ⇒ eff tracks "
              f"domain size)")
        print(f"  dx ordering: G3({g3['dx_km']:.1f}km)={g3['eff_y']:+.3f}  "
              f"G0/G1(15.6km)={g0['eff_y']:+.3f}/{g1['eff_y']:+.3f}  "
              f"G2({g2['dx_km']:.1f}km)={g2['eff_y']:+.3f}")
        print(f"  (monotonic rise toward finer dx ⇒ eff tracks dx — the "
              f"known trend)")

    print("\n" + "=" * 78)
    print("VERDICT:")
    if delta >= 0.12:
        print(f"  GRID ARTIFACT CONFIRMED DOMINANT (Δ={delta:+.3f} ≥ 0.12).")
        print("  Ivan's −8.3 h is substantially grid-borne; phase-2 pattern "
              "names the")
        print("  mechanism class.  Production options: fix the scaling, or "
              "pin Ivan to")
        print("  256²/4000 km (documented taper graze) as a bridge.")
    elif delta >= 0.05:
        print(f"  PARTIAL grid contribution (Δ={delta:+.3f}).  Some of Ivan's "
              f"excess is")
        print("  grid-borne but not all → expect BOTH a grid fix AND a "
              "storm-state term.")
    else:
        print(f"  GRID CLEAN at the Ivan config (Δ={delta:+.3f} < 0.05).")
        print("  A balanced vortex translates faithfully on 320²/5000 km → "
              "Ivan's excess")
        print("  is STORM-STATE physics (re-intensification / vortex-flow "
              "interaction).")
        print("  That is the deep branch — design the intensifying-vortex "
              "harness next.")
    print(f"\nWall time: {time.time()-t0:.0f}s")


def gate_drift_only():
    """GATE companion — the zonal-drift channel on production grids.

    (a) W/E asymmetry at the IVAN grid (Vmax=120, |V'|→70 — matches the
        historical effx rows for comparability): has the ~0.6 m/s eastward
        numerics drift changed magnitude on the new grid?
    (b) CAP-OFF pair (Vmax=55, no cap) at BOTH production grids — the
        long-pending numerics-drift isolation, now run where it matters.
    Wall estimate ≈ 1.9 h (two 320² + one 256² + one 320² rows at 10 h).
    Usage:  python run_translation_test.py gate-drift
    """
    t0 = time.time()
    print("=" * 78)
    print("GATE-DRIFT  (f-plane β OFF, taper ON, r_env=500km — production "
          "physics)")
    print("=" * 78)
    print("\n[a] W/E asymmetry at the Ivan grid (Vmax=120, |V'|→70):")
    for ue, note in ((-5.0, "pure W"), (+5.0, "pure E")):
        d = run_translation(120, u_env=ue, v_env=0.0, v_cap=70.0,
                            wind_taper=True, r_env=500e3,
                            nx=320, dom=5_000_000.0, f_ref=IVAN_F_REF)
        print(_grow(f"{note} @ Ivan grid", d))
    print("  READ vs history (untapered 128²/2000): W 0.87 / E 1.13 ⇒ "
          "+0.6 m/s east drift.")
    print("\n[b] CAP-OFF numerics isolation (Vmax=55, no cap):")
    for nxr, domr, lbl in ((256, 4_000_000.0, "Hugo/Katrina grid"),
                           (320, 5_000_000.0, "Ivan grid")):
        for ue, note in ((-5.0, "W"), (+5.0, "E")):
            d = run_translation(55, u_env=ue, v_env=0.0, v_cap=None,
                                wind_taper=True, r_env=500e3,
                                nx=nxr, dom=domr, f_ref=IVAN_F_REF)
            print(_grow(f"cap-off {note} @ {lbl}", d))
    print("  READ: asymmetry persisting with NO cap ⇒ drift is advection/"
          "projection/")
    print("        center-finder (report as limitation); asymmetry gone ⇒ "
          "cap implicated")
    print("        after all (then reconcile with the cap-free V8.6 storm "
          "evidence).")
    print(f"\nWall time: {time.time()-t0:.0f}s")


def gate_beta_only():
    """GATE-BETA (V8.7) — clean static beta-drift: MAGNITUDE + DIRECTION.

    beta-plane, zero background: motion = pure beta self-propagation.  Tests the
    two hypotheses the intensification ladder raised:
      (1) MIS-ORIENTATION.  J1 gave net drift (-0.32 W, +1.89 N) ~= heading 350
          (nearly due north) vs theory's NW band.  A too-poleward beta-gyre is a
          direct cause of the poleward landfall bias.
      (2) INTENSITY (in)DEPENDENCE.  Fiorino & Elsberry (1989): adiabatic beta-
          drift is set by the OUTER wind structure, ~independent of Vmax.  64 vs
          120 matching -> R_env/taper is the lever; diverging -> the model
          couples beta-drift to intensity (a model artifact).
    48 h so the gyre matures; mature drift read from a LATE window (30-48h) after
    the Holland-init adjustment settles, not the whole-run net.
    Usage:  python run_translation_test.py gate-beta
    """
    import math
    t0 = time.time()
    TH_SPD = (1.5, 2.5)        # m/s   theory magnitude band
    TH_HDG = (290.0, 335.0)    # deg   theory heading band (toward, NW)
    print("=" * 78)
    print("GATE-BETA  (beta-plane, u=v=0, taper ON, r_env=500km, Ivan grid "
          "320^2/5000km, 48h)")
    print(f"  theory: |drift| {TH_SPD[0]}-{TH_SPD[1]} m/s, toward "
          f"{TH_HDG[0]:.0f}-{TH_HDG[1]:.0f} deg (NW), ~intensity-independent")
    print("=" * 78)

    def _window(track, ta_h, tb_h):
        tt, xs, ys, vm = track
        ia = min(range(len(tt)), key=lambda k: abs(tt[k] - ta_h * 3600.0))
        ib = min(range(len(tt)), key=lambda k: abs(tt[k] - tb_h * 3600.0))
        dts = tt[ib] - tt[ia]
        return (xs[ib] - xs[ia]) / dts, (ys[ib] - ys[ia]) / dts

    results = []
    for vm, cap, lbl in ((64, 70.0, "Ivan-strength (64)"),
                         (120, 70.0, "capped-max (120->70)")):
        d = run_translation(vm, u_env=0.0, v_env=0.0, v_cap=cap, beta=True,
                            wind_taper=True, r_env=500e3, nx=320,
                            dom=5_000_000.0, hours=48.0, f_ref=IVAN_F_REF)
        track = d["track"]
        ux, vy = _window(track, 30.0, 48.0)          # MATURE drift
        spd = math.hypot(ux, vy)
        hdg = math.degrees(math.atan2(ux, vy)) % 360  # toward; 0=N, 90=E
        results.append((lbl, spd, hdg, ux, vy))
        print(f"\n  {lbl}:")
        print(f"    net (0-48h)      ({d['drift_x']:+5.2f},{d['drift_y']:+5.2f}) m/s")
        print(f"    MATURE (30-48h)  |{spd:.2f}| m/s  toward {hdg:5.1f} "
              f"({_compass(hdg)})   ({ux:+5.2f},{vy:+5.2f})   "
              f"Vmax_end={d['vmax_end']:.1f}")
        tt, xs, ys, vmt = track
        print(f"    {'t(h)':>5}  {'|drift|':>7}  {'toward':>6}  {'Vmax':>5}")
        step = max(1, len(tt) // 8)
        for k in range(step, len(tt), step):
            dts = tt[k] - tt[k - step]
            ix = (xs[k] - xs[k - step]) / dts
            iy = (ys[k] - ys[k - step]) / dts
            print(f"    {tt[k]/3600:5.1f}  {math.hypot(ix, iy):7.2f}  "
                  f"{math.degrees(math.atan2(ix, iy)) % 360:6.1f}  {vmt[k]:5.1f}")

    print("\n" + "=" * 78)
    print("READ:")
    lbl0, spd0, hdg0, ux0, vy0 = results[0]
    in_hdg = TH_HDG[0] <= hdg0 <= TH_HDG[1]
    in_spd = TH_SPD[0] <= spd0 <= TH_SPD[1]
    hdg_verdict = ("IN BAND" if in_hdg else
                   ("TOO POLEWARD" if (hdg0 > TH_HDG[1] or hdg0 < 90) else "OFF"))
    spd_verdict = ("IN BAND" if in_spd else
                   ("TOO STRONG" if spd0 > TH_SPD[1] else "TOO WEAK"))
    print(f"  direction: mature heading {hdg0:.0f} vs {TH_HDG[0]:.0f}-"
          f"{TH_HDG[1]:.0f}  ->  {hdg_verdict}")
    print(f"  magnitude: {spd0:.2f} m/s vs {TH_SPD[0]}-{TH_SPD[1]}  ->  "
          f"{spd_verdict}")
    ds = abs(results[1][1] - results[0][1])
    dh = abs(results[1][2] - results[0][2])
    print(f"  intensity coupling (64 vs 120): d|drift|={ds:.2f} m/s, "
          f"dheading={dh:.0f} deg")
    print("     small -> outer-structure-controlled (R_env is the lever); "
          "large -> model couples beta to intensity")
    print(f"  decomposition: plug MATURE drift_y={vy0:+.2f} into "
          f"eff_y(J2) ~= 1 + drift_y/v_env_bar")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_renv():
    """GATE-BETA R_env SWEEP (V8.7) — is the beta-drift a controlled function of
    the outer-wind bound?

    gate-beta: Ivan-strength beta-drift is TOO STRONG (2.69), TOO POLEWARD
    (rotates 324->357 over 48h), NON-STATIONARY (still climbing at 48h), and it
    strengthens as the vortex DECAYS/SPREADS -> outer circulation, not Vmax.
    Sweep R_env at fixed taper_start_frac to test: tighter bound -> magnitude
    into 1.5-2.5, heading back to NW, and a PLATEAU.  Single intensity (64); 48h.
    Usage:  python run_translation_test.py gate-beta-renv
    """
    import math
    t0 = time.time()
    TH_SPD = (1.5, 2.5)
    TH_HDG = (290.0, 335.0)
    RENV_KM = (300, 400, 500, 650, 800)      # sweep; 500 = current anchor
    print("=" * 78)
    print("GATE-BETA R_env SWEEP  (beta-plane, u=v=0, Vmax=64, cap 70, "
          "taper_start_frac=0.5, 48h)")
    print(f"  theory: |drift| {TH_SPD[0]}-{TH_SPD[1]} m/s, toward "
          f"{TH_HDG[0]:.0f}-{TH_HDG[1]:.0f} (NW); want mag IN + hdg IN + plateau")
    print("=" * 78)

    def _win(track, ta_h, tb_h):
        tt, xs, ys, vm = track
        ia = min(range(len(tt)), key=lambda k: abs(tt[k] - ta_h * 3600.0))
        ib = min(range(len(tt)), key=lambda k: abs(tt[k] - tb_h * 3600.0))
        dts = tt[ib] - tt[ia]
        ux = (xs[ib] - xs[ia]) / dts
        vy = (ys[ib] - ys[ia]) / dts
        return math.hypot(ux, vy), math.degrees(math.atan2(ux, vy)) % 360, ux, vy

    rows = []
    for rkm in RENV_KM:
        d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=70.0, beta=True,
                            wind_taper=True, r_env=rkm * 1e3, nx=320,
                            dom=5_000_000.0, hours=48.0, f_ref=IVAN_F_REF)
        track = d["track"]
        spd, hdg, ux, vy = _win(track, 30.0, 48.0)         # mature drift
        s_early = _win(track, 30.0, 39.0)[0]
        s_late = _win(track, 39.0, 48.0)[0]
        plat = s_late - s_early                            # ~0 => plateaued
        mag_v = ("IN" if TH_SPD[0] <= spd <= TH_SPD[1]
                 else ("STRONG" if spd > TH_SPD[1] else "WEAK"))
        hdg_v = ("IN" if TH_HDG[0] <= hdg <= TH_HDG[1]
                 else ("POLEWARD" if (hdg > TH_HDG[1] or hdg < 90) else "OFF"))
        plat_v = "YES" if abs(plat) < 0.30 else "NO"
        rows.append((rkm, spd, hdg, mag_v, hdg_v, plat_v, vy, d["vmax_end"]))
        tt, xs, ys, vmt = track
        pts = []
        for th in (12, 24, 36, 48):
            i = min(range(len(tt)), key=lambda k: abs(tt[k] - th * 3600.0))
            j = max(0, i - max(1, len(tt) // 8))
            dts = tt[i] - tt[j]
            s = math.hypot((xs[i] - xs[j]) / dts, (ys[i] - ys[j]) / dts)
            h = math.degrees(math.atan2((xs[i] - xs[j]) / dts,
                                        (ys[i] - ys[j]) / dts)) % 360
            pts.append(f"t{th}:{s:.2f}@{h:.0f}")
        print(f"\n  R_env={rkm}km:  " + "   ".join(pts))
        print(f"    MATURE(30-48) |{spd:.2f}| toward {hdg:.0f} "
              f"({_compass(hdg)})  plateauD={plat:+.2f}  "
              f"Vmax_end={d['vmax_end']:.1f}")
        print(f"    verdict: mag {mag_v}   hdg {hdg_v}   plateau {plat_v}")

    print("\n" + "=" * 78)
    print("SUMMARY:")
    print(f"  {'R_env':>6}  {'|drift|':>7}  {'toward':>6}  {'mag':>6}  "
          f"{'hdg':>8}  {'plateau':>7}")
    for rkm, spd, hdg, mv, hv, pv, vy, vmend in rows:
        print(f"  {rkm:5d}k  {spd:7.2f}  {hdg:6.0f}  {mv:>6}  {hv:>8}  {pv:>7}")
    hits = [r for r in rows if r[3] == "IN" and r[4] == "IN"]
    if hits:
        best = min(hits, key=lambda r: abs(r[1] - 2.0))
        print(f"\n  SWEET SPOT: R_env={best[0]}km lands mag+hdg in band "
              f"(|{best[1]:.2f}| @ {best[2]:.0f}); plug drift_y={best[6]:+.2f} "
              f"into eff_y ~= 1 + drift_y/v_env_bar")
    else:
        print("\n  No R_env lands BOTH bands.  If |drift| tracks R_env but the "
              "heading stays poleward regardless, the rotation is NOT a pure "
              "size effect -> check domain/sponge next (rerun 64 on a larger "
              "domain).")
    print(f"\nWall time: {time.time()-t0:.0f}s")


def ladder_only():
    """STORM-INGREDIENT LADDER (V8.6.1) — the deep-branch discriminator.

    The gate showed a balanced vortex translates faithfully on the Ivan grid
    (eff_y ≈ 1.04 both grids), yet storm-Ivan out-ran its own background by
    +1.4…+2.9 m/s.  Whatever causes that lives in the ingredients the storm
    runs have and the clean harness doesn't.  This mode adds them ONE AT A
    TIME at the exact Ivan grid; the rung where eff_y jumps names the
    mechanism.

      L0  f-plane, constant 5 m/s N background      (same-session anchor; ≈
          the gate G1 row)
      L1  + β-plane                                 (β-drift × steering
          interaction; pure superposition predicts eff_y ≈ L0 + drift_N/5)
      L2  + lockstep ramp                           (time-varying steering
          machinery itself under test: targets ramp u −1.1→+0.6,
          v 3.9→6.6 — Ivan's DLM history compressed ~2.5× into 12 h, which
          AMPLIFIES any ramp-rate-dependent error; eff is measured against
          the ∫ of the APPLIED background, exactly the storm diagnostic)
      L3  + retained warm core (keep_theta)         (θ′ + Newtonian cooling
          present, as in every storm run)

    All rungs: 320²/5000 km, taper ON, r_env=500 km, Vmax=64, |V'|→70 cap,
    drag ON, 12 h.  Wall ≈ 37 min/rung ≈ 2.5 h total.
    READ:  Δeff_y ≈ +0.28 corresponds to +1.4 m/s on a 5 m/s background —
    the LOW end of storm-Ivan's excess.  A rung contributing ≥ +0.10 beyond
    superposition is implicated; if all rungs stay near superposition the
    remaining suspect is the re-intensification coupling (→ gate-beta's
    intensity contrast, then an intensifying-vortex harness).
    Usage:  python run_translation_test.py ladder
    """
    t0 = time.time()
    print("=" * 78)
    print("STORM-INGREDIENT LADDER  (Ivan grid 320²/5000km, taper ON, "
          "r_env=500km,")
    print("  Vmax=64, |V'|→70, drag ON, 12 h — adds storm ingredients one "
          "rung at a time)")
    print("=" * 78)
    common = dict(nx=320, dom=5_000_000.0, v_cap=70.0,
                  f_ref=IVAN_F_REF,
                  wind_taper=True, r_env=500e3, hours=12.0)
    ramp = (-1.1, 3.9, +0.6, 6.6)           # Ivan DLM history, compressed

    print("\n  L0  f-plane, constant background:")
    l0 = run_translation(64, u_env=0.0, v_env=5.0, **common)
    print(_grow("L0 anchor", l0))

    print("\n  L1  + β-plane (β × steering interaction):")
    l1 = run_translation(64, u_env=0.0, v_env=5.0, beta=True, **common)
    print(_grow("L1 beta", l1))

    print("\n  L2  + lockstep ramp (time-varying steering machinery):")
    l2 = run_translation(64, u_env=ramp[0], v_env=ramp[1], beta=True,
                         steer_ramp=ramp, **common)
    print(_grow("L2 beta+ramp", l2))

    print("\n  L3  + retained warm core (θ′ + Newtonian cooling live):")
    l3 = run_translation(64, u_env=ramp[0], v_env=ramp[1], beta=True,
                         steer_ramp=ramp, keep_theta=True, **common)
    print(_grow("L3 beta+ramp+theta", l3))

    print("\n" + "=" * 78)
    print("LADDER READ  (mean applied v ≈ 5 m/s ⇒ Δeff_y 0.10 ≈ +0.5 m/s "
          "excess):")
    print(f"  L1 − L0 = {l1['eff_y']-l0['eff_y']:+.3f}   β interaction "
          f"(superposition ≈ +drift_N/5 ≈ +0.12…+0.20)")
    print(f"  L2 − L1 = {l2['eff_y']-l1['eff_y']:+.3f}   lockstep-ramp "
          f"machinery (should be ≈ 0 if the shift is clean)")
    print(f"  L3 − L2 = {l3['eff_y']-l2['eff_y']:+.3f}   warm-core/"
          f"Newtonian-cooling coupling (should be ≈ 0; θ′ is passive)")
    print("  Any rung ≥ +0.10 beyond its expectation → that ingredient is "
          "implicated.")
    print("  All rungs clean → re-intensification coupling is the last "
          "suspect standing")
    print("  (run gate-beta, then design the intensifying-vortex harness).")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_domain():
    """GATE-BETA DOMAIN CHECK (V8.7) — is the poleward ROTATION the reversed-beta
    taper, or intrinsic?

    R_env sweep: |drift| is R_env-controlled but heading is poleward at every
    R_env, arriving as a TIME rotation (NW at t12 -> ~350 by t48).  The harness
    beta-plane uses periodic_taper=True: the beta-deviation is damped to f0 in the
    outer 20% of ny at each y-edge.  On 5000km/320 the north taper starts at
    y=4000km; the 650km vortex's gyre N lobe reaches it as it drifts north ->
    asymmetric beta -> poleward rotation.  Hold dx=15.625km, enlarge the domain so
    the clean interior (0.6*Ly) grows.  Rotation shrinks/delays with domain ->
    taper artifact (re-measure beta-drift clean; mature -> NW early value ->
    R_env-only is the full fix).  Persists -> intrinsic -> taper_start_frac next.
    Usage:  python run_translation_test.py gate-beta-domain
    """
    import math
    t0 = time.time()
    DOMS = ((5000, 320), (6500, 416), (8000, 512))    # dx = 15.625 km fixed
    print("=" * 78)
    print("GATE-BETA DOMAIN CHECK  (beta-plane, u=v=0, Vmax=64, R_env=650km, "
          "cap 70, dx=15.625km, 48h)")
    print("  reversed-beta taper = outer 20% of ny each side; clean interior = "
          "0.6*Ly")
    print("  rotation shrinks/delays with domain -> taper artifact; unchanged -> "
          "intrinsic")
    print("=" * 78)

    def _win(track, ta_h, tb_h):
        tt, xs, ys, vm = track
        ia = min(range(len(tt)), key=lambda k: abs(tt[k] - ta_h * 3600.0))
        ib = min(range(len(tt)), key=lambda k: abs(tt[k] - tb_h * 3600.0))
        dts = tt[ib] - tt[ia]
        ux = (xs[ib] - xs[ia]) / dts
        vy = (ys[ib] - ys[ia]) / dts
        return math.hypot(ux, vy), math.degrees(math.atan2(ux, vy)) % 360

    rows = []
    for dkm, nx in DOMS:
        d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=70.0, beta=True,
                            wind_taper=True, r_env=650e3, nx=nx,
                            dom=dkm * 1e3, hours=48.0)
        track = d["track"]
        tt, xs, ys, vmt = track
        e_spd, e_hdg = _win(track, 6.0, 18.0)          # early (pre-rotation)
        m_spd, m_hdg = _win(track, 30.0, 48.0)         # mature
        rot = ((m_hdg - e_hdg + 180) % 360) - 180      # signed; + = poleward
        taper_N_km = 0.80 * dkm                        # north taper onset
        ctr0 = dkm / 2.0
        y_late_km = ys[-1] / 1e3
        clear_km = taper_N_km - (y_late_km + 650.0)    # gyre-lobe proxy clearance
        rows.append((dkm, e_hdg, m_hdg, rot, m_spd, clear_km))
        print(f"\n  DOM={dkm}km/{nx}^2:  taperN@{taper_N_km:.0f}km  "
              f"vortex_y(0)@{ctr0:.0f}km")
        pts = []
        for th in (12, 24, 36, 48):
            i = min(range(len(tt)), key=lambda k: abs(tt[k] - th * 3600.0))
            j = max(0, i - max(1, len(tt) // 8))
            dts = tt[i] - tt[j]
            s = math.hypot((xs[i] - xs[j]) / dts, (ys[i] - ys[j]) / dts)
            h = math.degrees(math.atan2((xs[i] - xs[j]) / dts,
                                        (ys[i] - ys[j]) / dts)) % 360
            pts.append(f"t{th}:{s:.2f}@{h:.0f}")
        print("    " + "   ".join(pts))
        print(f"    early(6-18) {e_hdg:.0f} ({_compass(e_hdg)})  ->  "
              f"MATURE(30-48) |{m_spd:.2f}| {m_hdg:.0f} ({_compass(m_hdg)})  "
              f"rotation={rot:+.0f}deg")
        print(f"    vortex_y(48)={y_late_km:.0f}km  gyre-lobe clearance to "
              f"taperN ~ {clear_km:+.0f}km")

    print("\n" + "=" * 78)
    print("SUMMARY:")
    print(f"  {'DOM':>6}  {'early':>6}  {'mature':>7}  {'rotation':>9}  "
          f"{'|drift|':>7}  {'clear':>7}")
    for dkm, eh, mh, rot, ms, cl in rows:
        print(f"  {dkm:5d}k  {eh:6.0f}  {mh:7.0f}  {rot:+8.0f}  {ms:7.2f}  "
              f"{cl:+6.0f}k")
    rot0, rotN = abs(rows[0][3]), abs(rows[-1][3])
    print("\nREAD:")
    if rotN < rot0 - 8:
        print(f"  Rotation SHRINKS with domain ({rot0:.0f}->{rotN:.0f} deg) -> "
              f"REVERSED-BETA TAPER artifact.")
        print(f"  True beta-drift = clean-domain heading {rows[-1][1]:.0f} "
              f"({_compass(rows[-1][1])}).  If NW, R_env-only IS the full fix; "
              f"re-measure the R_env sweep on the largest domain.")
    else:
        print(f"  Rotation PERSISTS across domains ({rot0:.0f}->{rotN:.0f} deg) "
              f"-> INTRINSIC, not the taper.")
        print("  Next: taper_start_frac sweep (the wind-taper SHAPE sets the "
              "outer vorticity gradient that orients the gyres).")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_nu4():
    """GATE-BETA NU4 PROBE (V8.7) — is the ~8° poleward aim FLOOR diffusion-limited?

    The β-drift aim bottoms at ~343° (taper-start 200 km, mag in-band) and is
    outer-structure-invariant.  Hypothesis: hyperdiffusion (NU4) smears the β-gyre
    asymmetry that produces the WESTWARD ventilation, tilting the drift poleward.
    Sweep NU4 DOWN at the in-band config and watch the mature heading: rotates NW
    (343 → <335) while stable → DIFFUSION-limited (lower NU4 / less-dissipative
    scheme).  Flat until grid-noise/instability → under-RESOLVED → resolution next.
    Baseline NU4=3.0e11 (stability max ~3.1e12 at dx=15.6/dt=30 → room to drop).
    Usage:  python run_translation_test.py gate-beta-nu4
    """
    import math
    global NU4
    t0 = time.time()
    TH_HDG = (290.0, 335.0)
    NU4_SWEEP = (3.0e11, 1.0e11, 3.0e10, 1.0e10, 0.0)   # baseline → off
    _orig_nu4 = NU4
    print("=" * 78)
    print("GATE-BETA NU4 PROBE  (β-plane, u=v=0, Vmax=64, R_env=500 km, "
          "taper-start 200 km, cap 70, 5000km/320, 48h)")
    print(f"  floor to beat: mature ~343° (NW band {TH_HDG[0]:.0f}-{TH_HDG[1]:.0f}); "
          f"baseline NU4={_orig_nu4:.1e}, stability max ~3.1e12")
    print("  heading rotates NW as NU4 drops (while stable) → diffusion-limited; "
          "flat → resolution next")
    print("=" * 78)

    def _win(track, ta_h, tb_h):
        tt, xs, ys, vm = track
        ia = min(range(len(tt)), key=lambda k: abs(tt[k] - ta_h * 3600.0))
        ib = min(range(len(tt)), key=lambda k: abs(tt[k] - tb_h * 3600.0))
        dts = tt[ib] - tt[ia]
        return (math.hypot((xs[ib]-xs[ia])/dts, (ys[ib]-ys[ia])/dts),
                math.degrees(math.atan2((xs[ib]-xs[ia])/dts,
                                        (ys[ib]-ys[ia])/dts)) % 360)

    rows = []
    try:
        for nu4 in NU4_SWEEP:
            NU4 = nu4   # run_translation reads the module global at call time
            try:
                d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=70.0,
                                    beta=True, wind_taper=True,
                                    taper_start_frac=0.40, r_env=500e3,
                                    nx=320, dom=5_000_000.0, hours=48.0,
                                    f_ref=IVAN_F_REF)
            except Exception as e:
                rows.append((nu4, float("nan"), float("nan"), float("nan"),
                             False, "CRASH"))
                print(f"\n  NU4={nu4:.1e}: ⚠ run raised ({type(e).__name__}) — "
                      "unstable at this NU4.")
                continue
            track = d["track"]
            tt, xs, ys, vmt = track
            finite = all(math.isfinite(v)
                         for v in (xs[-1], ys[-1], d["vmax_end"]))
            m_spd, m_hdg = (_win(track, 30.0, 48.0) if finite
                            else (float("nan"), float("nan")))
            stable = finite and d["vmax_end"] < 120.0
            hdg_v = ("IN" if (finite and TH_HDG[0] <= m_hdg <= TH_HDG[1])
                     else ("POLEWARD" if finite else "UNSTABLE"))
            rows.append((nu4, m_hdg, m_spd, d["vmax_end"], stable, hdg_v))
            print(f"\n  NU4={nu4:.1e}:")
            if finite:
                pts = []
                for th in (12, 24, 36, 48):
                    i = min(range(len(tt)), key=lambda k: abs(tt[k]-th*3600.0))
                    j = max(0, i - max(1, len(tt)//8))
                    dts = tt[i]-tt[j]
                    s = math.hypot((xs[i]-xs[j])/dts, (ys[i]-ys[j])/dts)
                    h = math.degrees(math.atan2((xs[i]-xs[j])/dts,
                                                (ys[i]-ys[j])/dts)) % 360
                    pts.append(f"t{th}:{s:.2f}@{h:.0f}")
                print("    " + "   ".join(pts))
                print(f"    MATURE(30-48) |{m_spd:.2f}| {m_hdg:.0f} "
                      f"({_compass(m_hdg)})  Vmax_end={d['vmax_end']:.1f}  "
                      f"{'STABLE' if stable else '⚠UNSTABLE'}  hdg {hdg_v}")
            else:
                print("    ⚠ NON-FINITE — blew up at this NU4 (grid-scale noise "
                      "unchecked).")
    finally:
        NU4 = _orig_nu4   # ALWAYS restore the module global

    print("\n" + "=" * 78)
    print("SUMMARY:")
    print(f"  {'NU4':>9}  {'mature_hdg':>10}  {'|drift|':>7}  {'Vmax_end':>8}  "
          f"{'verdict':>9}")
    for nu4, mh, ms, ve, st, hv in rows:
        mhs = f"{mh:.0f}" if mh == mh else "nan"
        mss = f"{ms:.2f}" if ms == ms else "nan"
        ves = f"{ve:.1f}" if ve == ve else "nan"
        print(f"  {nu4:9.1e}  {mhs:>10}  {mss:>7}  {ves:>8}  {hv:>9}")
    stable_rows = [r for r in rows if r[4] and r[1] == r[1]]
    base_hdg = next((r[1] for r in rows if r[0] == _orig_nu4 and r[1] == r[1]),
                    float("nan"))
    print("\nREAD:")
    if stable_rows:
        best = min(stable_rows, key=lambda r: r[1])   # smallest hdg = most NW
        if base_hdg == base_hdg and best[1] < base_hdg - 6:
            print(f"  Heading ROTATES NW as NU4↓ ({base_hdg:.0f} → {best[1]:.0f} at "
                  f"NU4={best[0]:.1e}) → the aim floor is DIFFUSION-LIMITED.")
            print("  → a lower NU4 (or less-dissipative scheme) recovers the westward "
                  "ventilation; re-confirm the winning NU4 at storm scale, watch "
                  "stability over the full 52 h + full physics.")
        else:
            print(f"  Heading ~FLAT across NU4 (best {best[1]:.0f} vs baseline "
                  f"{base_hdg:.0f}) → NOT diffusion-limited.")
            print("  → the gyre is under-RESOLVED → resolution sweep next.")
    else:
        print("  No stable rows below baseline → NU4 can't drop without grid-noise "
              "blow-up at this dx → aim is RESOLUTION-bound → resolution sweep.")
    print(f"\nWall time: {time.time()-t0:.0f}s")


def trace_only():
    """L0 TRACE + GALILEAN CONTROL (V8.6.1) — tracker artifact vs
    frame-breaking acceleration.

    The ladder anchor (L0, 12 h) read eff_y=1.230 where the gate (identical
    config, 10 h) read 1.042 → implied ~10.8 m/s motion on a 5 m/s
    background in hours 10-12.  Two possible worlds:

      (a) REAL late-onset acceleration → frame-breaking numerics.  A uniform
          background on an f-plane is Galilean-removable (every component is
          env-referenced), so genuine excess can only come from a frame-
          dependent discretization error — advection dispersion is the
          candidate.  Signature: SMOOTH per-interval speed ramp in run A,
          quiet run B.
      (b) CENTER-FINDER ARTIFACT → the vorticity-center hops late in the
          run.  Signature: one or two discrete jumps in A's interval speeds;
          intervals otherwise ≈ 5 m/s.

    Run A: the L0 config, 14 h, full track printed.  The 10 h checkpoint
    must reproduce the gate's +1.042 — if it doesn't, the harness EDIT is
    implicated, not the physics (closing the old-code/new-code loop).
    Run B: ZERO background, otherwise identical, 14 h — the stationary
    vortex must not move; its track is the center-finder noise floor plus
    any spontaneous drift.

    Wall ≈ 90 min (two 14 h rows at 320²).
    Usage:  python run_translation_test.py trace
    """
    t0 = time.time()
    print("=" * 78)
    print("L0 TRACE + GALILEAN CONTROL  (Ivan grid 320²/5000km, taper ON, "
          "r_env=500km,")
    print("  Vmax=64, |V'|→70, drag ON, f-plane, 14 h)")
    print("=" * 78)
    common = dict(nx=320, dom=5_000_000.0, v_cap=70.0,
                  f_ref=IVAN_F_REF,
                  wind_taper=True, r_env=500e3, hours=14.0)

    def _trace(d, v_bg, label):
        tt, xs, ys, vm = d["track"]
        print(f"\n  {label} — per-interval center motion "
              f"(30-min intervals, hourly rows):")
        vmax_hdr = "Vmax'"
        print(f"  {'t(h)':>5}  {'x(km)':>8}  {'y(km)':>8}  "
              f"{'v_inst(m/s)':>11}  {'vy_inst':>8}  {'eff_y_cum':>9}  "
              f"{vmax_hdr:>6}")
        jumps = []
        for k in range(1, len(tt)):
            dt_s = tt[k] - tt[k - 1]
            vx = (xs[k] - xs[k - 1]) / dt_s
            vy = (ys[k] - ys[k - 1]) / dt_s
            sp = (vx * vx + vy * vy) ** 0.5
            jumps.append((sp, tt[k] / 3600.0, vx, vy))
            if k % 2 == 0:   # hourly print (diag every 30 min)
                eff = ((ys[k] - ys[0]) / (v_bg * tt[k])
                       if v_bg > 0 else float("nan"))
                ef = f"{eff:+9.3f}" if v_bg > 0 else "      -- "
                print(f"  {tt[k]/3600:5.1f}  {xs[k]/1e3:8.1f}  "
                      f"{ys[k]/1e3:8.1f}  {sp:11.2f}  {vy:+8.2f}  {ef}  "
                      f"{vm[k]:6.1f}")
        jumps.sort(reverse=True)
        print(f"  largest interval speeds: " + ",  ".join(
            f"{s:.1f} m/s @ t={th:.1f}h" for s, th, _, _ in jumps[:3]))
        if v_bg > 0:
            for hh in (8.0, 10.0, 12.0, 14.0):
                k = min(range(len(tt)), key=lambda i: abs(tt[i] - hh * 3600))
                eff = (ys[k] - ys[0]) / (v_bg * tt[k])
                tag = "  ← must match gate G1 (+1.042)" if hh == 10.0 else ""
                print(f"    eff_y({hh:4.1f}h) = {eff:+.3f}{tag}")
        return jumps

    print("\n[A] moving frame (v_env=+5 N):")
    da = run_translation(64, u_env=0.0, v_env=5.0, **common)
    ja = _trace(da, 5.0, "RUN A")

    print("\n[B] Galilean control (u_env=v_env=0 — the vortex must stand "
          "still):")
    db = run_translation(64, u_env=0.0, v_env=0.0, **common)
    jb = _trace(db, 0.0, "RUN B")
    ttb, xb, yb, _ = db["track"]
    net = ((xb[-1] - xb[0]) ** 2 + (yb[-1] - yb[0]) ** 2) ** 0.5
    print(f"  control net drift over 14 h: {net/1e3:.1f} km "
          f"({net/ttb[-1]:.2f} m/s)")

    print("\n" + "=" * 78)
    print("READ:")
    print("  A smooth speed ramp toward ~10 m/s + B quiet  → REAL frame-"
          "breaking")
    print("    acceleration (advection dispersion suspect) — the eff>1 "
          "family, the W/E")
    print("    asymmetry, and Ivan's −8.3 h are plausibly ONE mechanism.")
    print("  A shows 1-2 discrete hops, else ≈5 m/s        → CENTER-FINDER "
          "artifact —")
    print("    inspect the vorticity field at the hop time; eff numbers "
          "after the hop")
    print("    are invalid, ladder rungs need re-reading.")
    print("  B drifts too                                  → spontaneous "
          "numerical drift,")
    print("    background-independent — a third mechanism class.")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def intensify_ladder():
    """EMERGENT INTENSIFICATION LADDER (V8.7) — the storm-matched test."""
    t0 = time.time()
    V0 = IVAN["Vmax_ms"]
    Bi = IVAN["B"]

    print("=" * 78)
    print(f"INTENSIFICATION LADDER  (Ivan grid 320^2/5000km, taper ON, "
          f"r_env=500km, init Vmax={V0:.0f}, B={Bi:.2f}, 52h)")
    print("  let the barotropic spin-up EMERGE; read eff_y(t) vs Vmax(t)")
    print("=" * 78)

    common = dict(
        nx=320,
        dom=5_000_000.0,
        f_ref=IVAN_F_REF,
        v_cap=70.0,
        wind_taper=True,
        r_env=500e3,
        hours=52.0,
        Rmax=75_000.0,
        B=Bi,
    )
    ramp = (-1.1, 3.9, +0.6, 6.6)

    def _trace_intensify(d, v_bg, label):
        tt, xs, ys, vm = d["track"]
        print(f"\n  {label} - Vmax(t) and translation:")
        print(f"  {'t(h)':>5}  {'Vmax':>6}  {'y(km)':>8}  {'x(km)':>8}  "
              f"{'vy(m/s)':>8}  {'eff_y_cum':>9}")

        step = max(1, len(tt) // 13)
        for k in range(step, len(tt), step):
            dt_s = tt[k] - tt[k - 1]
            vy = (ys[k] - ys[k - 1]) / dt_s
            eff = ((ys[k] - ys[0]) / (v_bg * tt[k])
                   if v_bg > 0 else float("nan"))
            ef = f"{eff:+9.3f}" if v_bg > 0 else "      -- "
            print(f"  {tt[k]/3600:5.1f}  {vm[k]:6.1f}  {ys[k]/1e3:8.1f}  "
                  f"{xs[k]/1e3:8.1f}  {vy:+8.2f}  {ef}")

        print(f"  Vmax: init {vm[0]:.0f} -> end {vm[-1]:.0f} "
              f"(peak {max(vm):.0f})")

    print("\n  J0  f-plane, constant 5 m/s N bg  (does it self-intensify?):")
    j0 = run_translation(V0, u_env=0.0, v_env=5.0, **common)
    print(_grow("J0 fplane-52h", j0))
    _trace_intensify(j0, 5.0, "J0")

    print("\n  J1  beta-plane, zero bg  (beta self-propagation vs intensity):")
    j1 = run_translation(V0, u_env=0.0, v_env=0.0, beta=True, **common)
    print(f"  J1 end beta-drift ({j1['drift_x']:+.2f},{j1['drift_y']:+.2f}) m/s")
    _trace_intensify(j1, 0.0, "J1")

    print("\n  J2  beta-plane + steering ramp  (full Ivan-like config):")
    j2 = run_translation(
        V0,
        u_env=ramp[0],
        v_env=ramp[1],
        beta=True,
        steer_ramp=ramp,
        **common,
    )
    print(_grow("J2 beta+ramp-52h", j2))
    _trace_intensify(j2, 0.0, "J2")

    print("\n" + "=" * 78)
    print("READ:")
    print("  J0 intensifies AND eff_y rises with Vmax(t)  -> intensity-coupled")
    print("     NUMERICS, structural.")
    print("  J0 does NOT intensify (settles ~49)          -> Ivan's spin-up needs")
    print("     beta/steering -> intensification is beta/steering-coupled, and the")
    print("     over-translation is plausibly the SAME coupling (J1/J2 confirm).")
    print("  J1 beta-drift grows faster than gate-beta static beta(Vmax)  ->")
    print("     transient beta enhancement during active intensification.")
    print("  J2 reproduces Ivan's +2.9 m/s / 44->84       -> harness captures Ivan;")
    print("     ablate to isolate.  Doesn't -> real ERA5/track specifics required.")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def _circ_delta(h, ref):
    """signed h-ref in (-180,180]; NEGATIVE = counter-clockwise = westward (toward NW)."""
    return ((h - ref + 180.0) % 360.0) - 180.0

def _mature_drift(track, ta_h=30.0, tb_h=48.0):
    """Return (speed m/s, heading deg [0=N,CW], westward-component m/s) over [ta,tb]."""
    import math
    tt, xs, ys, vm = track
    ia = min(range(len(tt)), key=lambda k: abs(tt[k] - ta_h * 3600.0))
    ib = min(range(len(tt)), key=lambda k: abs(tt[k] - tb_h * 3600.0))
    dts = tt[ib] - tt[ia]
    spd = math.hypot((xs[ib]-xs[ia])/dts, (ys[ib]-ys[ia])/dts)
    hdg = math.degrees(math.atan2((xs[ib]-xs[ia])/dts, (ys[ib]-ys[ia])/dts)) % 360
    west = -spd * math.sin(math.radians(hdg))   # +ve = westward
    return spd, hdg, west

def gate_beta_diffform():
    """GATE-BETA DIFFUSION-FORM PROBE (V8.7, Arm A). Does the aim depend on the
    diffusion OPERATOR/strength?  NU4 couldn't be LOWERED (load-bearing); this tests
    operator form + the MORE-dissipation direction: ∇⁴ ref vs ∇² at 3 nu_H.  ∇² damps
    the resolved gyre far more than ∇⁴, so: aim moves (heading Δ, west-comp) → aim is
    dissipation-sensitive → try a SHARPER filter (∇⁶/spectral) for the less-dissipation
    fix.  aim FLAT while |drift| drops → decoupled from gyre amplitude → not diffusion.
    Requires the 3 harness edits above.  Usage: python run_translation_test.py gate-beta-diffform
    """
    import math
    t0 = time.time()
    TH_HDG = (290.0, 335.0); CAP = 70.0
    RUNS = [
        ("nabla4  NU4=3.0e11 (ref)",          dict(diff_form="hyper")),
        ("nabla2  nu_H=1.2e4 (grid-matched)", dict(diff_form="laplacian", nu_H=1.2e4)),
        ("nabla2  nu_H=5.0e4",                dict(diff_form="laplacian", nu_H=5.0e4)),
        ("nabla2  nu_H=2.0e5 (recommended)",  dict(diff_form="laplacian", nu_H=2.0e5)),
    ]
    print("=" * 78)
    print("GATE-BETA DIFFUSION-FORM PROBE (Arm A)  (β, u=v=0, Vmax=64, R_env=500, "
          "taper-start 200, cap 70, 320/5000km, 48h)")
    print(f"  ref ∇⁴ ~350°; NW band {TH_HDG[0]:.0f}-{TH_HDG[1]:.0f}. aim rotates WEST "
          "(heading↓, west-comp↑) → form/dissipation matters; flat → it doesn't")
    print("=" * 78)
    rows = []
    for label, kw in RUNS:
        try:
            d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=CAP, beta=True,
                                wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                                nx=320, dom=5_000_000.0, hours=48.0,
                                f_ref=IVAN_F_REF, **kw)
        except Exception as e:
            rows.append((label, float("nan"), float("nan"), float("nan"),
                         float("nan"), "CRASH"))
            print(f"\n  {label}: ⚠ raised ({type(e).__name__}) — unstable.")
            continue
        track = d["track"]; tt, xs, ys, vmt = track; ve = d["vmax_end"]
        finite = all(math.isfinite(v) for v in (xs[-1], ys[-1], ve))
        if finite:
            spd, hdg, west = _mature_drift(track)
        else:
            spd = hdg = west = float("nan")
        contam = finite and ve > CAP * 1.10
        verdict = ("CONTAM" if contam else ("OK" if finite else "UNSTABLE"))
        rows.append((label, hdg, spd, west, ve, verdict))
        print(f"\n  {label}:")
        if finite:
            print(f"    MATURE(30-48) |{spd:.2f}| hdg {hdg:.0f} ({_compass(hdg)})  "
                  f"WEST-comp {west:+.2f} m/s  Vmax_end={ve:.1f}  "
                  f"{'⚠CONTAM' if contam else 'clean'}")
        else:
            print("    ⚠ NON-FINITE.")
    print("\n" + "=" * 78); print("SUMMARY:")
    print(f"  {'form / coeff':<34}{'hdg':>5}{'Δvs∇⁴':>7}{'|drift|':>8}{'west':>7}"
          f"{'Vmax':>7}  verdict")
    ref_hdg = next((r[1] for r in rows if "ref" in r[0] and r[1] == r[1]), float("nan"))
    for label, hdg, spd, west, ve, vd in rows:
        hs = f"{hdg:.0f}" if hdg == hdg else "nan"
        dl = (f"{_circ_delta(hdg, ref_hdg):+.0f}"
              if (hdg == hdg and ref_hdg == ref_hdg) else "  —")
        ss = f"{spd:.2f}" if spd == spd else "nan"
        ws = f"{west:+.2f}" if west == west else "nan"
        vs = f"{ve:.0f}" if ve == ve else "nan"
        print(f"  {label:<34}{hs:>5}{dl:>7}{ss:>8}{ws:>7}{vs:>7}  {vd}")
    print("\nREAD:")
    lap = [r for r in rows if "nabla2" in r[0] and r[5] == "OK" and r[1] == r[1]]
    if lap and ref_hdg == ref_hdg:
        most = min(lap, key=lambda r: _circ_delta(r[1], ref_hdg))
        dl = _circ_delta(most[1], ref_hdg)
        if dl <= -6:
            print(f"  ∇² rotates aim WEST ({ref_hdg:.0f}→{most[1]:.0f}, {dl:+.0f}°) → aim IS "
                  "operator-sensitive. Surprising (∇² damps the gyre MORE) but actionable: "
                  "build a SHARPER filter (∇⁶/spectral) to push the less-dissipation fix.")
        elif dl >= 6:
            print(f"  ∇² rotates aim EAST/poleward ({ref_hdg:.0f}→{most[1]:.0f}, {dl:+.0f}°) → "
                  "more gyre dissipation worsens the aim → westward ventilation lives in the "
                  "resolved gyre; ∇⁴ is the less-damaging choice and the fix direction is LESS "
                  "dissipation → build a SHARPER filter (∇⁶/spectral) next.")
        else:
            print(f"  aim FLAT (Δ {dl:+.0f}° ≤ 6°) even as ∇² damps the gyre (check |drift|/west) "
                  "→ aim DECOUPLED from gyre dissipation → not a diffusion artifact → seat is the "
                  "Coriolis/projection path (see Arm B, then true-β/projection probe).")
    else:
        print("  ∇² runs unstable/contaminated — inspect Vmax; grid-matched ∇² may need a "
              "higher nu_H floor to stay clean before the comparison is meaningful.")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_divdamp():
    """GATE-BETA DIVERGENCE-DAMPER PROBE (V8.7, Arm B). Does the aim depend on the
    Helmholtz divergence damper?  Reviewer hypothesis: it suppresses the divergent
    flow that communicates the β-effect (geostrophic adjustment / Rossby radiation) →
    too-poleward aim.  Sweep epsilon 0.5 (baseline) → 0 (off; harness drops the
    component).  aim rotates WEST as eps↓ → damper WAS suppressing the β-communicating
    divergent flow (reduce it / gentler form = fix).  flat → not the driver.  Unstable
    at low eps → damper load-bearing → next test is measuring the vorticity the discrete
    Helmholtz solve injects (reviewer's 'not exact on a grid').  NO harness edit needed.
    Usage:  python run_translation_test.py gate-beta-divdamp
    """
    import math
    t0 = time.time()
    TH_HDG = (290.0, 335.0); CAP = 70.0
    EPS_SWEEP = (0.5, 0.25, 0.1, 0.0)   # baseline → off
    print("=" * 78)
    print("GATE-BETA DIVERGENCE-DAMPER PROBE (Arm B)  (β, u=v=0, Vmax=64, R_env=500, "
          "taper-start 200, cap 70, 320/5000km, 48h)")
    print(f"  baseline eps=0.5 ~350°; NW band {TH_HDG[0]:.0f}-{TH_HDG[1]:.0f}. aim rotates "
          "WEST as eps↓ → divergent flow carries the β-effect; flat → it doesn't")
    print("  eps=0 drops the damper entirely (harness: if epsilon > 0)")
    print("=" * 78)
    rows = []
    for eps in EPS_SWEEP:
        try:
            d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=CAP, beta=True,
                                wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                                nx=320, dom=5_000_000.0, hours=48.0,
                                epsilon=eps, f_ref=IVAN_F_REF)
        except Exception as e:
            rows.append((eps, float("nan"), float("nan"), float("nan"),
                         float("nan"), "CRASH"))
            print(f"\n  eps={eps:.2f}: ⚠ raised ({type(e).__name__}) — unstable.")
            continue
        track = d["track"]; tt, xs, ys, vmt = track; ve = d["vmax_end"]
        finite = all(math.isfinite(v) for v in (xs[-1], ys[-1], ve))
        if finite:
            spd, hdg, west = _mature_drift(track)
        else:
            spd = hdg = west = float("nan")
        contam = finite and ve > CAP * 1.10
        verdict = ("CONTAM" if contam else ("OK" if finite else "UNSTABLE"))
        rows.append((eps, hdg, spd, west, ve, verdict))
        print(f"\n  eps={eps:.2f} {'(damper OFF)' if eps == 0 else ''}:")
        if finite:
            print(f"    MATURE(30-48) |{spd:.2f}| hdg {hdg:.0f} ({_compass(hdg)})  "
                  f"WEST-comp {west:+.2f} m/s  Vmax_end={ve:.1f}  "
                  f"{'⚠CONTAM' if contam else 'clean'}")
        else:
            print("    ⚠ NON-FINITE.")
    print("\n" + "=" * 78); print("SUMMARY:")
    print(f"  {'eps':>5}{'hdg':>6}{'Δvs0.5':>8}{'|drift|':>8}{'west':>7}{'Vmax':>7}  verdict")
    ref_hdg = next((r[1] for r in rows if r[0] == 0.5 and r[1] == r[1]), float("nan"))
    for eps, hdg, spd, west, ve, vd in rows:
        hs = f"{hdg:.0f}" if hdg == hdg else "nan"
        dl = (f"{_circ_delta(hdg, ref_hdg):+.0f}"
              if (hdg == hdg and ref_hdg == ref_hdg) else "  —")
        ss = f"{spd:.2f}" if spd == spd else "nan"
        ws = f"{west:+.2f}" if west == west else "nan"
        vs = f"{ve:.0f}" if ve == ve else "nan"
        print(f"  {eps:>5.2f}{hs:>6}{dl:>8}{ss:>8}{ws:>7}{vs:>7}  {vd}")
    print("\nREAD:")
    clean = [r for r in rows if r[5] == "OK" and r[1] == r[1]]
    if clean and ref_hdg == ref_hdg:
        most = min(clean, key=lambda r: _circ_delta(r[1], ref_hdg))
        dl = _circ_delta(most[1], ref_hdg)
        if dl <= -6:
            print(f"  aim rotates WEST as eps↓ ({ref_hdg:.0f}→{most[1]:.0f} at eps={most[0]:.2f}, "
                  f"{dl:+.0f}°) → the damper WAS suppressing the β-communicating divergent flow. "
                  "Reduce eps / gentler form → re-confirm at storm scale.")
        else:
            print(f"  aim ~FLAT across eps ({ref_hdg:.0f}→{most[1]:.0f}, {dl:+.0f}°) → divergence "
                  "damping is NOT the aim driver.")
        bad = [r[0] for r in rows if r[5] in ("UNSTABLE", "CRASH", "CONTAM")]
        if bad:
            print(f"  NOTE low-eps rows unstable/contaminated {bad} → damper partly load-bearing; "
                  "if it can't be lowered cleanly, next test is measuring the vorticity the discrete "
                  "Helmholtz ψ-solve injects (∇×ψ-correction ≠ 0 on a discrete grid).")
    else:
        print("  No clean rows — damper load-bearing for stability; measure injected vorticity "
              "instead of reducing eps.")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_taper():
    """GATE-BETA TAPER-SHAPE SWEEP (V8.7) — does the wind-taper SHAPE reorient the
    beta-drift, or only resize it?

    Boundary exonerated (domain check: +18deg rotation domain-invariant); the
    beta-drift settles too poleward (~350) intrinsically.  The cosine wind-taper
    lays a negative-vorticity ring between taper_start_frac*R_env and R_env; that
    ring has its own beta-response on top of the main gyre.  Sweep where the taper
    starts (fixed R_env=650km), watch the MATURE heading:
      heading moves toward NW at some frac -> SHAPE is the DIRECTION knob;
      heading flat ~350 while only magnitude moves -> shape only resizes (like
        R_env), poleward bias is FUNDAMENTAL (gyre resolution / NU4) -> accept
        ~10deg and go to R34->storms, or probe NU4/resolution.
    Usage:  python run_translation_test.py gate-beta-taper
    """
    import math
    t0 = time.time()
    TH_SPD = (1.5, 2.5)
    TH_HDG = (290.0, 335.0)
    FRACS = (0.30, 0.50, 0.65, 0.80)        # 0.50 = current anchor
    print("=" * 78)
    print("GATE-BETA TAPER-SHAPE SWEEP  (beta-plane, u=v=0, Vmax=64, R_env=650km, "
          "cap 70, 5000km/320, 48h)")
    print(f"  taper runs [frac*R_env -> R_env]; theory hdg {TH_HDG[0]:.0f}-"
          f"{TH_HDG[1]:.0f} (NW), |drift| {TH_SPD[0]}-{TH_SPD[1]} m/s")
    print("  heading moves with frac -> SHAPE is the direction knob; flat ~350 "
          "-> bias is fundamental")
    print("=" * 78)

    def _win(track, ta_h, tb_h):
        tt, xs, ys, vm = track
        ia = min(range(len(tt)), key=lambda k: abs(tt[k] - ta_h * 3600.0))
        ib = min(range(len(tt)), key=lambda k: abs(tt[k] - tb_h * 3600.0))
        dts = tt[ib] - tt[ia]
        ux = (xs[ib] - xs[ia]) / dts
        vy = (ys[ib] - ys[ia]) / dts
        return math.hypot(ux, vy), math.degrees(math.atan2(ux, vy)) % 360

    rows = []
    for frac in FRACS:
        d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=70.0, beta=True,
                            wind_taper=True, taper_start_frac=frac, r_env=650e3,
                            nx=320, dom=5_000_000.0, hours=48.0,
                            f_ref=IVAN_F_REF)
        track = d["track"]
        tt, xs, ys, vmt = track
        e_spd, e_hdg = _win(track, 6.0, 18.0)
        m_spd, m_hdg = _win(track, 30.0, 48.0)
        rot = ((m_hdg - e_hdg + 180) % 360) - 180
        mag_v = ("IN" if TH_SPD[0] <= m_spd <= TH_SPD[1]
                 else ("STRONG" if m_spd > TH_SPD[1] else "WEAK"))
        hdg_v = ("IN" if TH_HDG[0] <= m_hdg <= TH_HDG[1]
                 else ("POLEWARD" if (m_hdg > TH_HDG[1] or m_hdg < 90) else "OFF"))
        taper0_km = frac * 650.0
        rows.append((frac, e_hdg, m_hdg, rot, m_spd, mag_v, hdg_v, taper0_km))
        print(f"\n  taper_start_frac={frac:.2f}  (taper {taper0_km:.0f}->650km):")
        pts = []
        for th in (12, 24, 36, 48):
            i = min(range(len(tt)), key=lambda k: abs(tt[k] - th * 3600.0))
            j = max(0, i - max(1, len(tt) // 8))
            dts = tt[i] - tt[j]
            s = math.hypot((xs[i] - xs[j]) / dts, (ys[i] - ys[j]) / dts)
            h = math.degrees(math.atan2((xs[i] - xs[j]) / dts,
                                        (ys[i] - ys[j]) / dts)) % 360
            pts.append(f"t{th}:{s:.2f}@{h:.0f}")
        print("    " + "   ".join(pts))
        print(f"    early(6-18) {e_hdg:.0f} -> MATURE(30-48) |{m_spd:.2f}| "
              f"{m_hdg:.0f} ({_compass(m_hdg)})  rot={rot:+.0f}  "
              f"Vmax_end={d['vmax_end']:.1f}")
        print(f"    verdict: mag {mag_v}   hdg {hdg_v}")

    print("\n" + "=" * 78)
    print("SUMMARY:")
    print(f"  {'frac':>5}  {'taper0':>7}  {'mature_hdg':>10}  {'|drift|':>7}  "
          f"{'mag':>6}  {'hdg':>8}")
    for frac, eh, mh, rot, ms, mv, hv, t0k in rows:
        print(f"  {frac:5.2f}  {t0k:6.0f}k  {mh:10.0f}  {ms:7.2f}  {mv:>6}  "
              f"{hv:>8}")
    hdg_span = max(r[2] for r in rows) - min(r[2] for r in rows)
    mag_span = max(r[4] for r in rows) - min(r[4] for r in rows)
    print("\nREAD:")
    print(f"  heading span across frac = {hdg_span:.0f} deg ;  "
          f"|drift| span = {mag_span:.2f} m/s")
    if hdg_span >= 12:
        best = min(rows, key=lambda r: abs(r[2] - 312.5))   # nearest NW-band center
        print(f"  SHAPE MOVES DIRECTION -> taper_start_frac is a direction knob. "
              f"Best aim: frac={best[0]:.2f} -> {best[2]:.0f} "
              f"({_compass(best[2])}); re-confirm at storm-relevant R_env (R34).")
    else:
        print("  Heading ~flat across frac -> shape only resizes, like R_env; the "
              "poleward bias is FUNDAMENTAL (gyre resolution / NU4 diffusion). "
              "Either accept ~10deg and go to R34->storm re-runs, or probe "
              "NU4/resolution next.")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_res():
    """GATE-BETA RESOLUTION SWEEP (V8.7) — is the ~8° aim floor under-RESOLVED or structural?

    The NU4 probe showed NU4 can't be lowered at dx=15.6 km without grid-noise
    blow-up (load-bearing for stability) → the floor is not diffusion-tunable.
    Hold NU4=3.0e11 + in-band physics, refine nx {320,480,640} on a fixed 5000 km
    domain (dx 15.6/10.4/7.8 km), watch the mature heading.  Rotates toward the NW
    band (signed delta NEGATIVE = westward) as dx shrinks, Vmax sane → under-RESOLVED
    (finer dx is the fix).  Flat → structural residual (paper: characterized ~8°).
    DT=30 fixed; advective CFL ~0.31 at nx=640, nu4 under the hyperdiff limit → stable.
    Usage:  python run_translation_test.py gate-beta-res
    """
    import math
    global NU4
    t0 = time.time()
    TH_HDG = (290.0, 335.0)
    CAP = 70.0
    NX_SWEEP = (320, 480, 640)
    DOM = 5_000_000.0
    _orig_nu4 = NU4
    NU4 = 3.0e11   # the only stable value; pin it explicitly

    def _circ_delta(h, ref):
        """signed h-ref in (-180,180]; NEGATIVE = counter-clockwise = westward (toward NW)."""
        return ((h - ref + 180.0) % 360.0) - 180.0

    def _win(track, ta_h, tb_h):
        tt, xs, ys, vm = track
        ia = min(range(len(tt)), key=lambda k: abs(tt[k] - ta_h * 3600.0))
        ib = min(range(len(tt)), key=lambda k: abs(tt[k] - tb_h * 3600.0))
        dts = tt[ib] - tt[ia]
        return (math.hypot((xs[ib]-xs[ia])/dts, (ys[ib]-ys[ia])/dts),
                math.degrees(math.atan2((xs[ib]-xs[ia])/dts,
                                        (ys[ib]-ys[ia])/dts)) % 360)

    print("=" * 78)
    print("GATE-BETA RESOLUTION SWEEP  (β-plane, u=v=0, Vmax=64, R_env=500 km, "
          "taper-start 200 km, cap 70, NU4=3.0e11, dom 5000 km, 48h)")
    print(f"  baseline (nx=320) ~350°; NW band {TH_HDG[0]:.0f}-{TH_HDG[1]:.0f}. "
          "heading drops (westward) as dx shrinks → under-resolved; flat → structural")
    print("  contamination flag: Vmax_end > cap+10% (~77) = grid-noise, reading suspect")
    print("=" * 78)

    rows = []
    try:
        for nx in NX_SWEEP:
            dx_km = DOM / nx / 1e3
            try:
                d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=CAP, beta=True,
                                    wind_taper=True, taper_start_frac=0.40,
                                    r_env=500e3, nx=nx, dom=DOM, hours=48.0,
                                    f_ref=IVAN_F_REF)
            except Exception as e:
                rows.append((nx, dx_km, float("nan"), float("nan"), float("nan"),
                             "CRASH"))
                print(f"\n  nx={nx} (dx={dx_km:.1f} km): ⚠ run raised "
                      f"({type(e).__name__}) — unstable at this resolution.")
                continue
            track = d["track"]
            tt, xs, ys, vmt = track
            ve = d["vmax_end"]
            finite = all(math.isfinite(v) for v in (xs[-1], ys[-1], ve))
            m_spd, m_hdg = (_win(track, 30.0, 48.0) if finite
                            else (float("nan"), float("nan")))
            contam = finite and ve > CAP * 1.10
            cfl = 80.0 * 30.0 / (DOM / nx)   # rough, Vmax~80 worst case
            verdict = ("CONTAM" if contam else
                       ("IN" if (finite and TH_HDG[0] <= m_hdg <= TH_HDG[1])
                        else ("POLEWARD" if finite else "UNSTABLE")))
            rows.append((nx, dx_km, m_hdg, m_spd, ve, verdict))
            print(f"\n  nx={nx}  dx={dx_km:.1f} km  (CFL~{cfl:.2f}):")
            if finite:
                pts = []
                for th in (12, 24, 36, 48):
                    i = min(range(len(tt)), key=lambda k: abs(tt[k]-th*3600.0))
                    j = max(0, i - max(1, len(tt)//8))
                    dts = tt[i]-tt[j]
                    s = math.hypot((xs[i]-xs[j])/dts, (ys[i]-ys[j])/dts)
                    h = math.degrees(math.atan2((xs[i]-xs[j])/dts,
                                                (ys[i]-ys[j])/dts)) % 360
                    pts.append(f"t{th}:{s:.2f}@{h:.0f}")
                print("    " + "   ".join(pts))
                print(f"    MATURE(30-48) |{m_spd:.2f}| {m_hdg:.0f} "
                      f"({_compass(m_hdg)})  Vmax_end={ve:.1f}  "
                      f"{'⚠CONTAM (grid noise — reading suspect)' if contam else 'clean'}")
            else:
                print(f"    ⚠ NON-FINITE — blew up at dx={dx_km:.1f} km "
                      "(scale NU4 ∝ dx⁴ if so).")
    finally:
        NU4 = _orig_nu4

    print("\n" + "=" * 78)
    print("SUMMARY:")
    print(f"  {'nx':>4}  {'dx_km':>6}  {'mature_hdg':>10}  {'Δ vs 320':>9}  "
          f"{'|drift|':>7}  {'Vmax_end':>8}  {'verdict':>9}")
    base = next((r for r in rows if r[0] == 320 and r[2] == r[2]), None)
    base_hdg = base[2] if base else float("nan")
    for nx, dxk, mh, ms, ve, vd in rows:
        mhs = f"{mh:.0f}" if mh == mh else "nan"
        dlt = (f"{_circ_delta(mh, base_hdg):+.0f}"
               if (mh == mh and base_hdg == base_hdg) else "  —")
        mss = f"{ms:.2f}" if ms == ms else "nan"
        ves = f"{ve:.1f}" if ve == ve else "nan"
        print(f"  {nx:>4}  {dxk:6.1f}  {mhs:>10}  {dlt:>9}  {mss:>7}  {ves:>8}  {vd:>9}")
    clean = [r for r in rows if r[5] not in ("CONTAM", "CRASH") and r[2] == r[2]]
    print("\nREAD:")
    if base and len(clean) >= 2:
        finest = min(clean, key=lambda r: r[1])           # smallest dx
        delta = _circ_delta(finest[2], base_hdg)          # negative = westward = good
        if delta <= -6.0 and finest[5] != "CONTAM":
            print(f"  Heading rotates WESTWARD as dx shrinks "
                  f"({base_hdg:.0f}° @320 → {finest[2]:.0f}° @{finest[0]} = {delta:+.0f}°) "
                  "→ gyre was UNDER-RESOLVED.")
            print(f"  → the ~8° floor is a resolution artifact; production wants dx≈"
                  f"{finest[1]:.0f} km. Re-confirm on Ivan full physics at that grid "
                  "(watch runtime + stability over 52 h).")
        elif abs(delta) < 6.0:
            print(f"  Heading ~FLAT across dx ({base_hdg:.0f}° → {finest[2]:.0f}°, "
                  f"{delta:+.0f}°) → the ~8° is STRUCTURAL to the β-gyre representation, "
                  "not resolution-curable.")
            print("  → paper line: 'characterized residual, ~8° poleward at the gyre "
                  "scale'; the storm cross-track E (~+125 km) is a known, bounded "
                  "systematic. Outer-structure + diffusion + resolution all exhausted.")
        else:
            print(f"  Heading moved EAST ({delta:+.0f}°) as dx shrank — wrong way / "
                  "noisy. Inspect the per-run traces and Vmax before trusting.")
    else:
        print("  Not enough clean rows to compare — check contamination/CRASH flags "
              "above; scale NU4 ∝ dx⁴ on any blown-up finer grid and rerun.")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_rmax():
    """GATE-BETA Rmax PROBE — does the aim depend on CORE SIZE (vortex preserved, ∇⁴)?
    Vary Rmax {31,46,75 km} at fixed Vmax=64, R_env=500 (core only).  aim rotates WEST
    as Rmax↑ → broader core → more NW → core structure is the lever, and the realistic
    compact Rmax (31-46 vs our 75) drifts MORE poleward (worse — reviewer's concern realized).
    flat → core size is NOT it → the Arm A ∇² rotation was broadening/decay, not size.
    nx=480 so the 31 km core is ~3·dx.  Usage: python run_translation_test.py gate-beta-rmax
    """
    import math
    t0 = time.time()
    TH_HDG = (290.0, 335.0); CAP = 70.0; NX = 480
    RMAX_SWEEP = (31_000.0, 46_000.0, 75_000.0)   # obs-range → our value (ref)
    print("=" * 78)
    print(f"GATE-BETA Rmax PROBE  (β, u=v=0, Vmax=64, R_env=500, taper 200, cap 70, "
          f"∇⁴, nx={NX}/5000km, 48h)")
    print(f"  ref Rmax=75km ~351°; NW band {TH_HDG[0]:.0f}-{TH_HDG[1]:.0f}. aim rotates "
          "WEST as Rmax↑ → core size is the lever (compact = worse); flat → it isn't")
    print("=" * 78)
    rows = []
    for rmax in RMAX_SWEEP:
        dx_km = 5_000_000.0 / NX / 1e3
        try:
            d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=CAP, beta=True,
                                wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                                nx=NX, dom=5_000_000.0, hours=48.0, Rmax=rmax,
                                f_ref=IVAN_F_REF)
        except Exception as e:
            rows.append((rmax, float("nan"), float("nan"), float("nan"),
                         float("nan"), "CRASH"))
            print(f"\n  Rmax={rmax/1e3:.0f}km: ⚠ raised ({type(e).__name__}).")
            continue
        track = d["track"]; tt, xs, ys, vmt = track; ve = d["vmax_end"]
        finite = all(math.isfinite(v) for v in (xs[-1], ys[-1], ve))
        if finite:
            spd, hdg, west = _mature_drift(track)
        else:
            spd = hdg = west = float("nan")
        res = rmax / (dx_km * 1e3)
        contam = finite and ve > CAP * 1.10
        verdict = ("CONTAM" if contam else ("OK" if finite else "UNSTABLE"))
        rows.append((rmax, hdg, spd, west, ve, verdict))
        print(f"\n  Rmax={rmax/1e3:.0f}km ({res:.1f}·dx):")
        if finite:
            print(f"    MATURE(30-48) |{spd:.2f}| hdg {hdg:.0f} ({_compass(hdg)})  "
                  f"WEST {west:+.2f}  Vmax_end={ve:.1f}  "
                  f"{'⚠CONTAM' if contam else 'clean'}"
                  f"{'  ⚠UNDER-RESOLVED core' if res < 3 else ''}")
        else:
            print("    ⚠ NON-FINITE.")
    print("\n" + "=" * 78); print("SUMMARY:")
    print(f"  {'Rmax_km':>7}{'·dx':>5}{'hdg':>6}{'Δvs75':>7}{'|drift|':>8}{'west':>7}"
          f"{'Vmax':>7}  verdict")
    ref = next((r[1] for r in rows if r[0] == 75_000.0 and r[1] == r[1]), float("nan"))
    for rmax, hdg, spd, west, ve, vd in rows:
        hs = f"{hdg:.0f}" if hdg == hdg else "nan"
        dl = (f"{_circ_delta(hdg, ref):+.0f}"
              if (hdg == hdg and ref == ref) else "  —")
        ss = f"{spd:.2f}" if spd == spd else "nan"
        ws = f"{west:+.2f}" if west == west else "nan"
        vs = f"{ve:.0f}" if ve == ve else "nan"
        rr = rmax / (5_000_000.0/NX)
        print(f"  {rmax/1e3:>7.0f}{rr:>5.1f}{hs:>6}{dl:>7}{ss:>8}{ws:>7}{vs:>7}  {vd}")
    print("\nREAD:")
    ok = [r for r in rows if r[5] == "OK" and r[1] == r[1]]
    if len(ok) >= 2 and ref == ref:
        broadest = max(ok, key=lambda r: r[0])
        narrowest = min(ok, key=lambda r: r[0])
        dl_n = _circ_delta(narrowest[1], ref)
        if _circ_delta(broadest[1], ref) - dl_n <= -6:
            print(f"  aim rotates WEST as Rmax↑ (compact {narrowest[0]/1e3:.0f}km "
                  f"{narrowest[1]:.0f}° → broad {broadest[0]/1e3:.0f}km {broadest[1]:.0f}°) → "
                  "CORE SIZE is the lever. ⚠ implication: the realistic compact Rmax (31-46 km) "
                  "drifts MORE poleward than our 75 km → the aim bias is WORSE at observed core "
                  "size, not better. Per-storm Rmax from EBTRK becomes a correctness issue, not "
                  "just skill. (Cross-check: does init-Rmax-31 ≈ ∇²-decayed aim at matched Vmax?)")
        else:
            print(f"  aim ~FLAT across Rmax (Δ within ~5°) → core SIZE is NOT the lever. → the "
                  "Arm A ∇² rotation was the WEAKENING/decay process, not size → strength (Arm 2) "
                  "+ the energy audit are the live threads.")
    else:
        print("  Too few clean rows — check the 31 km row's resolution (needs ≥3·dx).")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_vmax():
    """GATE-BETA Vmax PROBE — does the aim depend on STRENGTH (preserved, ∇⁴)?
    Vary init Vmax {64,50,35,21} at fixed Rmax=75, ∇⁴.  Compare (Vmax_end, aim) to Arm A's
    ∇²-decay curve (42→350,35→347,21→338,7→318): ON-curve → aim is f(Vmax) intrinsic strength;
    OFF-curve (weak init still ~350°) → Arm A rotation was broadening/decay, not strength.
    ⚠ if self-intensification pulls weak inits back toward the cap, Vmax_end won't track init —
    that implicates the leak (→ energy audit).  Usage: python run_translation_test.py gate-beta-vmax
    """
    import math
    t0 = time.time()
    TH_HDG = (290.0, 335.0); CAP = 70.0
    VMAX_SWEEP = (64.0, 50.0, 35.0, 21.0)
    print("=" * 78)
    print("GATE-BETA Vmax PROBE  (β, u=v=0, Rmax=75, R_env=500, taper 200, cap 70, "
          "∇⁴, 320/5000km, 48h)")
    print("  compare (Vmax_end, aim) to Arm A ∇²-decay: (42,350)(35,347)(21,338)(7,318)")
    print("  ON-curve → intrinsic strength; OFF → broadening/decay was the Arm A driver")
    print("=" * 78)
    rows = []
    for vm0 in VMAX_SWEEP:
        try:
            d = run_translation(vm0, u_env=0.0, v_env=0.0, v_cap=CAP, beta=True,
                                wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                                nx=320, dom=5_000_000.0, hours=48.0,
                                f_ref=IVAN_F_REF)
        except Exception as e:
            rows.append((vm0, float("nan"), float("nan"), float("nan"),
                         float("nan"), "CRASH"))
            print(f"\n  init Vmax={vm0:.0f}: ⚠ raised ({type(e).__name__}).")
            continue
        track = d["track"]; tt, xs, ys, vmt = track; ve = d["vmax_end"]
        finite = all(math.isfinite(v) for v in (xs[-1], ys[-1], ve))
        if finite:
            spd, hdg, west = _mature_drift(track)
        else:
            spd = hdg = west = float("nan")
        verdict = ("OK" if finite else "UNSTABLE")
        rows.append((vm0, hdg, spd, west, ve, verdict))
        print(f"\n  init Vmax={vm0:.0f}:")
        if finite:
            drift = "self-intensified" if ve > vm0 + 3 else (
                    "decayed" if ve < vm0 - 3 else "preserved")
            print(f"    Vmax_end={ve:.1f} ({drift})  MATURE |{spd:.2f}| hdg {hdg:.0f} "
                  f"({_compass(hdg)})  WEST {west:+.2f}")
        else:
            print("    ⚠ NON-FINITE.")
    print("\n" + "=" * 78); print("SUMMARY:")
    print(f"  {'init':>5}{'Vmax_end':>9}{'hdg':>6}{'|drift|':>8}{'west':>7}"
          f"{'ArmA@Vend':>10}  verdict")
    # rough Arm-A curve for the comparison column
    armA = [(42, 350), (35, 347), (21, 338), (7, 318)]
    def _armA_at(v):
        if v != v:
            return float("nan")
        pts = sorted(armA)
        if v <= pts[0][0]:
            return pts[0][1]
        if v >= pts[-1][0]:
            return pts[-1][1]
        for (a, ha), (b, hb) in zip(pts, pts[1:]):
            if a <= v <= b:
                return ha + (hb - ha) * (v - a) / (b - a)
        return float("nan")
    for vm0, hdg, spd, west, ve, vd in rows:
        hs = f"{hdg:.0f}" if hdg == hdg else "nan"
        ss = f"{spd:.2f}" if spd == spd else "nan"
        ws = f"{west:+.2f}" if west == west else "nan"
        vs = f"{ve:.1f}" if ve == ve else "nan"
        aa = _armA_at(ve); aas = f"{aa:.0f}" if aa == aa else "—"
        print(f"  {vm0:>5.0f}{vs:>9}{hs:>6}{ss:>8}{ws:>7}{aas:>10}  {vd}")
    print("\nREAD:")
    ok = [r for r in rows if r[5] == "OK" and r[1] == r[1]]
    if ok:
        # are the preserved-vortex aims near Arm A's curve at the same Vmax_end?
        diffs = [abs(_circ_delta(r[1], _armA_at(r[4]))) for r in ok if _armA_at(r[4]) == _armA_at(r[4])]
        on = diffs and max(diffs) <= 6
        intens = [r for r in ok if r[4] > r[0] + 3]
        if intens:
            print(f"  ⚠ self-intensification: inits {[int(r[0]) for r in intens]} ended HIGHER than "
                  "init → the leak is real and resists making a weak vortex. The energy audit is now "
                  "the priority; structure-by-strength can't be cleanly tested while the leak refills it.")
        if on:
            print("  Preserved-vortex aims track Arm A's curve at matched Vmax_end → the aim is an "
                  "intrinsic function of vortex STRENGTH (path-independent). The bias is real "
                  "intense-vortex β-drift physics in this model, not a dissipation artifact.")
        else:
            print("  Preserved-vortex aims do NOT match Arm A at matched Vmax_end (weak inits stay "
                  "more poleward) → the Arm A rotation was the BROADENING/decay, not strength → "
                  "core size (Arm 1) and/or the leak are the seat.")
    else:
        print("  No clean rows.")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_fdecomp():
    """GATE-BETA f-PLANE DECOMPOSITION (V8.7) — isolate spurious drift from true β-drift.
    aim=f(Vmax_end), decomp WEST≈const~0.4 / NORTH∝Vmax.  On an f-plane the true drift
    is identically 0, so any f-plane drift is ARTIFACT; β−f = true β-drift.  Run each
    Vmax β ON and OFF.  Predict f-plane NORTH ∝ Vmax, WEST≈0 → spurious strength-coupled
    poleward self-translation.  f-plane ≈0 → β-plane north IS β-drift, mis-oriented.
    No harness edit (beta is a flag).  Usage: python run_translation_test.py gate-beta-fdecomp
    """
    import math
    t0 = time.time()
    CAP = 70.0
    VMAX_SWEEP = (64.0, 35.0, 21.0)

    def _nw(track):
        spd, hdg, west = _mature_drift(track)
        north = spd * math.cos(math.radians(hdg))
        return north, west, spd, hdg

    print("=" * 78)
    print("GATE-BETA f-PLANE DECOMPOSITION  (u=v=0, R_env=500, taper 200, cap 70, ∇⁴, "
          "320/5000km, 48h)")
    print("  f-plane true drift = 0 → any f-plane drift is ARTIFACT;  β − f = true β-drift")
    print("  predict: f-plane NORTH ∝ Vmax, WEST≈0;  true β-drift = the residual")
    print("=" * 78)
    rows = []
    for vm0 in VMAX_SWEEP:
        out = {}
        for beta_on in (False, True):
            try:
                d = run_translation(vm0, u_env=0.0, v_env=0.0, v_cap=CAP, beta=beta_on,
                                    wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                                    nx=320, dom=5_000_000.0, hours=48.0,
                                    f_ref=IVAN_F_REF)
                track = d["track"]; ve = d["vmax_end"]
                finite = all(math.isfinite(v)
                             for v in (track[1][-1], track[2][-1], ve))
                out[beta_on] = (_nw(track) if finite else None, ve)
            except Exception as e:
                out[beta_on] = (None, float("nan"))
                print(f"  Vmax {vm0:.0f} beta={beta_on}: raised {type(e).__name__}")
        rows.append((vm0, out))
        fp = out[False][0]; bp = out[True][0]
        print(f"\n  init Vmax={vm0:.0f}:")
        if fp:
            print(f"    f-plane (β off):  NORTH {fp[0]:+.2f}  WEST {fp[1]:+.2f}  "
                  f"|{fp[2]:.2f}|@{fp[3]:.0f}  Vend={out[False][1]:.1f}   ← SPURIOUS")
        if bp:
            print(f"    β-plane (β on):   NORTH {bp[0]:+.2f}  WEST {bp[1]:+.2f}  "
                  f"|{bp[2]:.2f}|@{bp[3]:.0f}  Vend={out[True][1]:.1f}")
        if fp and bp:
            print(f"    β − f (TRUE β-drift): NORTH {bp[0]-fp[0]:+.2f}  "
                  f"WEST {bp[1]-fp[1]:+.2f}")
    print("\n" + "=" * 78); print("SUMMARY (true β-drift = β − f):")
    print(f"  {'Vmax':>5}{'fp_N':>7}{'fp_W':>7}{'bp_N':>7}{'bp_W':>7}{'true_N':>8}{'true_W':>8}")
    for vm0, out in rows:
        fp = out[False][0]; bp = out[True][0]
        if fp and bp:
            print(f"  {vm0:>5.0f}{fp[0]:>7.2f}{fp[1]:>7.2f}{bp[0]:>7.2f}{bp[1]:>7.2f}"
                  f"{bp[0]-fp[0]:>8.2f}{bp[1]-fp[1]:>8.2f}")
        else:
            print(f"  {vm0:>5.0f}  (incomplete)")
    print("\nREAD:")
    good = [(vm0, out) for vm0, out in rows if out[False][0] and out[True][0]]
    if good:
        vms = [vm0 for vm0, _ in good]
        fns = [out[False][0][0] for _, out in good]
        fws = [out[False][0][1] for _, out in good]
        print("  f-plane NORTH vs Vmax: " +
              ", ".join(f"{v:.0f}:{n:+.2f}" for v, n in zip(vms, fns)))
        print("  f-plane WEST  vs Vmax: " +
              ", ".join(f"{v:.0f}:{w:+.2f}" for v, w in zip(vms, fws)))
        if max(abs(n) for n in fns) > 0.3:
            scales = (abs(fns[0]) > abs(fns[-1]) + 0.2) if len(fns) >= 2 else False
            print("  → f-plane vortex DRIFTS with β OFF = pure ARTIFACT (true f-plane drift is 0).")
            if scales:
                print("    NORTH grows with Vmax → spurious strength-coupled POLEWARD self-"
                      "translation confirmed. The poleward aim is NOT β-drift. TRUE β-drift = the "
                      "β−f residual (check it: likely weak, ~westward). Source hunt next: advection "
                      "asymmetry (does it flip sign if the vortex spins the other way?) vs the "
                      "self-intensification leak (energy/momentum audit).")
            else:
                print("    (NORTH not clearly Vmax-scaling — inspect; could be a constant grid drift.)")
        else:
            print("  → f-plane drift ≈ 0 → the β-plane north IS β-drift, mis-oriented poleward → a "
                  "gyre-orientation problem (e.g. the Coriolis term's gyre response), not a spurious "
                  "translation. Different hunt: instrument the β-gyre structure directly.")
    else:
        print("  Incomplete — check for crashes above.")
    print(f"\nWall time: {time.time()-t0:.0f}s")


def gate_beta_timeevol():
    """GATE-BETA TIME-EVOLUTION (V8.7) — does the β-drift heading rotate NW→N over the run
    (gyre over-rotation), and is the rotation steeper for stronger vortices?  f-plane decomp
    showed the drift is genuine β-drift mis-oriented poleward.  Over-rotation predicts early
    NW marching to N, steeper for Vmax 64 than 21.  Fixed-direction-from-the-start → static
    asymmetry instead.  Reports drift per window.  Usage: python run_translation_test.py gate-beta-timeevol
    """
    import math
    t0 = time.time()
    CAP = 70.0
    VMAX_SWEEP = (64.0, 21.0)
    WINDOWS = [(6, 12), (12, 24), (24, 36), (36, 48)]
    print("=" * 78)
    print("GATE-BETA TIME-EVOLUTION  (β, u=v=0, R_env=500, taper 200, cap 70, ∇⁴, "
          "320/5000km, 48h)")
    print("  over-rotation → heading climbs NW→N over the run, steeper for Vmax 64")
    print("  fixed direction from the start → static asymmetry (different seat)")
    print("=" * 78)
    summary = {}
    for vm0 in VMAX_SWEEP:
        try:
            d = run_translation(vm0, u_env=0.0, v_env=0.0, v_cap=CAP, beta=True,
                                wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                                nx=320, dom=5_000_000.0, hours=48.0,
                                f_ref=IVAN_F_REF)
        except Exception as e:
            print(f"\n  Vmax {vm0:.0f}: raised {type(e).__name__}")
            continue
        track = d["track"]
        print(f"\n  Vmax={vm0:.0f} (Vend={d['vmax_end']:.1f}):")
        hdgs = []
        for (ta, tb) in WINDOWS:
            spd, hdg, west = _mature_drift(track, ta, tb)
            north = spd * math.cos(math.radians(hdg))
            hdgs.append(hdg)
            print(f"    t{ta:02d}-{tb:02d}h: |{spd:.2f}| hdg {hdg:.0f} "
                  f"({_compass(hdg)})  N {north:+.2f}  W {west:+.2f}")
        # net rotation over the run (signed circular, + = poleward/clockwise toward N)
        rot = _circ_delta(hdgs[-1], hdgs[0])
        summary[vm0] = (hdgs, rot)
        print(f"    → net heading change {hdgs[0]:.0f}→{hdgs[-1]:.0f} = {rot:+.0f}° "
              f"({'POLEWARD/over-rotating' if rot > 4 else 'westward' if rot < -4 else 'flat'})")
    print("\n" + "=" * 78); print("READ:")
    if len(summary) == 2:
        r64 = summary[64.0][1]; r21 = summary[21.0][1]
        h64 = summary[64.0][0]; h21 = summary[21.0][0]
        if r64 > 4 and r64 > r21 + 3:
            print(f"  Heading ROTATES POLEWARD over the run (Vmax64 {h64[0]:.0f}→{h64[-1]:.0f} "
                  f"= {r64:+.0f}°) and STEEPER than Vmax21 ({r21:+.0f}°) → GYRE OVER-ROTATION "
                  "CONFIRMED. Seat = gyre EQUILIBRATION (vortex swirl out-runs it). Next: instrument "
                  "the wavenumber-1 gyre tilt directly, then ask what should arrest the rotation "
                  "(Rossby-wave dispersion radiating gyre vorticity) and whether it's impeded here.")
        elif abs(r64) <= 4:
            print(f"  Heading ~FIXED from the start (Vmax64 change {r64:+.0f}°) → NOT over-rotation "
                  "→ a STATIC asymmetry (vortex profile / advection) aims the ventilation poleward "
                  "from t=0. Seat = the asymmetric forcing, not equilibration → instrument the "
                  "initial gyre formation, not its time-evolution.")
        else:
            print(f"  Mixed: Vmax64 rotates {r64:+.0f}°, Vmax21 {r21:+.0f}° — inspect the per-window "
                  "trace; the Vmax-dependence of the rotation is the discriminator.")
    else:
        print("  Incomplete — check for crashes.")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_longrun():
    """GATE-BETA LONG-RUN TRIAGE (V8.7) — does the over-rotating β-drift heading PLATEAU
    (gyres equilibrate) or RUN AWAY past north?  Extend to 96h, Vmax 64; track heading AND
    Vmax per 12h window.  Climbs past N despite decay → RUNAWAY.  Plateaus while Vmax still
    drops → real equilibration (too-poleward).  Turns W → decay confound (→ instrument gyre).
    No harness edit.  Usage:  python run_translation_test.py gate-beta-longrun
    """
    import math
    t0 = time.time()
    CAP = 70.0; HOURS = 96.0
    WINDOWS = [(6, 12), (12, 24), (24, 36), (36, 48),
               (48, 60), (60, 72), (72, 84), (84, 96)]
    print("=" * 78)
    print(f"GATE-BETA LONG-RUN TRIAGE  (β, u=v=0, Vmax=64, R_env=500, taper 200, cap 70, "
          f"∇⁴, 320/5000km, {HOURS:.0f}h)")
    print("  climbs past N → RUNAWAY;  plateaus (Vmax still dropping) → equilibration;  "
          "turns W → decay confound")
    print("=" * 78)
    try:
        d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=CAP, beta=True,
                            wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                            nx=320, dom=5_000_000.0, hours=HOURS,
                            f_ref=IVAN_F_REF)
    except Exception as e:
        print(f"  ⚠ run raised {type(e).__name__} — likely late blow-up; "
              "try HOURS=72.")
        return
    track = d["track"]; tt = track[0]; vmt = track[3]

    def _vmax_at(th):
        i = min(range(len(tt)), key=lambda k: abs(tt[k] - th * 3600.0))
        return vmt[i]

    print(f"\n  Vmax_end({HOURS:.0f}h)={d['vmax_end']:.1f}\n")
    rows = []
    prev_h = prev_c = None
    for (ta, tb) in WINDOWS:
        spd, hdg, west = _mature_drift(track, ta, tb)
        north = spd * math.cos(math.radians(hdg))
        ctr = (ta + tb) / 2.0
        vm = _vmax_at(ctr)
        rate_s = ""
        if prev_h is not None and hdg == hdg:
            dh = _circ_delta(hdg, prev_h)
            rate_s = f"  rate {dh/(ctr-prev_c):+.2f}°/h"
        rows.append((ctr, hdg, spd, north, west, vm))
        hs = f"{hdg:.0f}" if hdg == hdg else "nan"
        print(f"    t{ta:02d}-{tb:02d}h: |{spd:.2f}| hdg {hs} ({_compass(hdg) if hdg==hdg else '--'})"
              f"  N {north:+.2f} W {west:+.2f}  Vmax {vm:.1f}{rate_s}")
        if hdg == hdg:
            prev_h = hdg; prev_c = ctr
    print("\n" + "=" * 78); print("READ:")
    fin = [r for r in rows if r[1] == r[1]]
    if len(fin) < 3:
        print("  Too few finite windows — late blow-up; rerun at HOURS=72.")
        print(f"\nWall time: {time.time()-t0:.0f}s"); return
    h0 = fin[0][1]; hL = fin[-1][1]
    net = _circ_delta(hL, h0)
    late_rate = _circ_delta(fin[-1][1], fin[-2][1]) / (fin[-1][0] - fin[-2][0])
    print(f"  net heading {h0:.0f}→{hL:.0f} = {net:+.0f}° ;  late rate {late_rate:+.2f}°/h ;  "
          f"Vmax {fin[0][5]:.0f}→{fin[-1][5]:.0f}")
    if hL >= 358 or hL < h0 - 30 or net > 32:
        print("  → heading reached/passed NORTH despite decay → RUNAWAY over-rotation: the gyres "
              "never equilibrate. The β-drift orientation is unbounded in this config → the bug is a "
              "missing/weak gyre-equilibration mechanism (Rossby-wave radiation of gyre vorticity). "
              "Gyre instrumentation next to confirm & locate.")
    elif abs(late_rate) < 0.10:
        print(f"  → heading PLATEAUED (~{hL:.0f}°) while Vmax still dropping → likely REAL "
              "equilibration at a too-poleward angle. The model HAS a steady β-drift, just mis-"
              "oriented → instrument the gyre to see why the steady tilt is wrong.")
    elif late_rate < -0.05:
        print("  → heading turned WESTWARD late → decay confound (weaker vortex → westward offset) "
              "overtaking → INCONCLUSIVE on equilibration → gyre instrumentation (decay-robust) required.")
    else:
        print(f"  → still climbing at {late_rate:+.2f}°/h, not yet at N → slow over-rotation, not "
              "clearly runaway or plateaued → extend to 120h OR go to gyre instrumentation.")
    print(f"\nWall time: {time.time()-t0:.0f}s")

def gate_beta_gyre(vmax=64.0):
    """GATE-BETA GYRE INSTRUMENTATION (V8.7) — SEE the β-gyre precess.
    Snapshots full state at t=12/24/36/48/60h, computes relative vorticity, extracts the
    wavenumber-1 (m=1) gyre asymmetry about the tracked center, reports PHASE (orientation)
    and AMPLITUDE vs time.  Expect amplitude saturates, phase precesses at ~the drift-heading
    rate → confirms a freely-precessing, non-phase-locked gyre.  Saves ζ fields to .npz for
    the figure.  Requires the 4 harness edits (snapshot hook).
    Usage:  python run_translation_test.py gate-beta-gyre
    """
    import numpy as np, math
    t0 = time.time()
    DOM = 5_000_000.0; NX = 320
    SNAP_H = [12.0, 24.0, 36.0, 48.0, 60.0]
    print("=" * 78)
    print(f"GATE-BETA GYRE INSTRUMENTATION  (β, u=v=0, Vmax={vmax:.0f}, R_env=500, taper 200, "
          f"cap 70, ∇⁴, 320/5000km, 60h)")
    print("  m=1 vorticity asymmetry: PHASE = gyre orientation, |A| = gyre amplitude")
    print("  expect amplitude SATURATES, phase PRECESSES at ~the drift-heading rate")
    print("=" * 78)
    d = run_translation(vmax, u_env=0.0, v_env=0.0, v_cap=70.0, beta=True,
                        wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                        nx=NX, dom=DOM, hours=60.0, snapshot_hours=SNAP_H,
                        f_ref=IVAN_F_REF)
    snaps = d.get("snapshots", [])
    if not snaps:
        print("  ⚠ no snapshots returned — the 4 harness edits aren't applied.")
        return
    dx = DOM / NX
    xs = (np.arange(NX) + 0.5) * dx
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    print(f"\n  {'t(h)':>5}{'gyre_phase':>11}{'Δ°/h':>8}{'gyre_amp':>11}"
          f"{'cx_km':>8}{'cy_km':>8}")
    rows = []
    prev_ph = prev_t = None
    for s in snaps:
        u2 = s["u"].mean(axis=2); v2 = s["v"].mean(axis=2)
        dvdx = (np.roll(v2, -1, 0) - np.roll(v2, 1, 0)) / (2 * dx)
        dudy = (np.roll(u2, -1, 1) - np.roll(u2, 1, 1)) / (2 * dx)
        zeta = dvdx - dudy
        RX = X - s["cx"]; RY = Y - s["cy"]
        R = np.hypot(RX, RY); TH = np.arctan2(RY, RX)
        band = (R > 75e3) & (R < 450e3)
        A = np.sum(zeta[band] * np.exp(-1j * TH[band]))
        amp = float(np.abs(A))
        phase = (-math.degrees(math.atan2(A.imag, A.real))) % 360.0
        th = s["t"] / 3600.0
        dps = ""
        if prev_ph is not None:
            dp = ((phase - prev_ph + 180.0) % 360.0) - 180.0
            dps = f"{dp/(th-prev_t):+.2f}"
        rows.append((th, phase, amp))
        print(f"  {th:>5.0f}{phase:>11.0f}{dps:>8}{amp:>11.2e}"
              f"{s['cx']/1e3:>8.0f}{s['cy']/1e3:>8.0f}")
        prev_ph = phase; prev_t = th
        np.savez(f"gyre_snap_v{int(vmax):02d}_t{int(th):02d}h.npz",
                 zeta=zeta, u=u2, v=v2, cx=s["cx"], cy=s["cy"], dx=dx)
    print("\n" + "=" * 78); print("READ:")
    if len(rows) >= 3:
        amps = [r[2] for r in rows]
        ph_rate = ((((rows[-1][1]-rows[0][1]+180) % 360)-180)
                   / (rows[-1][0]-rows[0][0]))
        amp_sat = (amps[-1] < amps[-2]*1.15) and (amps[-2] < amps[-3]*1.3)
        tt, xt, yt, vmt = d["track"]
        def _hdg(ta, tb):
            ia = min(range(len(tt)), key=lambda k: abs(tt[k]-ta*3600))
            ib = min(range(len(tt)), key=lambda k: abs(tt[k]-tb*3600))
            return math.degrees(math.atan2(xt[ib]-xt[ia], yt[ib]-yt[ia])) % 360
        drift_rate = ((((_hdg(48,60)-_hdg(12,24)+180) % 360)-180) / 36.0)
        print(f"  gyre m=1 PHASE rate {ph_rate:+.2f}°/h ;  drift HEADING rate "
              f"{drift_rate:+.2f}°/h  (magnitudes should match; sign differs by convention)")
        print(f"  gyre AMPLITUDE {amps[0]:.2e} → {amps[-1]:.2e} "
              f"({'SATURATING' if amp_sat else 'still growing'})")
        if abs(ph_rate) > 0.15 and amp_sat and abs(abs(ph_rate)-abs(drift_rate)) < 0.25:
            print("  → CONFIRMED: m=1 gyre PRECESSES (phase rotates) at ~the drift-heading rate "
                  "while its AMPLITUDE saturates → a freely-precessing, NON-phase-locked gyre IS "
                  "the drift bias. The .npz fields are the paper figure (the dipole rotating). "
                  "Next: WHY no lock — inspect the saved ζ far field for a radiating Rossby wave "
                  "train (gyre vorticity propagating away instead of staying co-located).")
        elif abs(ph_rate) > 0.15 and amp_sat:
            print(f"  → gyre precesses + amplitude saturates, but the phase rate ({ph_rate:+.2f}) "
                  f"doesn't match the drift rate ({drift_rate:+.2f}) — the m=1 band (75-450 km) may "
                  "be catching the wrong structure; inspect the .npz and adjust the annulus.")
        else:
            print("  → phase ~steady or amplitude still growing — inspect the .npz fields and the "
                  "band cuts; may need a wider/narrower annulus or a single level vs vertical mean.")
    tag = f"v{int(vmax):02d}"
    print(f"\n  saved {len(rows)} zeta fields (gyre_snap_{tag}_t*.npz) - "
          "plot zeta over (x-cx, y-cy) to see the dipole; "
          f"reanalyze with: python oracle_v8/reanalyze_gyre.py --tag {tag}")
    print(f"Wall time: {time.time()-t0:.0f}s")

if __name__ == "__main__":
    import sys
    arg = sys.argv[1].lower() if len(sys.argv) > 1 else ""
    if arg in ("7", "disc", "discriminator"):
        discriminator_only()
    elif arg in ("beta", "betadrift", "8"):
        betadrift_only()
    elif arg in ("effx", "xtrans", "9"):
        effx_only()
    elif arg in ("gate", "10"):
        gate_only()
    elif arg in ("gate-drift", "gatedrift", "gdrift", "11"):
        gate_drift_only()
    elif arg in ("gate-beta", "gatebeta", "gbeta", "12"):
        gate_beta_only()
    elif arg in ("ladder", "13"):
        ladder_only()
    elif arg in ("trace", "14"):
        trace_only()
    elif arg in ("intensify-ladder", "iladder", "16"):
        intensify_ladder()
    elif arg in ("gate-beta-renv", "gbeta-renv", "renv", "17"):
        gate_beta_renv()
    elif arg in ("gate-beta-domain", "gbeta-dom", "domain", "18"):
        gate_beta_domain()
    elif arg in ("gate-beta-taper", "gbeta-taper", "taper", "19"):
        gate_beta_taper()
    elif arg in ("gate-beta-nu4", "gbeta-nu4", "nu4", "20"):
        gate_beta_nu4()
    elif arg in ("gate-beta-res", "gbeta-res", "res", "21"):
        gate_beta_res()
    elif arg in ("gate-beta-diffform", "gbeta-diffform", "diffform", "22"):
        gate_beta_diffform()
    elif arg in ("gate-beta-divdamp", "gbeta-divdamp", "divdamp", "23"):
        gate_beta_divdamp()
    elif arg in ("gate-beta-rmax", "gbeta-rmax", "rmax", "24"):
        gate_beta_rmax()
    elif arg in ("gate-beta-vmax", "gbeta-vmax", "vmax", "25"):
        gate_beta_vmax()
    elif arg in ("gate-beta-fdecomp", "gbeta-fdecomp", "fdecomp", "26"):
        gate_beta_fdecomp()
    elif arg in ("gate-beta-timeevol", "gbeta-timeevol", "timeevol", "27"):
        gate_beta_timeevol()
    elif arg in ("gate-beta-longrun", "gbeta-longrun", "longrun", "28"):
        gate_beta_longrun()
    elif arg in ("gate-beta-gyre", "gbeta-gyre", "gyre", "29"):
        _vm = float(sys.argv[2]) if len(sys.argv) > 2 else 64.0
        gate_beta_gyre(_vm)
    else:
        main()
