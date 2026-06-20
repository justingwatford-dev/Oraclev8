"""
Oracle V8 — Unified storm runner
================================
One runner for any HURDAT2 Atlantic TC.  Replaces run_hugo.py / run_katrina.py /
run_ivan.py.  The physics config is the single shared production_config; the only
per-storm inputs are the HURDAT2 id + init time (everything else read from the best
track), and the domain size / run length are derived from the init→landfall geometry.

    python -m oracle_v8.run_storm Ivan
    python -m oracle_v8.run_storm Hugo
    python -m oracle_v8.run_storm AL142018 --init "2018-10-09 12Z"   # any storm, by id

Before a new storm's first run:  python -m oracle_v8.era5_steering --download --storm <name>
(REQUIRE_ERA5=True — a missing ERA5 file aborts rather than silently downgrading.)

Mirrors the validated run_ivan.py loop verbatim; the only changes are that the storm
dict, domain, and run length come from hurdat2 + production_config instead of being
hand-set per file.
"""
from __future__ import annotations

import sys as _sys
if hasattr(_sys.stdout, "reconfigure"):
    _sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import sys
import time
from dataclasses import replace as dc_replace

import numpy as np

from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import IntensityCapComponent, RK3Integrator
from oracle_v8.backend import xp, wrap_base
from oracle_v8.storm_tracker import find_storm_centre as _track_centre
from oracle_v8.era5_steering import ERA5Steering
from oracle_v8.landfall_verify import landfall_report

from oracle_v8 import hurdat2
from oracle_v8.production_config import (
    DT, DIAG_EVERY, N_PREBAL, NZ, LZ,
    RMAX_RUN_M, R_ENV_M, TAPER_START_FRAC, WIND_TAPER,
    VMAX_CAP_MS, TAU_CAP, REQUIRE_ERA5, TIME_VARYING_STEER, TAU_STEER,
    choose_domain, n_steps_for, build_base_state,
    build_prebal_config, build_production_config,
)

try:
    from oracle_v8 import diagnostics as _vdiag
    _HAVE_DIAG = True
except Exception:
    _HAVE_DIAG = False


def centre_to_latlon(x_c, y_c, lat0, lon0, Lx, Ly):
    R = 6_371_000.0
    dlat = np.degrees((y_c - Ly / 2.0) / R)
    dlon = np.degrees((x_c - Lx / 2.0) / (R * np.cos(np.radians(lat0))))
    return lat0 + dlat, lon0 + dlon


def main(name: str, init_override: str | None = None) -> int:
    # ---- Storm init + verification, straight from HURDAT2 ------------------
    s = hurdat2.load_storm(name, init_date=init_override)

    # ---- Geometry-derived domain & run length ------------------------------
    nx, Lx = choose_domain(s["lat0_deg"], s["threshold_lat"])
    ny, Ly = nx, Lx
    nz, Lz = NZ, LZ
    dx = Lx / nx
    N_STEPS = n_steps_for(s["landfall_time_h"])

    zc, rho0_arr, theta0_arr = build_base_state(nz, Lz)

    class Base:
        z = zc
        rho0 = rho0_arr
        theta0 = theta0_arr

    print("=" * 68)
    print(f"ORACLE V8 — HURRICANE {s['name'].upper()} ({s['year']}) TRACK  [unified run_storm]")
    print("=" * 68)
    print(f"  HURDAT2 id:   {s['hurdat2_id']}   init {s['init_date']}  (source: {s['init_source']})")
    print(f"  Position:     {s['lat0_deg']:.1f}°N, {abs(s['lon0_deg']):.1f}°W")
    print(f"  Vmax:         {s['Vmax_ms']:.1f} m/s ({s['Vmax_kt']:.0f} kt)   P_min {s['P_min_mb']} mb")
    print(f"  Rmax run:     {RMAX_RUN_M/1000:.0f} km (5×dx)   B {s['B']} (frozen)")
    print(f"  Domain:       {Lx/1e3:.0f} km, nx={nx} (dx={dx/1e3:.3f} km) — geometry-derived")
    print(f"  Wind taper:   {'ON' if WIND_TAPER else 'OFF'}  R_env={R_ENV_M/1000:.0f} km  "
          f"taper-start frac {TAPER_START_FRAC:.2f}")
    print(f"  Run:          {N_STEPS} steps = {N_STEPS*DT/3600:.0f} h  (dt {DT:.0f}s, "
          f"CFL {s['Vmax_ms']*DT/dx:.3f})")
    print(f"  Obs landfall: {s['landfall_lat']}°N, {abs(s['landfall_lon'])}°W  "
          f"t+{s['landfall_time_h']:.2f}h   threshold {s['threshold_lat']}°N")

    # ---- ERA5 DLM steering — load FIRST so t=0 seeds the vortex ------------
    print("\nLoading ERA5 DLM steering winds...")
    try:
        era5 = ERA5Steering.load(storm=s["name"].lower())
        era5.print_summary()
        USE_ERA5 = True
        if TIME_VARYING_STEER:
            u_env_t0, v_env_t0 = era5.get_dlm(0.0, s["lat0_deg"], abs(s["lon0_deg"]))
            print(f"  ERA5 steering: ACTIVE  (time-varying — init at t=0 DLM "
                  f"u={u_env_t0:+.2f} v={v_env_t0:+.2f}, relax τ={TAU_STEER/3600:.1f} h)")
        else:
            u_env_t0, v_env_t0 = era5.get_track_mean_dlm()
            print("  ERA5 steering: ACTIVE  (frozen track-mean — A/B only)")
    except FileNotFoundError as exc:
        if REQUIRE_ERA5:
            print(f"\n  FATAL: {exc}")
            print("  REQUIRE_ERA5=True — refusing to run with fallback constant steering")
            print(f"  (not the production config).  Download ERA5 for '{s['name'].lower()}'")
            print("  or set REQUIRE_ERA5=False in production_config for an explicit A/B run.")
            return 1
        print(f"  WARNING: {exc}\n  Falling back to HURDAT2 constant steering.")
        u_env_t0, v_env_t0 = s["u_env_ms"], s["v_env_ms"]
        era5, USE_ERA5 = None, False

    # ---- Init vortex -------------------------------------------------------
    print("\nInitializing vortex...")
    t0 = time.time()
    init = HollandVortexInit(
        Vmax=s["Vmax_ms"], Rmax=RMAX_RUN_M, B=s["B"], f=s["f"],
        R_env=R_ENV_M, wind_taper=WIND_TAPER, taper_start_frac=TAPER_START_FRAC,
        u_env=u_env_t0, v_env=v_env_t0,
    )
    state = init.build_state(nx, ny, nz, Lx, Ly, Base())
    state = dc_replace(state, theta_prime=xp.zeros_like(state.theta_prime))
    print(f"  {time.time()-t0:.1f}s  {init.summary()}")

    # ---- Pre-balance: remove grid-scale divergence -------------------------
    print(f"\nPre-balancing ({N_PREBAL} projection-only iterations)...")
    prebal_integrator = RK3Integrator(
        config=build_prebal_config(nx, ny, nz, Lx, Ly, Lz), base=Base())
    for i in range(N_PREBAL):
        state, pdiag = prebal_integrator.step(state, dt=1.0, step_number=i)
        phi_rms = float(xp.sqrt(xp.mean(state.projection_potential ** 2)))
        print(f"  iter {i+1}: phi_rms={phi_rms:.2f} Pa  max|u|={pdiag.max_u:.3f} m/s")
    print("  Pre-balance complete — initial divergence removed.")

    # ---- Production config (the one shared, storm-agnostic physics) --------
    config = build_production_config(nx, ny, nz, Lx, Ly, Lz, s["f"], u_env_t0, v_env_t0)
    integrator = RK3Integrator(config=config, base=Base())
    base_w = wrap_base(Base())

    cap_comp = (IntensityCapComponent(v_cap=VMAX_CAP_MS, tau=TAU_CAP,
                                      u_env=u_env_t0, v_env=v_env_t0)
                if VMAX_CAP_MS is not None else None)
    if cap_comp is not None:
        print(f"  Intensity cap: ACTIVE  |V'| → {VMAX_CAP_MS:.0f} m/s "
              f"(τ={TAU_CAP:.0f}s, perturbation-relative)")

    # ---- Run ---------------------------------------------------------------
    print(f"\n  {'Time(h)':>6}  {'Lat(°N)':>7}  {'Lon(°W)':>7}  "
          f"{'φ_min(Pa)':>12}  {'max|u|(m/s)':>12}")
    print(f"  {'------':>6}  {'-------':>7}  {'-------':>7}  "
          f"{'----------':>12}  {'-----------':>12}")

    track_t, track_lat, track_lon, track_phi = [], [], [], []
    t_run = time.time()
    any_nan = False
    last_pos = None
    last_raw_pos = None
    u_env, v_env = u_env_t0, v_env_t0
    alpha_steer = min(1.0, (DIAG_EVERY * DT) / TAU_STEER)

    for n in range(N_STEPS):
        state, diag = integrator.step(state, dt=DT, step_number=n)

        if cap_comp is not None:
            _ct = cap_comp.compute_tendency(state, None, None, base_w, DT)
            state = dc_replace(state, u=state.u + DT * _ct.du_dt,
                               v=state.v + DT * _ct.dv_dt)

        if n % DIAG_EVERY == 0 and np.isnan(float(diag.max_u)):
            print(f"\n  ✗ NaN at step {n} (t={n*DT/3600:.2f}h)")
            any_nan = True
            break

        if n % DIAG_EVERY == 0:
            t_h = n * DT / 3600.0
            (x_c, y_c), trk_method, trk_diag = _track_centre(
                state, base_w, nx, ny, dx,
                last_pos=last_pos,
                v_env_ms=v_env,
                u_env_ms=u_env,
                dt_output_s=DIAG_EVERY * DT,
                return_diag=True,
                last_raw_pos=last_raw_pos,
            )
            last_pos = (x_c, y_c)
            last_raw_pos = trk_diag['theta_raw_m']
            lat_c, lon_c = centre_to_latlon(x_c, y_c, s["lat0_deg"], s["lon0_deg"], Lx, Ly)
            track_t.append(t_h); track_lat.append(lat_c)
            track_lon.append(lon_c); track_phi.append(diag.surface_phi_min)

            if USE_ERA5:
                u_tgt, v_tgt = era5.get_dlm(t_h, lat_c, abs(lon_c))
                if TIME_VARYING_STEER and n > 0:
                    du = (u_tgt - u_env) * alpha_steer
                    dv = (v_tgt - v_env) * alpha_steer
                    state = dc_replace(state, u=state.u + du, v=state.v + dv)
                    u_env += du
                    v_env += dv
                    config.coriolis.set_env(u_env, v_env)
                    if config.surface_drag is not None:
                        config.surface_drag.set_env(u_env, v_env)
                    if cap_comp is not None:
                        cap_comp.set_env(u_env, v_env)
                    steer_tag = f" [{u_env:+.1f},{v_env:+.1f}→{u_tgt:+.1f},{v_tgt:+.1f}]"
                else:
                    steer_tag = f" [{u_tgt:+.1f},{v_tgt:+.1f}]"
            else:
                steer_tag = ""

            tag = " [spin-up]" if n < 120 else ""
            _TRK_TAG = {'theta_prime': 'the', 'theta_reacq': 'rcq',
                        'vorticity': 'vor', 'extrapolated': 'ext'}
            trk_tag = f" [{_TRK_TAG.get(trk_method, trk_method[:3])}]"

            xr, yr = trk_diag['theta_raw_m']
            lat_raw, lon_raw = centre_to_latlon(xr, yr, s["lat0_deg"], s["lon0_deg"], Lx, Ly)
            raw_tag = f"  θraw={lat_raw:5.2f}N/{abs(lon_raw):5.2f}W"

            print(f"  {t_h:>6.1f}  {lat_c:>7.2f}  {abs(lon_c):>7.2f}  "
                  f"{diag.surface_phi_min:>12.1f}  {diag.max_u:>12.4e}"
                  f"{steer_tag}{trk_tag}{tag}{raw_tag}")

            if _HAVE_DIAG and n > 0:
                try:
                    _dd = _vdiag.compute_diagnostics(
                        state, Base(), x_c, y_c, dx, dx,
                        u_env=u_env, v_env=v_env,
                        theta_xy_m=trk_diag.get("theta_raw_m", (x_c, y_c)))
                    print(f"          {_vdiag.format_diag_line(_dd)}")
                except Exception:
                    pass

    wall = time.time() - t_run
    print(f"\nWall time: {wall:.1f}s")

    if any_nan:
        print("Run FAILED"); return 1

    landfall_report(track_t, track_lat, track_lon, s)
    return 0


if __name__ == "__main__":
    args = [a for a in sys.argv[1:]]
    init_override = None
    if "--init" in args:
        k = args.index("--init")
        init_override = args[k + 1]
        del args[k:k + 2]
    storm = args[0] if args else "Ivan"
    sys.exit(main(storm, init_override))
