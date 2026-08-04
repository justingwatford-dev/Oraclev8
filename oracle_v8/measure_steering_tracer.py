"""
Oracle V8 — steering-only tracer (the steering share of the landfall skill)
===========================================================================
Registered analysis 2 of PAPER2_TIGHTENING_predictions.md (P-S1..P-S4).

A point with no vortex is advected by the annulus DLM sampled at the
tracer's own position from the cached per-storm ERA5 file (production
3-7 degree annulus), RK2 at dt = 0.25 h, from the HURDAT2 init fix to the
observed landfall time; scored with the manuscripts' landfall-fix
along/cross conventions. The model-minus-tracer difference is an
independent read of the self-propagation footprint.

Registered caveat: the tracer feels the DLM along its own (diverged) path
and omits the 3-h relaxation lag — a first-order steering share, not an
exact ablation.

No model runs. Inputs: hurdat2.txt + cached {storm}_era5_steering.nc.

Run:  python -m oracle_v8.measure_steering_tracer
"""
from __future__ import annotations

import math
import os
from datetime import datetime, timedelta

from .hurdat2 import parse_hurdat2
from .era5_steering import ERA5Steering
from .measure_persistence import (KM_PER_DEG_LAT, HURDAT2_PATH, STORMS,
                                  MODEL_CONTROL, km_offset, motion_kmh)

DT_H = 0.25
ATTENUATED_DRIFT = 0.42 * 2.49  # m/s — section-5.3 accounting, for P-S4


def advect(steer, lat0, lon0_signed, T):
    """RK2 tracer under annulus DLM sampled at tracer position. Signed lon (W<0)."""
    lat, lon = lat0, lon0_signed
    t = 0.0
    while t < T - 1e-9:
        dt = min(DT_H, T - t)
        u1, v1 = steer.get_dlm(t, lat, abs(lon))
        lat_m = lat + v1 * dt * 3600 / 2 / (KM_PER_DEG_LAT * 1000)
        lon_m = lon + u1 * dt * 3600 / 2 / (KM_PER_DEG_LAT * 1000 * math.cos(math.radians(lat)))
        u2, v2 = steer.get_dlm(t + dt / 2, lat_m, abs(lon_m))
        lat += v2 * dt * 3600 / (KM_PER_DEG_LAT * 1000)
        lon += u2 * dt * 3600 / (KM_PER_DEG_LAT * 1000 * math.cos(math.radians(lat)))
        t += dt
    return lat, lon


def main():
    storms = parse_hurdat2(HURDAT2_PATH)
    print("=" * 78)
    print("STEERING-ONLY TRACER — annulus DLM, no vortex, RK2 dt=0.25h")
    print("(along + = ahead of obs; cross + = right of obs motion; model = control run)")
    print("=" * 78)
    print(f"{'storm':<9} {'trc along':>10} {'trc cross':>10} {'trc total':>10} "
          f"{'m-t along':>10} {'|m - t|':>9} {'P-S4 ref':>9}")

    ok, rows = [], {}
    for name, sid, init, T in STORMS:
        tr = storms[sid]
        try:
            steer = ERA5Steering.load(storm=name.lower())
            f0 = tr.fix_at(init)
            plat, plon = advect(steer, f0.lat, f0.lon, T)
        except Exception as exc:            # domain exit, missing file, ...
            print(f"{name:<9} UNSCOREABLE: {exc}")
            continue

        landfall_dt = init + timedelta(hours=T)
        obs = tr.fix_at(landfall_dt)
        oe, on = motion_kmh(tr, landfall_dt, 3.0)
        norm = math.hypot(oe, on)
        ue, un = oe / norm, on / norm

        e_e, e_n = km_offset(plat, plon, obs.lat, obs.lon)
        along = e_e * ue + e_n * un
        cross = e_e * un - e_n * ue
        total = math.hypot(e_e, e_n)

        m_along, m_cross = MODEL_CONTROL[name]
        d_along, d_cross = m_along - along, m_cross - cross
        sep = math.hypot(d_along, d_cross)          # |model - tracer| at landfall
        ref = ATTENUATED_DRIFT * T * 3.6            # P-S4 reference footprint (km)

        ok.append(name)
        rows[name] = (along, cross, total, d_along, sep, ref)
        print(f"{name:<9} {along:>+10.1f} {cross:>+10.1f} {total:>10.1f} "
              f"{d_along:>+10.1f} {sep:>9.1f} {ref:>9.1f}")

    # -- scoring -------------------------------------------------------------
    print("-" * 78)
    print(f"P-S1 (scoreable 6/6): {'CONFIRMED' if len(ok) == 6 else f'FAILED ({len(ok)}/6)'}")
    if ok:
        cross_rms = math.sqrt(sum(rows[n][1] ** 2 for n in ok) / len(ok))
        model_cross_rms = math.sqrt(sum(MODEL_CONTROL[n][1] ** 2 for n in ok) / len(ok))
        print(f"tracer cross RMS {cross_rms:.1f} km vs model {model_cross_rms:.1f} km "
              f"(P-S2 threshold: <= {1.5 * model_cross_rms:.0f} km): "
              f"{'CONFIRMED' if cross_rms <= 1.5 * model_cross_rms else 'FAILED'}")
        ahead = [n for n in ok if rows[n][3] > 0]
        print(f"P-S3 (model ahead of tracer along-track on >=4/6): {len(ahead)}/6 "
              f"({', '.join(ahead)}) — {'CONFIRMED' if len(ahead) >= 4 else 'FAILED'}")
        gc = [n for n in ("Katrina", "Fran", "Michael", "Laura") if n in ok]
        within = [n for n in gc if 0.5 <= rows[n][4] / rows[n][5] <= 2.0]
        print(f"P-S4 (|model-tracer| within 2x of attenuated footprint, guard-clean): "
              f"{len(within)}/{len(gc)} ({', '.join(within)}) — "
              f"{'CONFIRMED' if len(within) >= 3 else 'FAILED'}")


if __name__ == "__main__":
    main()
