"""
Oracle V8 — persistence baseline for the six-storm landfall record
==================================================================
Registered analysis 1 of PAPER2_TIGHTENING_predictions.md (P-B1..P-B3).

For each of the six storms: take the observed HURDAT2 motion at init
(centered difference over init +/- 6 h), extrapolate linearly for the
storm's transit time, and score the persistence forecast at the observed
landfall fix with the manuscripts' along/cross conventions
(along + = ahead of obs; cross + = right of obs motion).

No model runs. Inputs: oracle_v8/hurdat2.txt only.

Run:  python -m oracle_v8.measure_persistence
"""
from __future__ import annotations

import math
import os
from datetime import datetime, timedelta

from .hurdat2 import parse_hurdat2

KM_PER_DEG_LAT = 111.32
HERE = os.path.dirname(__file__)
HURDAT2_PATH = os.path.join(HERE, "hurdat2.txt")

# (storm, hurdat2 id, init, transit hours) — from the checked-in run-log headers
STORMS = [
    ("Hugo",    "AL111989", datetime(1989, 9, 21, 0),  28.00),
    ("Katrina", "AL122005", datetime(2005, 8, 28, 0),  35.17),
    ("Ivan",    "AL092004", datetime(2004, 9, 14, 12), 42.83),
    ("Fran",    "AL061996", datetime(1996, 9, 5, 0),   24.50),
    ("Michael", "AL142018", datetime(2018, 10, 9, 12), 29.50),
    ("Laura",   "AL132020", datetime(2020, 8, 26, 6),  24.00),
]

# Model control (compact taper) landfall-fix decomposition — committed Table 2 values
MODEL_CONTROL = {
    "Hugo":    (+23.3, +110.2),
    "Katrina": (+76.5, +124.6),
    "Ivan":    (+249.3, +126.3),
    "Fran":    (-45.5, +7.7),
    "Michael": (+123.8, -98.7),
    "Laura":   (+37.4, -31.9),
}


def km_offset(lat, lon, lat_ref, lon_ref):
    """(east_km, north_km) of (lat, lon) relative to (lat_ref, lon_ref)."""
    east = (lon - lon_ref) * KM_PER_DEG_LAT * math.cos(math.radians(lat_ref))
    north = (lat - lat_ref) * KM_PER_DEG_LAT
    return east, north


def motion_kmh(track, dt_center, half_window_h):
    """Centered-difference observed motion (east, north) km/h around dt_center."""
    a = track.fix_at(dt_center - timedelta(hours=half_window_h))
    b = track.fix_at(dt_center + timedelta(hours=half_window_h))
    e, n = km_offset(b.lat, b.lon, a.lat, a.lon)
    return e / (2 * half_window_h), n / (2 * half_window_h)


def main():
    storms = parse_hurdat2(HURDAT2_PATH)
    print("=" * 78)
    print("PERSISTENCE BASELINE — observed motion at init, extrapolated to landfall")
    print("(along + = ahead of obs; cross + = right of obs motion; model = control run)")
    print("=" * 78)
    print(f"{'storm':<9} {'T(h)':>6} {'pers along':>11} {'pers cross':>11} "
          f"{'pers total':>11} {'model total':>12}")

    rows = {}
    for name, sid, init, T in STORMS:
        tr = storms[sid]
        landfall_dt = init + timedelta(hours=T)

        # persistence: init position + init motion * T
        f0 = tr.fix_at(init)
        me, mn = motion_kmh(tr, init, 6.0)
        pred_e_total, pred_n_total = me * T, mn * T
        pred_lat = f0.lat + pred_n_total / KM_PER_DEG_LAT
        pred_lon = f0.lon + pred_e_total / (KM_PER_DEG_LAT * math.cos(math.radians(f0.lat)))

        # observed landfall fix + observed motion there
        obs = tr.fix_at(landfall_dt)
        oe, on = motion_kmh(tr, landfall_dt, 3.0)
        norm = math.hypot(oe, on)
        ue, un = oe / norm, on / norm

        # error vector, persistence minus observed, decomposed
        e_e, e_n = km_offset(pred_lat, pred_lon, obs.lat, obs.lon)
        along = e_e * ue + e_n * un
        cross = e_e * un - e_n * ue
        total = math.hypot(e_e, e_n)

        m_along, m_cross = MODEL_CONTROL[name]
        m_total = math.hypot(m_along, m_cross)
        rows[name] = (along, cross, total, m_total)
        print(f"{name:<9} {T:>6.2f} {along:>+11.1f} {cross:>+11.1f} "
              f"{total:>11.1f} {m_total:>12.1f}")

    # -- scoring -------------------------------------------------------------
    crosses = [rows[n][1] for n, *_ in [(s[0],) for s in STORMS]]
    cross_rms = math.sqrt(sum(c * c for c in crosses) / len(crosses))
    model_cross_rms = math.sqrt(sum(MODEL_CONTROL[n][1] ** 2 for n in rows) / len(rows))
    beats = [n for n in rows if rows[n][2] < rows[n][3]]
    recurver_totals = {n: rows[n][2] for n in ("Ivan", "Michael")}

    print("-" * 78)
    print(f"persistence cross-track RMS: {cross_rms:.1f} km   "
          f"(model control: {model_cross_rms:.1f} km; P-B1 threshold >190 km)")
    print(f"P-B1: {'CONFIRMED' if cross_rms > 190 else 'FAILED'}")
    print(f"P-B2 (persistence total > model total on all six): "
          f"{'CONFIRMED' if not beats else 'FAILED — persistence beats model on: ' + ', '.join(beats)}")
    print(f"P-B3 (Ivan & Michael persistence totals > 300 km): "
          f"Ivan {recurver_totals['Ivan']:.0f}, Michael {recurver_totals['Michael']:.0f} — "
          f"{'CONFIRMED' if all(v > 300 for v in recurver_totals.values()) else 'FAILED'}")


if __name__ == "__main__":
    main()
