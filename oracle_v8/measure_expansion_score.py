"""
Oracle V8 — expansion A/B scoring (Charley, Florence, Ida)
==========================================================
Scores P-N1..P-N5 of PAPER2_EXPANSION_predictions.md from the six GPU run
logs. Observed landfall-fix decompositions are transcribed from the runs'
landfall_report blocks (provenance: oracle_v8/Logs/{Storm}/*.txt, this
commit). The per-storm decomposition (P-N4) reuses the registered
machinery of measure_transmission_decomp verbatim.

Run:  python -m oracle_v8.measure_expansion_score
"""
from __future__ import annotations

import math

import numpy as np

from .measure_transmission_decomp import (env_vec, cmp_vec, spinup,
                                          vmax_history, find_log,
                                          DELTA_MATURE, cross_proj)

# storm -> (registered heading, transit h, dominant axis)
EXP = {
    "Charley":  (15.0,  25.75, "cross"),
    "Florence": (290.0, 35.25, "along"),
    "Ida":      (330.0, 28.92, "cross"),
}

# landfall-fix decomposition (along, cross) km, from the run logs
OBS = {
    "Charley":  {"control": (-160.7, -31.8), "envelope": (-172.8, -45.0)},
    "Florence": {"control": (-124.3, +62.7), "envelope": (-78.3, +65.7)},
    "Ida":      {"control": (-14.0, +36.1),  "envelope": (+0.6, +13.1)},
}

# same-latitude crossing timing (h, model minus obs) per arm, from the logs
TIMING = {
    "Charley":  {"control": +4.4, "envelope": +5.2},
    "Florence": {"control": -3.4, "envelope": -3.0},
    "Ida":      {"control": +0.1, "envelope": +0.1},
}


def along_proj(vec, heading_deg):
    th = math.radians(heading_deg)
    return vec[0] * math.sin(th) + vec[1] * math.cos(th)


def main():
    print("=" * 78)
    print("EXPANSION A/B SCORING — Charley, Florence, Ida (GPU runs)")
    print("=" * 78)
    print(f"{'storm':<9} {'axis':<6} {'obs dAlong':>10} {'obs dCross':>10} "
          f"{'strong':>8} {'obs ratio':>9} {'pred ratio':>10} {'dVmax':>6} {'dT(h)':>6}")

    results = {}
    for storm, (hdg, T, axis) in EXP.items():
        c_along, c_cross = OBS[storm]["control"]
        e_along, e_cross = OBS[storm]["envelope"]
        d_along, d_cross = e_along - c_along, e_cross - c_cross

        proj = cross_proj if axis == "cross" else along_proj
        strong = proj(DELTA_MATURE, hdg) * T * 3.6
        obs_shift = d_cross if axis == "cross" else d_along
        obs_ratio = obs_shift / strong

        # per-storm decomposition prediction on the dominant axis
        t_c, v_c = vmax_history(find_log(storm, ["*_Agnostic*.txt", "*_agnostic*.txt"]))
        t_e, v_e = vmax_history(find_log(storm, ["*auss_envelope*.txt"]))
        acc = 0.0
        for t in np.arange(0.0, T + 1e-9, 0.5):
            Vc = float(np.interp(t, t_c, v_c))
            Ve = float(np.interp(t, t_e, v_e))
            d = (env_vec(Ve)[0] - cmp_vec(Vc)[0], env_vec(Ve)[1] - cmp_vec(Vc)[1])
            acc += spinup(t) * proj(d, hdg) * 0.5 * 3.6
        pred_ratio = acc / strong

        # guards: max |Vmax difference| between arms over the transit; A/B timing
        grid = np.arange(0.0, T + 1e-9, 0.5)
        dv = max(abs(float(np.interp(t, t_c, v_c)) - float(np.interp(t, t_e, v_e)))
                 for t in grid)
        dt_ab = abs(TIMING[storm]["envelope"] - TIMING[storm]["control"])
        results[storm] = dict(d_along=d_along, d_cross=d_cross, strong=strong,
                              obs_ratio=obs_ratio, pred_ratio=pred_ratio,
                              dv=dv, dt=dt_ab, axis=axis)
        print(f"{storm:<9} {axis:<6} {d_along:>+10.1f} {d_cross:>+10.1f} "
              f"{strong:>+8.0f} {obs_ratio:>9.2f} {pred_ratio:>10.2f} "
              f"{dv:>6.1f} {dt_ab:>6.1f}")

    # -- scoring -------------------------------------------------------------
    print("-" * 78)
    clean = [s for s in results if results[s]["dv"] <= 10.0 and results[s]["dt"] <= 3.0]
    print(f"guards: clean = {clean} "
          f"(intensity <=10 m/s and A/B timing <=3 h)")

    fix_cross = {"Charley": (-31.8, -45.0), "Florence": (62.7, 65.7),
                 "Ida": (36.1, 13.1)}
    ok150 = [s for s, (c, e) in fix_cross.items() if abs(c) <= 150 and abs(e) <= 150]
    print(f"P-N1 (scoreable, >=2/3 with |cross|<=150 both arms): {len(ok150)}/3 "
          f"({', '.join(ok150)}) — {'CONFIRMED' if len(ok150) >= 2 else 'FAILED'}")

    sign_ok = []
    for s in results:
        r = results[s]
        want_neg = r["axis"] == "cross"          # west (Charley, Ida)
        obs = r["d_cross"] if r["axis"] == "cross" else r["d_along"]
        expected_sign = -1 if want_neg else +1   # Florence: along FORWARD (+)
        if obs * expected_sign > 0:
            sign_ok.append(s)
    print(f"P-N2 (dominant-axis sign, all three): {len(sign_ok)}/3 "
          f"({', '.join(sign_ok)}) — {'CONFIRMED' if len(sign_ok) == 3 else 'FAILED'}")

    in_band = {s: 0.15 <= results[s]["obs_ratio"] <= 0.55 for s in clean}
    print(f"P-N3 (guard-clean ratios in [0.15, 0.55]): "
          + ", ".join(f"{s} {results[s]['obs_ratio']:.2f}" for s in clean)
          + f" — {'CONFIRMED' if all(in_band.values()) else 'FAILED'}")

    close = {s: abs(results[s]["pred_ratio"] - results[s]["obs_ratio"]) <= 0.15
             for s in clean}
    print(f"P-N4 (decomp within +/-0.15, guard-clean, out of sample): "
          + ", ".join(f"{s} pred {results[s]['pred_ratio']:.2f} obs "
                      f"{results[s]['obs_ratio']:.2f}" for s in clean)
          + f" — {'CONFIRMED' if all(close.values()) else 'FAILED'}")

    fl = results["Florence"]
    print(f"P-N5 (Florence |along| > |cross|): |{fl['d_along']:.1f}| vs "
          f"|{fl['d_cross']:.1f}| — "
          f"{'CONFIRMED' if abs(fl['d_along']) > abs(fl['d_cross']) else 'FAILED'}")


if __name__ == "__main__":
    main()
