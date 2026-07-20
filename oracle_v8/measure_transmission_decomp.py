"""
Oracle V8 — per-storm transmission decomposition
================================================
Registered analysis 3 of PAPER2_TIGHTENING_predictions.md (P-D1..P-D4).

The section-5.3 aggregate accounting (intensity scaling x gyre spin-up
~ 0.42) made per-storm: family drift-vs-intensity laws from the committed
testbed tables, the committed envelope spin-up curve s(t), and each run's
own Vmax'(t) history parsed from its checked-in log. Predicted per-storm
shift = transit integral of the landfall-heading cross-projection of
s(t) * [env_vec(V_env(t)) - cmp_vec(V_cmp(t))]; predicted ratio = that
integral over the strong-form (mature-Delta x T) projection.

No model runs. Inputs: committed tables + checked-in run logs.

Run:  python -m oracle_v8.measure_transmission_decomp
"""
from __future__ import annotations

import glob
import math
import os
import re

import numpy as np

HERE = os.path.dirname(__file__)
LOGS = os.path.join(HERE, "Logs")

# storm -> (heading deg cw from N, transit h, observed ratio, guard-clean?)
STORMS = {
    "Hugo":    (325.0, 28.00, 0.34, False),
    "Katrina": (350.0, 35.17, 0.43, True),
    "Ivan":    (340.0, 42.83, 0.62, False),
    "Fran":    (335.0, 24.50, 0.46, True),
    "Michael": (10.0,  29.50, 0.20, True),
    "Laura":   (350.0, 24.00, 0.27, True),
}

# ---- family drift laws (committed testbed tables) --------------------------
# compact family: west ~ constant 0.41 across V; north proportional to V
#   anchored at the production control (V_end 42.2 -> north 2.45).
# envelope family (r_d = 420 ladder): (V_end, drift, west): least-squares
#   linear fits in V for west and north (north = sqrt(drift^2 - west^2)).
ENV_LADDER = [(39.7, 2.30, 1.20), (25.7, 1.63, 0.92), (17.2, 1.23, 0.72)]
V_CLAMP = (15.0, 45.0)


def _linfit(xs, ys):
    A = np.vstack([xs, np.ones(len(xs))]).T
    m, b = np.linalg.lstsq(A, ys, rcond=None)[0]
    return float(m), float(b)


_env_V = [v for v, _, _ in ENV_LADDER]
_env_w = [w for _, _, w in ENV_LADDER]
_env_n = [math.sqrt(d * d - w * w) for _, d, w in ENV_LADDER]
W_M, W_B = _linfit(_env_V, _env_w)
N_M, N_B = _linfit(_env_V, _env_n)


def cmp_vec(V):
    V = min(max(V, *[V_CLAMP[0]]), V_CLAMP[1])
    return (-0.41, (2.45 / 42.2) * V)          # (east, north): west 0.41 -> east -0.41


def env_vec(V):
    V = min(max(V, V_CLAMP[0]), V_CLAMP[1])
    return (-(W_M * V + W_B), N_M * V + N_B)


# mature Delta reproduces the committed (-0.78, -0.49) at V_end 39.7 / 42.2
DELTA_MATURE = (env_vec(39.7)[0] - cmp_vec(42.2)[0],
                env_vec(39.7)[1] - cmp_vec(42.2)[1])

# ---- spin-up curve (committed envelope trace, r_d = 420) -------------------
S_T = [0.0, 12.0, 24.0, 36.0, 48.0]
S_V = [0.0, 1.11 / 2.37, 1.73 / 2.37, 2.20 / 2.37, 1.0]


def spinup(t):
    return float(np.interp(t, S_T, S_V))


# ---- Vmax'(t) parsing from run logs ----------------------------------------
TIME_RE = re.compile(r"^\s+(\d+\.\d)\s+\d+\.\d+\s+\d+\.\d+\s")
VMAX_RE = re.compile(r"Vmax'=\s*([\d.]+)@")


def vmax_history(log_path):
    ts, vs, t_cur = [], [], None
    with open(log_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = TIME_RE.match(line)
            if m:
                t_cur = float(m.group(1))
                continue
            m = VMAX_RE.search(line)
            if m and t_cur is not None:
                ts.append(t_cur)
                vs.append(float(m.group(1)))
    if not ts:
        raise ValueError(f"no Vmax' history in {log_path}")
    return np.array(ts), np.array(vs)


def find_log(storm, patterns):
    for pat in patterns:
        hits = glob.glob(os.path.join(LOGS, storm, pat))
        if hits:
            return sorted(hits)[0]
    raise FileNotFoundError(f"{storm}: none of {patterns}")


def cross_proj(vec, heading_deg):
    th = math.radians(heading_deg)
    return vec[0] * math.cos(th) - vec[1] * math.sin(th)


def main():
    print("=" * 78)
    print("PER-STORM TRANSMISSION DECOMPOSITION — intensity history x gyre spin-up")
    print(f"(mature Delta = ({DELTA_MATURE[0]:+.2f} E, {DELTA_MATURE[1]:+.2f} N) m/s; "
          "committed: (-0.78, -0.49))")
    print("=" * 78)
    print(f"{'storm':<9} {'pred shift':>10} {'strong':>8} {'pred ratio':>10} "
          f"{'obs ratio':>9}  flag")

    pred = {}
    for storm, (hdg, T, obs_ratio, clean) in STORMS.items():
        t_c, v_c = vmax_history(find_log(storm, ["*_Agnostic*.txt", "*_agnostic*.txt"]))
        t_e, v_e = vmax_history(find_log(storm, ["*auss_envelope*.txt"]))

        dt = 0.5
        grid = np.arange(0.0, T + 1e-9, dt)
        shift = 0.0
        for t in grid:
            Vc = float(np.interp(t, t_c, v_c))
            Ve = float(np.interp(t, t_e, v_e))
            d = (env_vec(Ve)[0] - cmp_vec(Vc)[0], env_vec(Ve)[1] - cmp_vec(Vc)[1])
            shift += spinup(t) * cross_proj(d, hdg) * dt * 3.6
        strong = cross_proj(DELTA_MATURE, hdg) * T * 3.6
        ratio = shift / strong
        pred[storm] = ratio
        print(f"{storm:<9} {shift:>+10.1f} {strong:>+8.1f} {ratio:>10.2f} "
              f"{obs_ratio:>9.2f}  {'' if clean else 'guard-flagged'}")

    # -- scoring -------------------------------------------------------------
    print("-" * 78)
    gc = {s: v for s, (h, t, o, c) in STORMS.items() if c for v in [pred[s]]}
    obs = {s: STORMS[s][2] for s in gc}
    top_pred = set(sorted(gc, key=gc.get, reverse=True)[:2])
    print(f"P-D1 (predicted top two of guard-clean = {{Katrina, Fran}}): "
          f"predicted {top_pred} — "
          f"{'CONFIRMED' if top_pred == {'Katrina', 'Fran'} else 'FAILED'}")
    within = {s: abs(gc[s] - obs[s]) <= 0.15 for s in gc}
    print(f"P-D2 (all four within +/-0.15): "
          + ", ".join(f"{s} {gc[s] - obs[s]:+.2f}" for s in gc)
          + f" — {'CONFIRMED' if all(within.values()) else 'FAILED'}")
    ivan_short = STORMS['Ivan'][2] - pred['Ivan']
    print(f"P-D3 (Ivan predicted short of 0.62 by >0.15): shortfall {ivan_short:+.2f} — "
          f"{'CONFIRMED' if ivan_short > 0.15 else 'FAILED'}")
    mean_gc = sum(gc.values()) / len(gc)
    print(f"P-D4 (guard-clean mean within +/-0.10 of 0.42): mean {mean_gc:.2f} — "
          f"{'CONFIRMED' if abs(mean_gc - 0.42) <= 0.10 else 'FAILED'}")

    # -- POST-HOC (not registered): along-axis transmission where the observed
    #    along-track shift is committed (section 5.2: Katrina 76.5->61.4,
    #    Michael 123.8->93.4, Laura 37.4->30.0). Tests whether a cross-axis
    #    outlier is axis-specific. Exploratory; label as such wherever used.
    print("-" * 78)
    print("POST-HOC along-axis transmission (observed / strong-form, same Delta):")
    OBS_ALONG = {"Katrina": 61.4 - 76.5, "Michael": 93.4 - 123.8, "Laura": 30.0 - 37.4}
    for s, d_obs in OBS_ALONG.items():
        hdg, T, _, _ = STORMS[s]
        th = math.radians(hdg)
        a_vel = DELTA_MATURE[0] * math.sin(th) + DELTA_MATURE[1] * math.cos(th)
        strong_along = a_vel * T * 3.6
        # predicted along ratio with the same intensity-history x spin-up integral
        t_c, v_c = vmax_history(find_log(s, ["*_Agnostic*.txt", "*_agnostic*.txt"]))
        t_e, v_e = vmax_history(find_log(s, ["*auss_envelope*.txt"]))
        acc = 0.0
        for t in np.arange(0.0, T + 1e-9, 0.5):
            Vc = float(np.interp(t, t_c, v_c))
            Ve = float(np.interp(t, t_e, v_e))
            d = (env_vec(Ve)[0] - cmp_vec(Vc)[0], env_vec(Ve)[1] - cmp_vec(Vc)[1])
            acc += spinup(t) * (d[0] * math.sin(th) + d[1] * math.cos(th)) * 0.5 * 3.6
        print(f"  {s:<9} obs {d_obs / strong_along:>5.2f}   pred {acc / strong_along:>5.2f}   "
              f"(strong-form along {strong_along:+.0f} km, observed {d_obs:+.1f} km)")


if __name__ == "__main__":
    main()
