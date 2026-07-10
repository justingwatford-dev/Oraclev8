"""Paper-1 figures. Provenance: every number cites its committed record.
Run from repo root:  python figures/make_paper1_figures.py
"""
import numpy as np
import matplotlib.pyplot as plt
from fig_style import COMPACT, ENVELOPE, OBS, REF, save

# ---- F1: the compensating-errors cascade (paper-1 §3.2) ----------------------
CASCADE = [  # (flattering state, probe that broke it, what was exposed)
    ("Intensity 'signal'",        "honest-Vmax check",      "numerical runaway → cap"),
    ("Timing 'skill' (+8.9 h)",   "domain enlargement",     "boundary β-taper was the brake"),
    ("Contained vortex",          "R_env inertness test",   "wind never bounded → taper"),
    ("Hugo's 48-km landfall",     "structural re-run",      "β-drift × frozen-steering cancel"),
    ("Katrina's +14-km landfall", "HURDAT2 re-scoring",     "init offset × drift cancel"),
    ("The input data itself",     "primary-source audit",   "inits/tracks never HURDAT2"),
    ("13% 'over-translation'",    "Galilean control",       "±23-km tracker flicker"),
]

def f1_cascade():
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    ax.set_axis_off()
    n = len(CASCADE)
    for i, (state, probe, exposed) in enumerate(CASCADE):
        y = n - 1 - i
        ax.text(0.00, y, f"{i+1}. {state}", ha="left", va="center",
                fontsize=9, color=OBS, fontweight="bold")
        ax.annotate("", xy=(0.52, y), xytext=(0.40, y),
                    arrowprops=dict(arrowstyle="->", color=REF, lw=1.2))
        ax.text(0.46, y + 0.22, probe, ha="center", va="bottom",
                fontsize=7.8, color=REF, style="italic")
        ax.text(0.54, y, exposed, ha="left", va="center",
                fontsize=9, color=COMPACT)
        if i < n - 1:
            ax.annotate("", xy=(0.06, y - 0.72), xytext=(0.06, y - 0.28),
                        arrowprops=dict(arrowstyle="->", color=OBS, lw=1.0))
    ax.set_xlim(0, 1.02); ax.set_ylim(-0.6, n - 0.2)
    ax.set_title("The compensating-errors cascade: each repair exposed the "
                 "error above it had hidden", loc="left")
    save(fig, "figures/p1_f1_cascade")

# ---- F2: six storms — the cluster dissolves (paper-1 §4.2 + projection test) -
# landfall-fix cross/along, control config (PAPER_track_error_characterization /
# BETA_DRIFT_PROJECTION_TEST.md tables); projection-test predicted cross.
STORMS  = ["Hugo", "Katrina", "Ivan", "Fran", "Michael", "Laura"]
CROSS   = [110.2, 124.6, 126.3, 7.7, -98.7, -31.9]
ALONG   = [23.3, 76.5, 249.3, -45.5, 123.8, 37.4]
PREDX   = [140, 146, 198, 117, 86, 100]      # single-bias-vector projection
DISC    = [True, True, True, False, False, False]

def f2_sixstorm():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.4))
    for i, s in enumerate(STORMS):
        m = "o" if DISC[i] else "s"
        a.plot(CROSS[i], ALONG[i], m, ms=8, mfc=(OBS if DISC[i] else "white"),
               mec=OBS, mew=1.2)
        a.annotate(s, (CROSS[i], ALONG[i]), textcoords="offset points",
                   xytext=(6, 5), fontsize=7.8)
    a.axvline(0, color=REF, lw=0.8); a.axhline(0, color=REF, lw=0.8)
    a.axvspan(95, 140, color=COMPACT, alpha=0.12)
    a.text(117, -80, "the\n'cluster'", ha="center", fontsize=7.8, color=COMPACT)
    a.set_xlabel("cross-track error at landfall fix (km, + = east)")
    a.set_ylabel("along-track error (km, + = ahead)")
    a.set_title("(a) discovery (●) vs test (□) storms", loc="left")
    for i, s in enumerate(STORMS):
        m = "o" if DISC[i] else "s"
        b.plot(PREDX[i], CROSS[i], m, ms=8, mfc=(OBS if DISC[i] else "white"),
               mec=OBS, mew=1.2)
        b.annotate(s, (PREDX[i], CROSS[i]), textcoords="offset points",
                   xytext=(5, 4), fontsize=7.8)
    lim = [-130, 220]
    b.plot(lim, lim, "--", color=REF, lw=1.0, label="perfect attribution")
    b.axhline(0, color=REF, lw=0.8)
    b.set_xlim(60, 220); b.set_ylim(*lim)
    b.set_xlabel("cross-track predicted by β-bias projection (km)")
    b.set_ylabel("observed cross-track (km)")
    b.set_title("(b) the bias predicts east for all six;\n"
                "the test storms miss west — bridge falsified", loc="left")
    b.legend(loc="lower right")
    fig.tight_layout(); save(fig, "figures/p1_f2_sixstorm")

# ---- F3: β-drift characterization (paper-1 §4.1) ------------------------------
# compact mature vector 2.49 m/s @ 350° (gate-beta-shape control row);
# canonical band 1.5–2.5 m/s toward 290–335° (paper-1 Table 1 refs);
# precession series: gate-beta-longrun clean windows t6–60: 337/342/347/351/356.
def f3_betadrift():
    fig = plt.figure(figsize=(7.2, 3.3))
    a = fig.add_subplot(121, projection="polar")
    a.set_theta_zero_location("N"); a.set_theta_direction(-1)
    th = np.radians(np.linspace(290, 335, 60))
    a.fill_between(th, 1.5, 2.5, color=ENVELOPE, alpha=0.18)
    a.text(np.radians(312), 2.9, "canonical\nβ-drift", ha="center",
           fontsize=8, color=ENVELOPE)
    a.annotate("", xy=(np.radians(350), 2.49), xytext=(0, 0),
               arrowprops=dict(arrowstyle="-|>", color=COMPACT, lw=2.2))
    a.text(np.radians(355), 2.75, "Oracle\n2.49 @ 350°", ha="center",
           fontsize=8, color=COMPACT)
    a.set_ylim(0, 3.2); a.set_yticks([1, 2, 3])
    a.set_thetamin(250); a.set_thetamax(40)
    a.set_title("(a) mature β-drift vector", loc="left", pad=18)
    b = fig.add_subplot(122)
    tw = [12, 23, 34, 44, 55]                       # window centres, t6–60 h
    hd = [337, 342, 347, 351, 356]
    b.plot(tw, hd, "-o", color=COMPACT, ms=5, label="Oracle (heading)")
    b.axhspan(290, 335, color=ENVELOPE, alpha=0.18)
    b.text(46, 312, "canonical band", fontsize=8, color=ENVELOPE)
    b.axhline(360, color=REF, lw=0.8, ls=":")
    b.text(8, 361, "due north", fontsize=7.5, color=REF)
    b.set_xlabel("time (h)"); b.set_ylabel("drift heading (° toward)")
    b.set_ylim(285, 372)
    b.set_title("(b) the heading never locks:\n~0.4° h⁻¹ through north",
                loc="left")
    fig.tight_layout(); save(fig, "figures/p1_f3_betadrift")

if __name__ == "__main__":
    f1_cascade(); f2_sixstorm(); f3_betadrift()
