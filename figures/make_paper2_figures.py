"""Paper-2 figures. Provenance: profiles computed live from vortex_init;
all other numbers cite OVERROTATION_CANDIDATES.md / ENVELOPE_INTENSIFICATION.md
tables and the committed run logs. Run from repo root:
    python figures/make_paper2_figures.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("ORACLE_GPU", "0")
import numpy as np
import matplotlib.pyplot as plt
from fig_style import COMPACT, ENVELOPE, OBS, REF, save
from oracle_v8.vortex_init import HollandVortexInit

# ---- F1: profiles — the cutoff vs the tail (computed from vortex_init) --------
def f1_profiles():
    r = np.linspace(1e3, 1000e3, 4000)
    kw = dict(Vmax=64.0, Rmax=75e3, B=1.5, f=5.7e-5, R_env=500e3)
    vc = HollandVortexInit(**kw, wind_taper=True,
                           taper_start_frac=0.40).tangential_wind(r)
    ve = HollandVortexInit(**kw, outer_envelope_m=420e3).tangential_wind(r)
    def zeta(v):                                    # (1/r) d(rV)/dr
        return np.gradient(r * v, r) / r
    fig, (a, b) = plt.subplots(2, 1, figsize=(5.4, 4.6), sharex=True)
    a.plot(r/1e3, vc, color=COMPACT, label="compact taper (200–500 km)")
    a.plot(r/1e3, ve, color=ENVELOPE, label="Gaussian envelope (r_d = 420 km)")
    a.axvspan(200, 500, color=REF, alpha=0.10)
    a.set_ylabel("tangential wind (m s$^{-1}$)"); a.legend()
    a.set_title("(a) the envelope keeps the outer tail the taper removes",
                loc="left")
    b.plot(r/1e3, zeta(vc)*1e4, color=COMPACT)
    b.plot(r/1e3, zeta(ve)*1e4, color=ENVELOPE)
    b.axhline(0, color=REF, lw=0.8)
    b.axvspan(200, 500, color=REF, alpha=0.10)
    b.annotate("taper's anticyclonic ring", xy=(330, -1.1), fontsize=8,
               color=COMPACT)
    b.set_xlabel("radius (km)"); b.set_ylabel(r"$\zeta$ (10$^{-4}$ s$^{-1}$)")
    b.set_xlim(0, 1000); b.set_ylim(-1.6, 3.0)
    b.set_title("(b) implied relative vorticity", loc="left")
    fig.tight_layout(); save(fig, "figures/p2_f1_profiles")

# ---- F2: phase lock vs free precession ----------------------------------------
# gauss traces: OVERROTATION_CANDIDATES.md phase-lock append (t12/24/36/48);
# compact precession: gate-beta-longrun clean windows (t6–60).
def f2_phaselock():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2))
    t4 = [12, 24, 36, 48]
    a.plot([12, 23, 34, 44, 55], [337, 342, 347, 351, 356], "-o",
           color=COMPACT, ms=5, label="compact (free precession)")
    a.plot(t4, [328, 327, 328, 329], "-o", color=ENVELOPE, ms=5,
           label="envelope r_d=420 (locked)")
    a.plot(t4, [329, 325, 325, 326], "-o", color=ENVELOPE, ms=5, mfc="white",
           label="envelope r_d=560")
    a.axhspan(290, 335, color=ENVELOPE, alpha=0.12)
    a.axhline(360, color=REF, lw=0.8, ls=":")
    a.set_xlabel("time (h)"); a.set_ylabel("gyre/drift heading (° toward)")
    a.set_ylim(285, 372); a.legend(loc="upper left", fontsize=7.6)
    a.set_title("(a) orientation: locked by t = 12 h,\nor never", loc="left")
    b.plot(t4, [1.11, 1.73, 2.20, 2.37], "-o", color=ENVELOPE, ms=5)
    b.plot(t4, [1.50, 2.25, 2.83, 2.81], "-o", color=ENVELOPE, ms=5, mfc="white")
    b.axhline(2.49, color=COMPACT, lw=1.2, ls="--")
    b.text(13, 2.55, "compact mature speed", fontsize=7.8, color=COMPACT)
    b.set_xlabel("time (h)"); b.set_ylabel("drift speed (m s$^{-1}$)")
    b.set_ylim(0, 3.2)
    b.set_title("(b) amplitude grows into a\nfixed orientation", loc="left")
    fig.tight_layout(); save(fig, "figures/p2_f2_phaselock")

# ---- F3: six-storm transmission (Stage-3 table) --------------------------------
STORMS = ["Hugo", "Katrina", "Ivan", "Fran", "Michael", "Laura"]
PRED   = [-93, -108, -139, -81, -73, -74]     # strong-form Δcross
OBSD   = [-31.5, -46.4, -86.1, -37.6, -14.5, -19.9]
GUARD  = [False, True, False, True, True, True]   # True = intensity-guard clean

def f3_transmission():
    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    lim = [-150, 10]
    ax.plot(lim, lim, "--", color=REF, lw=1.0, label="1:1 (linear projection)")
    x = np.array(lim)
    ax.plot(x, 0.34 * x, "-", color=ENVELOPE, lw=1.4,
            label="fitted transmission ≈ 0.34")
    for i, s in enumerate(STORMS):
        ax.plot(PRED[i], OBSD[i], "o", ms=8, mec=OBS, mew=1.2,
                mfc=(OBS if GUARD[i] else "white"))
        ax.annotate(s, (PRED[i], OBSD[i]), textcoords="offset points",
                    xytext=(6, 5), fontsize=7.8)
    ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect("equal")
    ax.set_xlabel("predicted westward shift, linear projection (km)")
    ax.set_ylabel("observed westward shift (km)")
    ax.set_title("Six for six on sign, one-third on size\n"
                  "(open = intensity guard exceeded)", loc="left")
    ax.legend(loc="upper left", fontsize=7.8)
    fig.tight_layout(); save(fig, "figures/p2_f3_transmission")

# ---- F4: the fix's price — delay, and the Ekman mechanism ---------------------
# panel a: gate-j2-profile Vmax traces (Run-3 log, t = 4..52 h);
# panel b: Run-4 moving-frame budget minBL300, pre-onset (AM-budget_gh.txt).
def f4_delay():
    t = np.arange(4, 53, 4)
    vc = [54.5, 52.0, 50.1, 48.1, 48.9, 49.7, 63.3, 78.2, 77.1, 79.9, 73.8,
          82.5, 75.6]
    ve = [53.0, 50.1, 47.6, 46.6, 46.3, 45.4, 45.0, 47.8, 49.1, 66.2, 70.5,
          74.1, 72.9]
    tb = np.arange(2, 25, 2)
    mc = [-7.24e7, -4.04e7, -1.22e7, 1.53e6, 1.40e7, 2.33e7, 2.80e7, 3.21e7,
          3.79e7, 4.25e7, 4.43e7, 4.77e7]
    me = [-4.99e7, -3.13e7, -1.78e7, -7.98e6, -2.96e6, 2.91e6, 4.99e6,
          8.96e6, 1.22e7, 1.48e7, 1.68e7, 1.83e7]
    fig, (a, b) = plt.subplots(2, 1, figsize=(5.6, 4.9),
                               gridspec_kw={"height_ratios": [1.5, 1]})
    a.plot(t, vc, "-o", color=COMPACT, ms=4, label="compact taper")
    a.plot(t, ve, "-o", color=ENVELOPE, ms=4, label="envelope r_d=420")
    a.annotate("onset 24 h", xy=(26, 62), fontsize=8, color=COMPACT)
    a.annotate("onset 32 h", xy=(35, 55), fontsize=8, color=ENVELOPE)
    a.set_ylabel("V$_{max}$ (m s$^{-1}$)"); a.legend(loc="upper left")
    a.set_title("(a) delayed, not damped: same equilibrium, 8 h later",
                loc="left")
    b.plot(tb, np.array(mc)/1e7, "-o", color=COMPACT, ms=4)
    b.plot(tb, np.array(me)/1e7, "-o", color=ENVELOPE, ms=4)
    b.axhline(0, color=REF, lw=0.8)
    b.set_xlabel("time (h)")
    b.set_ylabel("BL mass inflow at 300 km\n(10$^{7}$ kg s$^{-1}$)")
    b.set_title("(b) the mechanism: mid-radius Ekman inflow,\n"
                "×2.5–3 stronger under the compact profile", loc="left")
    fig.tight_layout(); save(fig, "figures/p2_f4_delay")

if __name__ == "__main__":
    f1_profiles(); f2_phaselock(); f3_transmission(); f4_delay()
