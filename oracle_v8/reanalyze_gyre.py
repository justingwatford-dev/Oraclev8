#!/usr/bin/env python3
"""
Re-extract beta-gyre diagnostics from saved snapshots.

The workspace can contain several snapshot families at once, for example:
  gyre_snap_t12h.npz          legacy untagged set
  gyre_snap_v64_t12h.npz      Vmax-tagged set

Pass --tag explicitly when more than one set is present:
  python oracle_v8/reanalyze_gyre.py --tag v64
  python oracle_v8/reanalyze_gyre.py --tag untagged
"""
import argparse
import glob
import math
import re

import numpy as np


def _args():
    p = argparse.ArgumentParser(description="Re-analyze saved beta-gyre snapshots.")
    p.add_argument(
        "--tag",
        default=None,
        help=(
            "Snapshot tag to load, e.g. v64 for gyre_snap_v64_t*.npz. "
            "Use 'untagged' for legacy gyre_snap_t*.npz."
        ),
    )
    return p.parse_args()


def _snapshot_set(tag):
    if tag in ("untagged", "legacy", "base"):
        return "untagged", "gyre_snap_t*.npz"
    if tag:
        return tag, f"gyre_snap_{tag}_t*.npz"

    tagged = sorted(glob.glob("gyre_snap_v*_t*.npz"))
    untagged = sorted(glob.glob("gyre_snap_t*.npz"))
    if tagged and untagged:
        tags = sorted({re.search(r"gyre_snap_(v\d+)_t", f).group(1) for f in tagged})
        raise SystemExit(
            "multiple gyre snapshot sets are present; pass --tag "
            + "|".join(tags + ["untagged"])
        )
    if tagged:
        tags = sorted({re.search(r"gyre_snap_(v\d+)_t", f).group(1) for f in tagged})
        if len(tags) > 1:
            raise SystemExit(
                "multiple tagged gyre snapshot sets are present; pass --tag "
                + "|".join(tags)
            )
        return tags[0], f"gyre_snap_{tags[0]}_t*.npz"
    return "untagged", "gyre_snap_t*.npz"


ARGS = _args()
SNAP_TAG, SNAP_PATTERN = _snapshot_set(ARGS.tag)
NPZ = sorted(
    glob.glob(SNAP_PATTERN),
    key=lambda f: int(re.search(r"t(\d+)h", f).group(1)),
)
if not NPZ:
    raise SystemExit(f"no {SNAP_PATTERN} in cwd - run gate-beta-gyre first")
OUT_PNG = (
    "gyre_precession.png"
    if SNAP_TAG == "untagged"
    else f"gyre_precession_{SNAP_TAG}.png"
)

BANDS_KM = [(75, 450), (150, 450), (250, 500), (150, 350)]
PROFILE_R_KM = [50, 100, 150, 200, 250, 300, 350, 400, 450]
CORE_WIN_KM = 120.0
DISK_KM = 200.0


def _grid(a, dx):
    nx, _ny = a.shape
    xs = (np.arange(nx) + 0.5) * dx
    return np.meshgrid(xs, xs, indexing="ij")


def to_compass(math_deg):
    return (90.0 - math_deg) % 360.0


def crate(a, b, dt):
    return ((((b - a + 180) % 360) - 180) / dt)


def refine_center(zeta, cx, cy, dx, win=CORE_WIN_KM * 1e3):
    X, Y = _grid(zeta, dx)
    R = np.hypot(X - cx, Y - cy)
    w = np.where((R < win) & (zeta > 0), zeta, 0.0)
    tot = w.sum()
    if tot <= 0:
        return cx, cy
    return float((w * X).sum() / tot), float((w * Y).sum() / tot)


def m1(zeta, cx, cy, dx, r_in, r_out):
    X, Y = _grid(zeta, dx)
    RX = X - cx
    RY = Y - cy
    R = np.hypot(RX, RY)
    TH = np.arctan2(RY, RX)
    band = (R > r_in) & (R < r_out)
    n = int(band.sum())
    if n == 0:
        return 0.0, float("nan")
    A = np.sum(zeta[band] * np.exp(-1j * TH[band])) / n
    # e^{-i theta} gives arg(A) = -(angle of the positive lobe).
    return float(abs(A)), -math.degrees(math.atan2(A.imag, A.real))


def phase_profile(zeta, cx, cy, dx, radii_km, hw=25.0):
    out = []
    for rk in radii_km:
        _, ph = m1(zeta, cx, cy, dx, (rk - hw) * 1e3, (rk + hw) * 1e3)
        out.append(to_compass(ph) if ph == ph else float("nan"))
    return out


def steering(u2, v2, cx, cy, dx, disk=DISK_KM * 1e3, nb=80):
    """Symmetric-swirl-removed asymmetric flow averaged over the inner disk."""
    X, Y = _grid(u2, dx)
    RX = X - cx
    RY = Y - cy
    R = np.hypot(RX, RY)
    TH = np.arctan2(RY, RX)
    ut = -u2 * np.sin(TH) + v2 * np.cos(TH)
    rb = np.clip((R / disk * nb).astype(int), 0, nb - 1)
    utm = np.zeros(nb)
    for k in range(nb):
        m = rb == k
        if m.any():
            utm[k] = ut[m].mean()
    ut_sym = utm[rb]
    us = -ut_sym * np.sin(TH)
    vs = ut_sym * np.cos(TH)
    m = R < disk
    return float((u2 - us)[m].mean()), float((v2 - vs)[m].mean())


snaps = []
for f in NPZ:
    d = np.load(f)
    z = d["zeta"]
    dx = float(d["dx"])
    cx0 = float(d["cx"])
    cy0 = float(d["cy"])
    th = int(re.search(r"t(\d+)h", f).group(1))
    rcx, rcy = refine_center(z, cx0, cy0, dx)
    off = math.hypot(rcx - cx0, rcy - cy0) / 1e3
    rows = {b: m1(z, rcx, rcy, dx, b[0] * 1e3, b[1] * 1e3) for b in BANDS_KM}
    prof = phase_profile(z, rcx, rcy, dx, PROFILE_R_KM)
    su, sv = steering(d["u"], d["v"], rcx, rcy, dx)
    sdir = math.degrees(math.atan2(su, sv)) % 360.0
    sspd = math.hypot(su, sv)
    snaps.append(
        dict(
            t=th,
            cx=cx0,
            cy=cy0,
            rcx=rcx,
            rcy=rcy,
            off=off,
            rows=rows,
            prof=prof,
            sdir=sdir,
            sspd=sspd,
            su=su,
            sv=sv,
            z=z,
            dx=dx,
        )
    )

print("=" * 96)
print(f"GYRE RE-ANALYSIS v2  set={SNAP_TAG}  (re-centered; compass deg, 0=N CW)")
print("=" * 96)
hdr = f"  {'t(h)':>4}{'off_km':>8}"
for b in BANDS_KM:
    hdr += f"  {f'phi[{b[0]}-{b[1]}]':>12}"
hdr += f"{'steer_dir':>10}{'steer_ms':>9}"
print(hdr)
for s in snaps:
    line = f"  {s['t']:>4}{s['off']:>8.1f}"
    for b in BANDS_KM:
        _, ph = s["rows"][b]
        line += f"  {to_compass(ph):>12.0f}" if ph == ph else f"  {'nan':>12}"
    line += f"{s['sdir']:>10.0f}{s['sspd']:>9.2f}"
    print(line)

print("\n  SPIRAL  phi(r) [compass deg] - rows=radius(km), cols=t(h):")
print("    r\\t " + "".join(f"{s['t']:>7}" for s in snaps))
for i, rk in enumerate(PROFILE_R_KM):
    vals = "".join(
        f"{s['prof'][i]:>7.0f}" if s["prof"][i] == s["prof"][i] else f"{'--':>7}"
        for s in snaps
    )
    print(f"   {rk:>4} {vals}")
print("   (phase increasing inward = inner gyre wound poleward vs outer = sheared spiral)")

print("\n" + "=" * 96)
print("READ (clean window t12-48, excluding the t60 breakdown):")
clean = [s for s in snaps if s["t"] <= 48]
if len(clean) >= 2 and len(snaps) >= 2:
    dt = clean[-1]["t"] - clean[0]["t"]

    def seg_head(a, b):
        return math.degrees(math.atan2(b["cx"] - a["cx"], b["cy"] - a["cy"])) % 360

    drift_early = seg_head(snaps[0], snaps[1])
    drift_late = seg_head(clean[-2], clean[-1])
    drift_rate = crate(drift_early, drift_late, clean[-1]["t"] - snaps[1]["t"])
    print(
        f"  drift heading (center motion): {drift_early:.0f} -> {drift_late:.0f}  "
        f"({drift_rate:+.2f} deg/h)   [ground truth]"
    )
    for b in BANDS_KM:
        p0 = to_compass(clean[0]["rows"][b][1])
        p1 = to_compass(clean[-1]["rows"][b][1])
        r = crate(p0, p1, dt)
        tag = "  <-stationary" if abs(r) < 0.15 else (
            "  <-matches drift" if abs(r - drift_rate) < 0.20 else ""
        )
        print(f"  gyre phi[{b[0]}-{b[1]}km]: {p0:.0f} -> {p1:.0f}  ({r:+.2f} deg/h){tag}")
    s0 = clean[0]["sdir"]
    s1 = clean[-1]["sdir"]
    print(
        f"  STEERING flow (swirl removed): {s0:.0f} -> {s1:.0f}  "
        f"({crate(s0, s1, dt):+.2f} deg/h)   [should match drift if extraction is clean]"
    )
    moff = max(s["off"] for s in snaps)
    print(
        f"\n  max center offset: {moff:.1f} km - "
        f"{'raw inner band was CONTAMINATED; trust outer bands + steering' if moff > 8 else 'small'}"
    )
    print("  Key questions: (1) is the OUTER gyre stationary (precession was contamination)?")
    print("  (2) does STEERING dir now match the drift (rigorous gyre->drift link)?")
    print("  (3) does phi(r) increase inward (sheared spiral, inner gyre too poleward)?")

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(snaps)
    fig, axes = plt.subplots(1, n, figsize=(3.1 * n, 3.4))
    if n == 1:
        axes = [axes]
    HALF = 500e3
    for ax, s in zip(axes, snaps):
        z = s["z"]
        dx = s["dx"]
        rcx, rcy = s["rcx"], s["rcy"]
        X, Y = _grid(z, dx)
        R = np.hypot(X - rcx, Y - rcy)
        nb = 60
        rb = np.clip((R / HALF * nb).astype(int), 0, nb + 1)
        zb = np.zeros(nb + 2)
        for k in range(nb + 1):
            mm = rb == k
            if mm.any():
                zb[k] = z[mm].mean()
        zasym = z - zb[rb]
        xs = (np.arange(z.shape[0]) + 0.5) * dx
        ix = np.abs(xs - rcx) < HALF
        iy = np.abs(xs - rcy) < HALF
        sub = zasym[np.ix_(ix, iy)]
        vmax = np.nanmax(np.abs(sub)) or 1.0
        ax.imshow(
            sub.T,
            origin="lower",
            extent=[-HALF / 1e3, HALF / 1e3, -HALF / 1e3, HALF / 1e3],
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            aspect="equal",
        )
        ax.set_title(f"t={s['t']}h  off={s['off']:.0f}km", fontsize=9)
        ax.axhline(0, color="k", lw=0.3)
        ax.axvline(0, color="k", lw=0.3)
        sc = 100.0
        ax.annotate(
            "",
            xy=(s["su"] * sc, s["sv"] * sc),
            xytext=(0, 0),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6),
        )
        ax.set_xlabel("E-W (km)", fontsize=8)
    axes[0].set_ylabel("N-S (km)", fontsize=8)
    fig.suptitle("beta-gyre m=1 asymmetry (azimuthal mean removed), re-centered", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130)
    print(f"\n  wrote {OUT_PNG}")
except Exception as e:
    print(f"\n  (figure skipped: {e})")
