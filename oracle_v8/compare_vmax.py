#!/usr/bin/env python3
"""
compare_vmax.py — Vmax-dependence of the β-gyre spiral (swirl-shear vs β-Rossby).

Reads gyre_snap_v{VMAX}_t{T}h.npz sets (from `gate-beta-gyre 64/35/21`), groups by Vmax, and for
each snapshot computes:
  Vmax_est (max wind in the core), inner-lobe bearing [75-250km], outer-lobe bearing [300-500km],
  WIND-UP (inner - outer angular separation = spiral tightness), swirl-removed STEERING bearing+speed.
READ (mature t24-48 means): does wind-up / steering poleward-offset SCALE with Vmax?
  scales with Vmax  -> SWIRL-SHEAR (spiral wound by the vortex's own rotation; bias scales w/ intensity)
  ~Vmax-independent -> not swirl-driven (β-Rossby response; matches the Vmax-indep drift rate)

Run as a plain script:  python compare_vmax.py
"""
import glob, math, re
import numpy as np

INNER = (75, 250); OUTER = (300, 500); DISK_KM = 200.0; CORE_WIN_KM = 120.0


def _grid(a, dx):
    nx, ny = a.shape
    xs = (np.arange(nx) + 0.5) * dx
    return np.meshgrid(xs, xs, indexing="ij")


def to_compass(m):
    return (90.0 - m) % 360.0


def cdelta(a, b):
    return (((b - a + 180) % 360) - 180)


def refine_center(z, cx, cy, dx, win=CORE_WIN_KM * 1e3):
    X, Y = _grid(z, dx)
    R = np.hypot(X - cx, Y - cy)
    w = np.where((R < win) & (z > 0), z, 0.0)
    tot = w.sum()
    return (cx, cy) if tot <= 0 else (float((w * X).sum() / tot), float((w * Y).sum() / tot))


def m1_bearing(z, cx, cy, dx, r_in, r_out):
    X, Y = _grid(z, dx)
    RX = X - cx; RY = Y - cy
    R = np.hypot(RX, RY); TH = np.arctan2(RY, RX)
    b = (R > r_in) & (R < r_out)
    n = int(b.sum())
    if n == 0:
        return float("nan")
    A = np.sum(z[b] * np.exp(-1j * TH[b])) / n
    lobe = -math.degrees(math.atan2(A.imag, A.real))    # +lobe angle, sign-corrected
    return to_compass(lobe)


def steering(u2, v2, cx, cy, dx, disk=DISK_KM * 1e3, nb=80):
    X, Y = _grid(u2, dx)
    RX = X - cx; RY = Y - cy
    R = np.hypot(RX, RY); TH = np.arctan2(RY, RX)
    ut = -u2 * np.sin(TH) + v2 * np.cos(TH)
    rb = np.clip((R / disk * nb).astype(int), 0, nb - 1)
    utm = np.zeros(nb)
    for k in range(nb):
        m = rb == k
        if m.any():
            utm[k] = ut[m].mean()
    uts = utm[rb]
    us = -uts * np.sin(TH); vs = uts * np.cos(TH)
    m = R < disk
    su = float((u2 - us)[m].mean()); sv = float((v2 - vs)[m].mean())
    return math.degrees(math.atan2(su, sv)) % 360, math.hypot(su, sv)


def vmax_est(u2, v2, cx, cy, dx):
    X, Y = _grid(u2, dx)
    core = np.hypot(X - cx, Y - cy) < 300e3
    return float(np.sqrt(u2 ** 2 + v2 ** 2)[core].max())


files = glob.glob("gyre_snap_v*_t*h.npz")
if not files:
    raise SystemExit("no gyre_snap_v*_t*h.npz — run `gate-beta-gyre 64/35/21` (parameterized mode)")
groups = {}
for f in files:
    mv = re.search(r"_v(\d+)_t(\d+)h", f)
    groups.setdefault(int(mv.group(1)), []).append((int(mv.group(2)), f))
for vm in groups:
    groups[vm].sort()

print("=" * 96)
print("Vmax-DEPENDENCE OF THE β-GYRE SPIRAL  (swirl-shear vs β-Rossby; COMPASS deg)")
print("=" * 96)
summary = {}
for vm in sorted(groups, reverse=True):
    print(f"\n  --- init Vmax {vm} ---")
    print(f"  {'t(h)':>4}{'Vmax_est':>9}{'inner':>8}{'outer':>8}{'windup':>8}{'steer':>7}{'ms':>6}")
    mature = []
    for th, f in groups[vm]:
        d = np.load(f); z = d["zeta"]; dx = float(d["dx"])
        rcx, rcy = refine_center(z, float(d["cx"]), float(d["cy"]), dx)
        ib = m1_bearing(z, rcx, rcy, dx, INNER[0] * 1e3, INNER[1] * 1e3)
        ob = m1_bearing(z, rcx, rcy, dx, OUTER[0] * 1e3, OUTER[1] * 1e3)
        wind = cdelta(ob, ib)
        sb, sm = steering(d["u"], d["v"], rcx, rcy, dx)
        ve = vmax_est(d["u"], d["v"], rcx, rcy, dx)
        print(f"  {th:>4}{ve:>9.1f}{ib:>8.0f}{ob:>8.0f}{wind:>8.0f}{sb:>7.0f}{sm:>6.2f}")
        if 24 <= th <= 48:
            mature.append((ve, wind, sb))
    if mature:
        ve_m = float(np.mean([m[0] for m in mature]))
        wind_m = float(np.mean([m[1] for m in mature]))
        sb_m = math.degrees(math.atan2(
            np.mean([math.sin(math.radians(m[2])) for m in mature]),
            np.mean([math.cos(math.radians(m[2])) for m in mature]))) % 360
        summary[vm] = (ve_m, wind_m, sb_m)

print("\n" + "=" * 96); print("READ (mature t24-48 means):")
print(f"  {'initVmax':>9}{'Vmax_est':>10}{'|windup|':>10}{'steer_brg':>11}{'poleward_of_NW':>16}")
for vm in sorted(summary, reverse=True):
    ve, wind, sb = summary[vm]
    print(f"  {vm:>9}{ve:>10.1f}{abs(wind):>10.0f}{sb:>11.0f}{cdelta(315.0, sb):>+16.0f}")
if len(summary) >= 2:
    vms = sorted(summary, reverse=True)
    w_hi, w_lo = abs(summary[vms[0]][1]), abs(summary[vms[-1]][1])
    p_hi, p_lo = cdelta(315.0, summary[vms[0]][2]), cdelta(315.0, summary[vms[-1]][2])
    print()
    if (w_hi - w_lo) > 25 or (p_hi - p_lo) > 15:
        print(f"  → wind-up &/or poleward-offset GROW with Vmax (Vmax {vms[-1]}→{vms[0]}: "
              f"|windup| {w_lo:.0f}→{w_hi:.0f}°, poleward {p_lo:+.0f}→{p_hi:+.0f}°) → SWIRL-SHEAR: "
              "the spiral is wound by the vortex's own rotation, so the poleward bias scales with "
              "intensity. Reconciles time-evol (drift RATE ~Vmax-indep, but the structural wind-up "
              "and steering-offset are not). Paper: stronger storms → larger poleward track bias.")
    else:
        print(f"  → wind-up & poleward-offset ~Vmax-INDEPENDENT (|windup| {w_lo:.0f}→{w_hi:.0f}°, "
              f"poleward {p_lo:+.0f}→{p_hi:+.0f}°) → NOT swirl-driven → β-Rossby gyre response, "
              "consistent with the Vmax-indep drift rate. The bias is intrinsic to the β-gyre, not "
              "the vortex strength.")
