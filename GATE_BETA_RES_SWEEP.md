# gate-beta resolution sweep (paste-safe)

New mode `gate_beta_res`. The NU4 probe came back **diffusion-NOT-tunable**: only the baseline
3.0e11 ran clean (350°, Vmax_end 42); 1e11 was already contaminated (Vmax_end **82** = grid-scale
noise inflating the wind max, "heading 7°" is junk and points *east* anyway), and ≤3e10 went NaN.
NU4 is load-bearing for stability at dx=15.6 km — it removes real grid-scale noise, so it can't be
lowered to test the gyre. That points at **resolution**, not diffusion.

> ⚠ The NU4 mode's auto-READ line was **wrong** — it reported "350 → 7 rotates NW" but 7° is NNE
> (east of north) and 350° is NNW; the comparison was fooled by the 360/0 wrap and by a too-loose
> `Vmax<120` stability guard that let the contaminated 1e11 row through. This mode fixes both:
> signed **circular** heading deltas, and a contamination flag tied to the cap.

Method: hold NU4=3.0e11 and the in-band physics (R_env 500, taper-start 200, Vmax 64, cap 70, β,
u=v=0, 48 h) and refine **nx {320, 480, 640}** on a fixed 5000 km domain (dx 15.6 → 10.4 → 7.8 km),
watching the mature heading.

- **Heading rotates toward the NW band (350 → <335, i.e. signed delta NEGATIVE/westward) as dx
  shrinks, Vmax_end staying ~40–65** → the gyre was **under-resolved**; production needs finer dx
  (and the ~8° floor is a resolution artifact, not fundamental).
- **Heading ~flat across dx (|delta| ≲ 5°)** → the ~8° is **structural** to the β-gyre representation
  itself, not resolution-curable → the honest paper line is "characterized residual, ~8° poleward at
  the gyre scale," and the storm cross-track E (~+125 km) is a known, bounded systematic.

The 320 row re-runs the baseline as a determinism/consistency check — it should reproduce ~350°.

## Drop-in

```python
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
                                    r_env=500e3, nx=nx, dom=DOM, hours=48.0)
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
```

Dispatch:

```python
    elif arg in ("gate-beta-res", "gbeta-res", "res", "21"):
        gate_beta_res()
```

## Run

```
python run_translation_test.py gate-beta-res
```

Three 48 h runs. `DT=30` is fixed, so step count is the same; cost scales with the cell count, so
roughly **320 ≈ 12 min, 480 ≈ 28 min, 640 ≈ 50 min → ~90 min total**. The 320 row is a re-confirm
of the NU4 baseline (~350°); if you're impatient you can comment 320 out of `NX_SWEEP` and lean on
the existing 350° as the reference — the SUMMARY's Δ column just won't compute (set `base_hdg`
manually to 350 if you do).

## What the readout decides

- **`mature_hdg` dropping (350 → 340 → 330…), Δ-column going negative, `Vmax_end` ~40–65** = the
  β-gyre was under-resolved. The ~8° floor is a grid artifact, and the real fix is production
  resolution. Next: re-run Ivan full-physics at the winning dx and check whether the storm's
  eastward residual finally moves more than the taper's ~15 km.
- **`mature_hdg` flat (Δ within ±5°) across all three** = the ~8° is structural to how the model
  represents the β-gyre at all — not outer structure, not diffusion, not resolution. That's a clean,
  defensible paper statement (a *characterized* bias, not a mystery), and it closes the track-error
  investigation: Oracle's residual is a bounded ~8° poleward aim → ~+125 km cross-track E, same sign
  across all three storms.
- **A finer grid goes NON-FINITE** = nu4=3e11 is under the hyperdiff CFL there, so a blow-up means
  genuine under-resolution biting through the noise sink; scale `NU4 ∝ dx⁴` for that row and rerun
  before concluding.

Either branch is a real finding. The interesting one is "flat" — it would mean we've chased the aim
error through outer structure, diffusion, *and* resolution and found it irreducible at the gyre
scale, which is exactly the kind of honest, bounded error statement a BAMS reviewer respects more
than a tuned-away number.
