# Formulation probe — two independent arms (paste-safe)

Peer review correctly flagged that "structural" is premature: we tested outer structure, diffusion
*magnitude* (NU4), and resolution, but never the **formulation**. This probe tests the two
formulation knobs the reviewer named, each in isolation so a null in one can't mask an effect in the
other. Both run at the in-band config (β, u=v=0, Vmax 64, R_env 500, taper-start 200, cap 70, nx 320,
5000 km, 48 h) and read the mature (30–48 h) drift with the corrected circular-delta + cap-tied
contamination logic from `gate-beta-res`.

**Arm A — diffusion FORM** (`gate-beta-diffform`): ∇⁴ baseline vs ∇² (Laplacian) at three nu_H.
**Arm B — divergence damper** (`gate-beta-divdamp`): sweep Helmholtz `epsilon` 0.5 → 0 (off).

Read both against the **westward component** (the honest metric), not magnitude. Reference ∇⁴/eps=0.5
aim ≈ 350° (≈ NNW, west-comp ~0.4 m/s); canonical NW drift would be ~315°, west-comp ~1.4 m/s.

### Orientation on what each outcome means
- **∇² (Arm A) is MORE dissipative at the gyre than ∇⁴**, not less. So Arm A tests the
  more-dissipation direction. If ∇² rotates the aim (either way) → the aim is dissipation-sensitive,
  and the *fix* direction (less dissipation than ∇⁴) is a SHARPER filter (∇⁶/spectral) — I'll build
  that only if this arm shows sensitivity. If ∇² is flat even while it eats the gyre amplitude
  (watch |drift|), the aim is decoupled from gyre dissipation.
- **If BOTH arms are flat**, diffusion and divergence are exonerated, and the seat is the Coriolis
  term itself or the projection's geostrophic adjustment — which would make the original "structural"
  read correct *but now properly earned*, and points the next (harder) probe at a true-β / non-periodic
  projection test.

---

## Shared helpers (paste ONCE, above both modes)

```python
def _circ_delta(h, ref):
    """signed h-ref in (-180,180]; NEGATIVE = counter-clockwise = westward (toward NW)."""
    return ((h - ref + 180.0) % 360.0) - 180.0

def _mature_drift(track, ta_h=30.0, tb_h=48.0):
    """Return (speed m/s, heading deg [0=N,CW], westward-component m/s) over [ta,tb]."""
    import math
    tt, xs, ys, vm = track
    ia = min(range(len(tt)), key=lambda k: abs(tt[k] - ta_h * 3600.0))
    ib = min(range(len(tt)), key=lambda k: abs(tt[k] - tb_h * 3600.0))
    dts = tt[ib] - tt[ia]
    spd = math.hypot((xs[ib]-xs[ia])/dts, (ys[ib]-ys[ia])/dts)
    hdg = math.degrees(math.atan2((xs[ib]-xs[ia])/dts, (ys[ib]-ys[ia])/dts)) % 360
    west = -spd * math.sin(math.radians(hdg))   # +ve = westward
    return spd, hdg, west
```

---

## ARM A — diffusion form. THREE harness edits, then the mode.

### Edit 1 — import (add next to HyperDiffusionComponent)
```
    HyperDiffusionComponent,
```
becomes
```
    HorizontalDiffusionComponent,
    HyperDiffusionComponent,
```

### Edit 2 — run_translation signature (add two kwargs)
```
                    Rmax=None, B=None,
```
becomes
```
                    Rmax=None, B=None, diff_form="hyper", nu_H=2.0e5,
```

### Edit 3 — make the diffusion component switchable
```
        horiz_diffusion=HyperDiffusionComponent(nu4=NU4, Lx=Lx, Ly=Ly,
                                                nx=nx, ny=ny),
```
becomes
```
        horiz_diffusion=(HorizontalDiffusionComponent(nu_H=nu_H, Lx=Lx, Ly=Ly,
                                                      nx=nx, ny=ny)
                         if diff_form == "laplacian" else
                         HyperDiffusionComponent(nu4=NU4, Lx=Lx, Ly=Ly,
                                                 nx=nx, ny=ny)),
```

### Drop-in mode
```python
def gate_beta_diffform():
    """GATE-BETA DIFFUSION-FORM PROBE (V8.7, Arm A). Does the aim depend on the
    diffusion OPERATOR/strength?  NU4 couldn't be LOWERED (load-bearing); this tests
    operator form + the MORE-dissipation direction: ∇⁴ ref vs ∇² at 3 nu_H.  ∇² damps
    the resolved gyre far more than ∇⁴, so: aim moves (heading Δ, west-comp) → aim is
    dissipation-sensitive → try a SHARPER filter (∇⁶/spectral) for the less-dissipation
    fix.  aim FLAT while |drift| drops → decoupled from gyre amplitude → not diffusion.
    Requires the 3 harness edits above.  Usage: python run_translation_test.py gate-beta-diffform
    """
    import math
    t0 = time.time()
    TH_HDG = (290.0, 335.0); CAP = 70.0
    RUNS = [
        ("nabla4  NU4=3.0e11 (ref)",          dict(diff_form="hyper")),
        ("nabla2  nu_H=1.2e4 (grid-matched)", dict(diff_form="laplacian", nu_H=1.2e4)),
        ("nabla2  nu_H=5.0e4",                dict(diff_form="laplacian", nu_H=5.0e4)),
        ("nabla2  nu_H=2.0e5 (recommended)",  dict(diff_form="laplacian", nu_H=2.0e5)),
    ]
    print("=" * 78)
    print("GATE-BETA DIFFUSION-FORM PROBE (Arm A)  (β, u=v=0, Vmax=64, R_env=500, "
          "taper-start 200, cap 70, 320/5000km, 48h)")
    print(f"  ref ∇⁴ ~350°; NW band {TH_HDG[0]:.0f}-{TH_HDG[1]:.0f}. aim rotates WEST "
          "(heading↓, west-comp↑) → form/dissipation matters; flat → it doesn't")
    print("=" * 78)
    rows = []
    for label, kw in RUNS:
        try:
            d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=CAP, beta=True,
                                wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                                nx=320, dom=5_000_000.0, hours=48.0, **kw)
        except Exception as e:
            rows.append((label, float("nan"), float("nan"), float("nan"),
                         float("nan"), "CRASH"))
            print(f"\n  {label}: ⚠ raised ({type(e).__name__}) — unstable.")
            continue
        track = d["track"]; tt, xs, ys, vmt = track; ve = d["vmax_end"]
        finite = all(math.isfinite(v) for v in (xs[-1], ys[-1], ve))
        if finite:
            spd, hdg, west = _mature_drift(track)
        else:
            spd = hdg = west = float("nan")
        contam = finite and ve > CAP * 1.10
        verdict = ("CONTAM" if contam else ("OK" if finite else "UNSTABLE"))
        rows.append((label, hdg, spd, west, ve, verdict))
        print(f"\n  {label}:")
        if finite:
            print(f"    MATURE(30-48) |{spd:.2f}| hdg {hdg:.0f} ({_compass(hdg)})  "
                  f"WEST-comp {west:+.2f} m/s  Vmax_end={ve:.1f}  "
                  f"{'⚠CONTAM' if contam else 'clean'}")
        else:
            print("    ⚠ NON-FINITE.")
    print("\n" + "=" * 78); print("SUMMARY:")
    print(f"  {'form / coeff':<34}{'hdg':>5}{'Δvs∇⁴':>7}{'|drift|':>8}{'west':>7}"
          f"{'Vmax':>7}  verdict")
    ref_hdg = next((r[1] for r in rows if "ref" in r[0] and r[1] == r[1]), float("nan"))
    for label, hdg, spd, west, ve, vd in rows:
        hs = f"{hdg:.0f}" if hdg == hdg else "nan"
        dl = (f"{_circ_delta(hdg, ref_hdg):+.0f}"
              if (hdg == hdg and ref_hdg == ref_hdg) else "  —")
        ss = f"{spd:.2f}" if spd == spd else "nan"
        ws = f"{west:+.2f}" if west == west else "nan"
        vs = f"{ve:.0f}" if ve == ve else "nan"
        print(f"  {label:<34}{hs:>5}{dl:>7}{ss:>8}{ws:>7}{vs:>7}  {vd}")
    print("\nREAD:")
    lap = [r for r in rows if "nabla2" in r[0] and r[5] == "OK" and r[1] == r[1]]
    if lap and ref_hdg == ref_hdg:
        most = min(lap, key=lambda r: _circ_delta(r[1], ref_hdg))
        dl = _circ_delta(most[1], ref_hdg)
        if dl <= -6:
            print(f"  ∇² rotates aim WEST ({ref_hdg:.0f}→{most[1]:.0f}, {dl:+.0f}°) → aim IS "
                  "operator-sensitive. Surprising (∇² damps the gyre MORE) but actionable: "
                  "build a SHARPER filter (∇⁶/spectral) to push the less-dissipation fix.")
        elif dl >= 6:
            print(f"  ∇² rotates aim EAST/poleward ({ref_hdg:.0f}→{most[1]:.0f}, {dl:+.0f}°) → "
                  "more gyre dissipation worsens the aim → westward ventilation lives in the "
                  "resolved gyre; ∇⁴ is the less-damaging choice and the fix direction is LESS "
                  "dissipation → build a SHARPER filter (∇⁶/spectral) next.")
        else:
            print(f"  aim FLAT (Δ {dl:+.0f}° ≤ 6°) even as ∇² damps the gyre (check |drift|/west) "
                  "→ aim DECOUPLED from gyre dissipation → not a diffusion artifact → seat is the "
                  "Coriolis/projection path (see Arm B, then true-β/projection probe).")
    else:
        print("  ∇² runs unstable/contaminated — inspect Vmax; grid-matched ∇² may need a "
              "higher nu_H floor to stay clean before the comparison is meaningful.")
    print(f"\nWall time: {time.time()-t0:.0f}s")
```

Dispatch:
```python
    elif arg in ("gate-beta-diffform", "gbeta-diffform", "diffform", "22"):
        gate_beta_diffform()
```

---

## ARM B — divergence damper. NO harness edit; epsilon is already a parameter.

```python
def gate_beta_divdamp():
    """GATE-BETA DIVERGENCE-DAMPER PROBE (V8.7, Arm B). Does the aim depend on the
    Helmholtz divergence damper?  Reviewer hypothesis: it suppresses the divergent
    flow that communicates the β-effect (geostrophic adjustment / Rossby radiation) →
    too-poleward aim.  Sweep epsilon 0.5 (baseline) → 0 (off; harness drops the
    component).  aim rotates WEST as eps↓ → damper WAS suppressing the β-communicating
    divergent flow (reduce it / gentler form = fix).  flat → not the driver.  Unstable
    at low eps → damper load-bearing → next test is measuring the vorticity the discrete
    Helmholtz solve injects (reviewer's 'not exact on a grid').  NO harness edit needed.
    Usage:  python run_translation_test.py gate-beta-divdamp
    """
    import math
    t0 = time.time()
    TH_HDG = (290.0, 335.0); CAP = 70.0
    EPS_SWEEP = (0.5, 0.25, 0.1, 0.0)   # baseline → off
    print("=" * 78)
    print("GATE-BETA DIVERGENCE-DAMPER PROBE (Arm B)  (β, u=v=0, Vmax=64, R_env=500, "
          "taper-start 200, cap 70, 320/5000km, 48h)")
    print(f"  baseline eps=0.5 ~350°; NW band {TH_HDG[0]:.0f}-{TH_HDG[1]:.0f}. aim rotates "
          "WEST as eps↓ → divergent flow carries the β-effect; flat → it doesn't")
    print("  eps=0 drops the damper entirely (harness: if epsilon > 0)")
    print("=" * 78)
    rows = []
    for eps in EPS_SWEEP:
        try:
            d = run_translation(64, u_env=0.0, v_env=0.0, v_cap=CAP, beta=True,
                                wind_taper=True, taper_start_frac=0.40, r_env=500e3,
                                nx=320, dom=5_000_000.0, hours=48.0, epsilon=eps)
        except Exception as e:
            rows.append((eps, float("nan"), float("nan"), float("nan"),
                         float("nan"), "CRASH"))
            print(f"\n  eps={eps:.2f}: ⚠ raised ({type(e).__name__}) — unstable.")
            continue
        track = d["track"]; tt, xs, ys, vmt = track; ve = d["vmax_end"]
        finite = all(math.isfinite(v) for v in (xs[-1], ys[-1], ve))
        if finite:
            spd, hdg, west = _mature_drift(track)
        else:
            spd = hdg = west = float("nan")
        contam = finite and ve > CAP * 1.10
        verdict = ("CONTAM" if contam else ("OK" if finite else "UNSTABLE"))
        rows.append((eps, hdg, spd, west, ve, verdict))
        print(f"\n  eps={eps:.2f} {'(damper OFF)' if eps == 0 else ''}:")
        if finite:
            print(f"    MATURE(30-48) |{spd:.2f}| hdg {hdg:.0f} ({_compass(hdg)})  "
                  f"WEST-comp {west:+.2f} m/s  Vmax_end={ve:.1f}  "
                  f"{'⚠CONTAM' if contam else 'clean'}")
        else:
            print("    ⚠ NON-FINITE.")
    print("\n" + "=" * 78); print("SUMMARY:")
    print(f"  {'eps':>5}{'hdg':>6}{'Δvs0.5':>8}{'|drift|':>8}{'west':>7}{'Vmax':>7}  verdict")
    ref_hdg = next((r[1] for r in rows if r[0] == 0.5 and r[1] == r[1]), float("nan"))
    for eps, hdg, spd, west, ve, vd in rows:
        hs = f"{hdg:.0f}" if hdg == hdg else "nan"
        dl = (f"{_circ_delta(hdg, ref_hdg):+.0f}"
              if (hdg == hdg and ref_hdg == ref_hdg) else "  —")
        ss = f"{spd:.2f}" if spd == spd else "nan"
        ws = f"{west:+.2f}" if west == west else "nan"
        vs = f"{ve:.0f}" if ve == ve else "nan"
        print(f"  {eps:>5.2f}{hs:>6}{dl:>8}{ss:>8}{ws:>7}{vs:>7}  {vd}")
    print("\nREAD:")
    clean = [r for r in rows if r[5] == "OK" and r[1] == r[1]]
    if clean and ref_hdg == ref_hdg:
        most = min(clean, key=lambda r: _circ_delta(r[1], ref_hdg))
        dl = _circ_delta(most[1], ref_hdg)
        if dl <= -6:
            print(f"  aim rotates WEST as eps↓ ({ref_hdg:.0f}→{most[1]:.0f} at eps={most[0]:.2f}, "
                  f"{dl:+.0f}°) → the damper WAS suppressing the β-communicating divergent flow. "
                  "Reduce eps / gentler form → re-confirm at storm scale.")
        else:
            print(f"  aim ~FLAT across eps ({ref_hdg:.0f}→{most[1]:.0f}, {dl:+.0f}°) → divergence "
                  "damping is NOT the aim driver.")
        bad = [r[0] for r in rows if r[5] in ("UNSTABLE", "CRASH", "CONTAM")]
        if bad:
            print(f"  NOTE low-eps rows unstable/contaminated {bad} → damper partly load-bearing; "
                  "if it can't be lowered cleanly, next test is measuring the vorticity the discrete "
                  "Helmholtz ψ-solve injects (∇×ψ-correction ≠ 0 on a discrete grid).")
    else:
        print("  No clean rows — damper load-bearing for stability; measure injected vorticity "
              "instead of reducing eps.")
    print(f"\nWall time: {time.time()-t0:.0f}s")
```

Dispatch:
```python
    elif arg in ("gate-beta-divdamp", "gbeta-divdamp", "divdamp", "23"):
        gate_beta_divdamp()
```

---

## Run (independent)

```
python run_translation_test.py gate-beta-diffform   # Arm A, 4 runs, ~48 min
python run_translation_test.py gate-beta-divdamp     # Arm B, 4 runs, ~48 min
```

Run them in either order, separately. Arm A needs the three harness edits first; Arm B runs as-is.

## Decision tree after both

| Arm A (∇²) | Arm B (eps) | Conclusion / next |
|---|---|---|
| aim moves | — | dissipation-sensitive → build ∇⁶/spectral sharp filter to test the less-dissipation FIX |
| — | aim moves | divergence-flow carries β → reduce/replace damper, re-confirm on a storm |
| flat | flat | diffusion + divergence exonerated → seat is Coriolis term or projection adjustment → true-β / non-periodic-projection probe (the "structural" read would then be earned) |
| unstable rows | unstable rows | knob is load-bearing (like NU4) → can't reduce → measure what it injects, don't infer from absence |

Whatever moves the aim — or the clean double-null — is a real result. A double-null is the one that
vindicates last night's "structural" conclusion, but only because we'll have actually tested the
formulation instead of assuming it.
