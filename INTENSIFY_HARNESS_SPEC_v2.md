# Intensifying-Vortex Harness — spec v2 (V8.7-pre)

**Supersedes** the `ramp_theta` / I3-thermo branch of v1. The run-script files
resolved the open question, and the answer rules out the thermodynamic path.

## What the files settled

Ivan re-intensifies through an **emergent barotropic** mechanism — angular-
momentum spin-up from the advective/anelastic dynamics, the "inviscid runaway"
(max|u|→177 m/s uncapped) that `SurfaceDrag` + `IntensityCap` *contain*
(run_ivan.py L318–328). θ′ is **zeroed at init** (L267, barotropic) and actively
suppressed by Newtonian cooling, which the code explicitly notes is irrelevant to
the intensification. There is no buoyancy and no heating term in the config.

Consequences:

1. **`ramp_theta=True` is anti-faithful — drop it.** Driving θ′ up is the
   opposite of what the storm does. Delete the thermo branch from I3.
2. **The imposed kinematic ramp (v1 I0/I1) is a Galilean control, not a storm
   replica.** It forces the wind field up by fiat; the storm lets the dynamics
   spin it up. Keep it — it bounds the *magnitude* effect — but read it as a
   control (see below), not a verdict.
3. **Ivan's structure == the harness's structure.** `RMAX_RUN_M = 75 km`,
   `R_ENV_M = 500 km`, `WIND_TAPER=True`, `dx = 15.6 km`, cap 70, drag — all
   identical to the ladder defaults. So the difference between ladder-L0
   (decays 64→49 in 12 h) and Ivan (spins 44→84 over 52 h) is **not** Rmax/R_env.
   It is some combination of: init Vmax, β-plane, time-varying steering, B, and
   **duration**.
4. **The original ladder was too short to test intensity at all.** Barotropic
   spin-up needs >24 h to develop; the ladder ran 12 h. "Re-intensification is
   the last suspect" is right, but the ladder couldn't have seen it either way.

## How to read the I0/I1 run currently going

`intensify` (v1) forces the perturbation winds up 44→84 on an f-plane and reads
eff_y. Interpret:

- **eff_y flat ≈ 1.00 through the ramp** → advection is faithful to forced
  *wind magnitude*. This exonerates the crudest (b): "stronger winds → numerical
  translation error." It does **not** exonerate the structural version
  (gradient sharpening / contraction), because a uniform rescale keeps the
  profile shape fixed. Bound, not verdict → run the emergent ladder.
- **eff_y rises with forced Vmax(t)** → there is a wind-magnitude-dependent
  frame-breaking term, independent of structure. Strong (b) signal; the emergent
  ladder then says whether the storm actually drives it that hard.

Either way the next experiment is the same.

---

## The intensification ladder — the storm-matched test

Same idea as the original ingredient ladder, but (i) start at Ivan's init Vmax,
(ii) run **52 h** so the barotropic spin-up can emerge, (iii) read eff_y(t) and
realized **Vmax(t) together**, (iv) add β then steering one rung at a time. The
rung where the eff excess appears — *and whether it appears with the
intensification or independently* — names the mechanism.

```
J0  f-plane, constant 5 m/s N bg, Ivan init   → does it self-intensify on an
                                                 f-plane?  eff_y(t) vs Vmax(t):
                                                 intensity-coupled NUMERICS channel
J1  β-plane, zero bg, Ivan init               → β self-propagation vs Vmax(t);
                                                 compare to gate-beta static β(Vmax):
                                                 intensity-coupled β-DRIFT channel
J2  β-plane + steering ramp, Ivan init        → full Ivan-like config minus the
                                                 real ERA5/track; does it reproduce
                                                 BOTH 44→84 AND the eff excess?
```

The decisive contrasts:

- **J0 intensifies and eff_y tracks Vmax(t)** → intensity-coupled numerics is
  real and structural (the thing I1's uniform rescale couldn't fully test).
- **J0 does NOT intensify (settles ~49 like the 12 h ladder, just slower)** →
  the f-plane structure+drag+cap equilibrium is ~49; Ivan's extra spin-up to 84
  **requires β/steering** → look at J1/J2. This would mean **the intensification
  itself is β/steering-coupled, and the over-translation is plausibly the same
  coupling** — the deepest possible finding, and it redirects the hunt from
  "advection numerics" to "β-vortex / steering-vortex interaction."
- **J1 β-drift grows faster than gate-beta's static β(Vmax)** → intensification-
  transient β enhancement (intensity-coupled (a), beyond superposition).
- **J2 reproduces Ivan's +2.9 m/s and 44→84** → the idealized harness captures
  Ivan; ablate from there (kill β, kill steering, freeze intensity) to isolate.
  **J2 does NOT reproduce it** → something storm-specific (real ERA5 ramp shape,
  real track curvature) is required — a different branch.

### Plumbing — `run_translation` (drop-in, additive)

Add `Rmax` and `B` overrides so a rung can match Ivan's vortex exactly:

```python
def run_translation(Vmax, u_env=0.0, v_env=5.0, epsilon=0.5,
                    drag_on=True, nx=128, hours=10.0, v_cap=None, dom=None,
                    r_env=None, beta=False, wind_taper=False,
                    taper_start_frac=0.5, keep_theta=False, steer_ramp=None,
                    Rmax=None, B=None,                       # <-- NEW (V8.7)
                    subcell=True, verbose=False):
    ...
    f0   = KATRINA["f"]
    renv = r_env if r_env is not None else R_ENV_M
    rmax = Rmax if Rmax is not None else RMAX_M             # <-- NEW
    bb   = B if B is not None else KATRINA["B"]             # <-- NEW
    init = HollandVortexInit(Vmax=Vmax, Rmax=rmax, B=bb, f=f0,
                             R_env=renv, u_env=u_env, v_env=v_env,
                             wind_taper=wind_taper,
                             taper_start_frac=taper_start_frac)
```

(No `vmax_ramp` / `ramp_theta` needed for the J-ladder — intensification is
emergent. Keep the v1 `vmax_ramp` plumbing only if you still want I1 as a
control; it's orthogonal.)

### New mode — `intensify_ladder()`

```python
def intensify_ladder():
    """EMERGENT INTENSIFICATION LADDER (V8.7) — the storm-matched test.

    Ivan's structure == this harness's structure (Rmax 75km, R_env 500km, taper,
    dx 15.6km, cap 70, drag).  The barotropic spin-up that took Ivan 44→84 needs
    >24h to emerge — the 12h ladder was too short to see it.  Run 52h at Ivan's
    init, add β then steering, and read eff_y(t) against EMERGENT Vmax(t).  The
    rung where the excess appears, and whether it rides the intensification or
    not, names the mechanism.
    Usage:  python run_translation_test.py intensify-ladder
    """
    from oracle_v8.storm_data import IVAN
    t0 = time.time()
    V0  = IVAN["Vmax_ms"]          # Ivan HURDAT2 t=0 intensity
    Bi  = IVAN["B"]
    print("=" * 78)
    print(f"INTENSIFICATION LADDER  (Ivan grid 320²/5000km, taper ON, "
          f"r_env=500km, init Vmax={V0:.0f}, B={Bi:.2f}, 52h)")
    print("  let the barotropic spin-up EMERGE; read eff_y(t) vs Vmax(t)")
    print("=" * 78)
    common = dict(nx=320, dom=5_000_000.0, v_cap=70.0, wind_taper=True,
                  r_env=500e3, hours=52.0, Rmax=75_000.0, B=Bi)
    ramp = (-1.1, 3.9, +0.6, 6.6)        # Ivan DLM endpoints (relaxation spans run)

    def _trace(d, v_bg, label):
        tt, xs, ys, vm = d["track"]
        print(f"\n  {label} — Vmax(t) and translation:")
        print(f"  {'t(h)':>5}  {'Vmax':>6}  {'y(km)':>8}  {'x(km)':>8}  "
              f"{'vy(m/s)':>8}  {'eff_y_cum':>9}")
        step = max(1, len(tt) // 13)          # ~13 rows
        for k in range(step, len(tt), step):
            dt_s = tt[k] - tt[k - 1]
            vy = (ys[k] - ys[k - 1]) / dt_s
            eff = ((ys[k] - ys[0]) / (v_bg * tt[k]) if v_bg > 0
                   else float("nan"))
            ef = f"{eff:+9.3f}" if v_bg > 0 else "      -- "
            print(f"  {tt[k]/3600:5.1f}  {vm[k]:6.1f}  {ys[k]/1e3:8.1f}  "
                  f"{xs[k]/1e3:8.1f}  {vy:+8.2f}  {ef}")
        print(f"  Vmax: init {vm[0]:.0f} → end {vm[-1]:.0f} "
              f"(peak {max(vm):.0f})")

    print("\n  J0  f-plane, constant 5 m/s N bg  (does it self-intensify?):")
    j0 = run_translation(V0, u_env=0.0, v_env=5.0, **common)
    print(_grow("J0 fplane-52h", j0)); _trace(j0, 5.0, "J0")

    print("\n  J1  β-plane, zero bg  (β self-propagation vs intensity):")
    j1 = run_translation(V0, u_env=0.0, v_env=0.0, beta=True, **common)
    print(f"  J1 end β-drift ({j1['drift_x']:+.2f},{j1['drift_y']:+.2f}) m/s")
    _trace(j1, 0.0, "J1")

    print("\n  J2  β-plane + steering ramp  (full Ivan-like config):")
    j2 = run_translation(V0, u_env=ramp[0], v_env=ramp[1], beta=True,
                         steer_ramp=ramp, **common)
    print(_grow("J2 beta+ramp-52h", j2)); _trace(j2, None or 0.0, "J2")

    print("\n" + "=" * 78)
    print("READ:")
    print("  J0 intensifies AND eff_y rises with Vmax(t)  → intensity-coupled")
    print("     NUMERICS, structural (the part I1's uniform rescale can't test).")
    print("  J0 does NOT intensify (settles ~49)          → Ivan's spin-up needs")
    print("     β/steering → intensification is β/steering-coupled, and the over-")
    print("     translation is plausibly the SAME coupling (J1/J2 confirm).")
    print("  J1 β-drift grows faster than gate-beta static β(Vmax)  → transient")
    print("     β enhancement during active intensification.")
    print("  J2 reproduces Ivan's +2.9 m/s / 44→84        → idealized harness")
    print("     captures Ivan; ablate to isolate.  Doesn't  → real ERA5/track")
    print("     specifics required, separate branch.")
    print(f"\nWall time: {time.time()-t0:.0f}s")
```

Dispatch:

```python
    elif arg in ("intensify-ladder", "iladder", "16"):
        intensify_ladder()
```

### Cost

52 h = 6240 steps. At the validated ~0.16 s/step that's ~17 min/rung, ~50 min
for J0–J2. Pair with `gate-beta` (the static β(Vmax) reference J1 is compared
against) and the whole intensity question closes in one GPU session.

### Sequencing with the run you have going

1. Let `intensify` (I0/I1) finish — it's the magnitude bound.
2. Run `gate-beta` — gives the static β(Vmax) curve J1 needs.
3. Run `intensify-ladder` — the verdict.

If J0 fails to intensify on the f-plane, that's not a null result — it's the
signal to stop hunting advection numerics and start looking at the β/steering-
vortex interaction as the shared root of *both* the intensification and the
over-translation.
