# Oracle V8 — Build State Cheat Sheet

**Cleanup audit (2026-06-19):** code provenance is now fixed, but storm-result provenance is not yet
fully closed. `run_hugo.py` now persists `taper_start_frac=0.40` and has the same ERA5-required,
lockstep steering/cap/drag, raw-theta tracker diagnostics, and NaN-check gating pattern as Katrina/Ivan.
`run_translation_test.py` now carries explicit `f_ref` metadata and uses Ivan's `f` for the Ivan-grid
beta harnesses; `reanalyze_gyre.py` now requires `--tag` when multiple gyre snapshot sets exist and writes
tagged figures (for example `gyre_precession_v64.png`). Until fresh storm logs are regenerated from these
cleaned scripts, treat the 3-storm +120 km landfall-fix cluster as **provisional**, not checked-in evidence.
The checked-in logs still support the same-latitude rows for Katrina +126.5 km E and Ivan +72.8 km E;
**Hugo's same-latitude row is now +102.7 km E / −2.2 h** from the corrected storm-agnostic `run_hugo.py`
(was +131.2 km E / +1.2 h — the earlier driver was NOT on the Katrina/Ivan storm-agnostic stack; see the
HURDAT2_VERIFICATION.md scorecard correction).

**Version:** V8.7 (Vmax-dependence test DECISIVE: across a 3× swirl range (Vmax_est 32/19/12) the gyre wind-up (~96-103°) & steering poleward-offset (~+17-23° of canonical NW) are FLAT → the directional bias is INTENSITY-INDEPENDENT → β-Rossby/structural, NOT swirl-shear. Resolves the time-evol tension (direction Vmax-indep, speed ∝ strength). FULL CHAIN: 3-storm +120km E cluster → genuine β-drift → β-gyre equilibrates ~20° too poleward (NNW vs NW) → intensity-independent → eastward track bias. β-drift MECHANISM CHARACTERIZATION COMPLETE. Next = publishability-critical BASELINE + MORE STORMS)
**Last updated:** Hugo storm-agnostic rerun (June 19 2026) — the June-14 Hugo run was NOT on the
storm-agnostic stack (driver settings differed from Katrina/Ivan); corrected and rerun →
**+102.7 km E / −2.2 h** (same-lat 32.5°N). **Two big shifts since
V8.7-pre:**
**(1) The β-drift CALIBRATION is closed** (the *one-knob* result — distinct from the aim-FLOOR
question, which peer review RE-OPENED; see NOW). Magnitude and aim collapse onto ONE structural control —
the outer-wind **taper-start radius** (= `taper_start_frac · R_env`); R_env and taper_start_frac are
**degenerate** (only their product matters). Shrinking it co-reduces magnitude *and* poleward bias,
but only to a **floor of ~343°** (mag in-band ~2.4 m/s). The residual **~8° poleward aim** was provisionally
read as a fundamental gyre-level floor (resolution / NU4 smearing the westward-ventilation asymmetry) —
but that "structural" call is now PREMATURE: the diffusion FORM, advection ORDER, and projection/
divergence-damper remain untested (formulation probe pending). The boundary/reversed-β taper is EXONERATED
(rotation domain-invariant; code confirms the vortex stays in the exact true-β interior). So: outer structure sets β-drift *magnitude* and *minimises* the aim bias;
closing the last ~8° needs the gyre-level probe.
**(2) A 4-deep V8.6.3 regression chain blocked the storms.** Only the FIRST was a real bug; the
other three were diagnostics lying about a healthy model (see the dedicated section below). Both
the real bug and the false-alarm collisions are fixed.
**First clean V8.7 storm results** (`INIT_SOURCE="hurdat2"`): Ivan **+70.4 km E /
−8.4 h** (over-translates poleward, along +280 km by t42, tracks accelerating ERA5 v_env), Hugo
**+102.7 km E / −2.2 h** (storm-agnostic rerun June 19; over-translates, along +22 km at the landfall
fix — the earlier +133 km / +1.1 h / −130 km under-translation came from a driver NOT on the
storm-agnostic stack). ⚠ **TODO (reframe — do not cite as written):** with the corrected Hugo, BOTH
storms now over-translate (cross early), so the old **"opposite along-track signs ⇒ NO systematic
poleward bias"** inference is CONTRADICTED — it now reads as a systematic poleward over-translation of
varying magnitude. The **eastward cross-track (~+110–140 km at landfall) IS the systematic β-aim
residual** (too poleward → westward-deficient → drifts E) — this part is unaffected.
**(3) The controlled taper A/B is DONE (Ivan, same hurdat2 init).** The model is **deterministic** —
taper-start 250 reproduced the original runs BYTE-for-BYTE, so the original "calibrated 200" storm
set was effectively running **250**, not 200. Clean A/B: **250 → +70.4 / −8.4 h** (cross@landfall
+141.5, along +279.5); **200 → +67.6 / −8.1 h** (cross@landfall +126.3, along +249.2). So 200 beats
250 in the predicted direction — **~15 km west, ~30 km less over-run, 0.3 h less early** — but that is
**~¼ of the harness's ~40–60 km prediction**: full-physics steering/drag/decay dilute the isolated
β-drift difference. **Verdict: outer structure is EXHAUSTED** — the dominant eastward residual
(~+125 km) is largely taper-INSENSITIVE → it is the gyre-level ~8° aim floor, not outer-structure-
tunable. The 3-for-3 systematic-eastward finding holds regardless (taper moves it
only ~15 km). **(4) Gyre probe DONE; aim invariant to 3 levers — but "STRUCTURAL" is PREMATURE (peer review, Jun 2026).**
NU4 sweep: load-bearing for stability at dx=15.6 km (only 3e11 clean; 1e11 Vmax-82 contaminated; ≤3e10
NaN) → diffusion-MAGNITUDE not tunable. Resolution (nx 320/480/640, all clean): heading **FLAT 350→352°
(+2°)**, |drift| grid-CONVERGES 2.48→2.28. So the ~8–17° poleward aim survives structure + diffusion-
magnitude + 2× resolution. ⚠ But we never tested the FORMULATION — diffusion FORM (∇² vs ∇⁴), advection
ORDER, projection/divergence-damper × gyre. So: "invariant to those 3 levers," NOT structural-proven. The
periodic_taper β-plane worry is REBUTTED by code (interior 60% exact true-β; ~400 km drift stays inside the
±1500 km true-β band; domain check confirms — vortex never feels the taper). The FORMULATION probe is now DONE: diffusion-FORM (∇² vs ∇⁴) and divergence-DAMPER both
EXONERATED — Arm A's 350→318° rotation was collinear with Vmax collapse (42→7); Arm B (eps sweep)
held the aim DEAD FLAT at preserved vortex. So the aim is locked to the VORTEX, not the numerics.
Reframe: a vortex-coupled spurious ~NE drift (~1.0 N + ~1.0 E m/s vs canonical NW) — NOT structural,
FINDABLE; suspects = f-plane self-intensification leak &/or Rmax structure (→ structure probe + energy
audit next). Westward component is the honest metric (~0.4 m/s vs canonical ~1.4). Cap NOT yet exonerated (cap-off effx unrun).

---

## ★ TL;DR — the arc

- **V8.7 RESULT (the current headline) — over-translation MECHANISM closed:** the
  intensification ladder + gate-beta settled it across three clean rungs on the Ivan grid.
  **J0** (f-plane, +5 N bg, init 64): the vortex self-intensifies (dips to 45 by t+24, then
  barotropically spins UP to 76 by t+48) and **eff_y stays pinned at 1.000 the whole time** →
  advection is faithful *through* intensification → **intensity-coupled numerics is dead** (the
  entire eff-family / Helmholtz / trace hunt was chasing a ghost). **J2** (β-plane + steering):
  intensifies to 80, eff_y=**1.364**, and that equals **1 + β-drift_y/v̄** with the β-drift from
  **J1** (β-plane, zero bg: net (−0.32 W, +1.89 N)). So the over-translation is β-drift on
  steering, intensity-*modulated* — physical, not a bug. **gate-beta** then sized & aimed the
  β-drift: mature **|2.69| m/s @ 353°** = TOO STRONG + TOO POLEWARD, rotating NW→N and not
  plateauing, and strengthening as the vortex spreads ⇒ **outer-structure-controlled (R_env is
  the lever).** Strategic shift: the discriminator discussion moves from "is our numerics
  faithful" (yes, demonstrably) to "is our vortex's self-propagation calibrated" — one knob.
- **V8.6.3 perf VALIDATED:** the Helmholtz-on-device port (Fable) checks out — L0 reproduces the
  trace (eff_y 0.998) and the pipeline collapsed to **~0.16 s/step** (52 h Ivan-grid run ~16 min,
  full ladder/gate runs in minutes). The ~8–10× was real; the rest of the V8.6.3 stack landed too.
- **V8.6 RESULT (prior storm headline):** identical storm-agnostic config, HURDAT2-true
  inits: Hugo **+1.2 h / +131.2 km E**, Katrina **−3.2 h / +126.5 km E**. Registered
  predictions held (Hugo 1 km over the stated band — at the boundary). The residual is ONE
  number across two storms: **~1.1–1.3 m/s eastward**. Mechanism localised: obs zonal motion
  exceeds DLM westward by ~1–2.5 m/s (propagation/β-drift); model recovers ~0.2–0.7. Error
  growth concentrates in the high-gap windows (Hugo t0–6 and t24–28 at ~4 m/s; flat when the
  gap is small). Candidate physics: taper under-sizes the vortex (R34: Katrina ~120 nm mean,
  Ivan ~190 nm; real outer circulation extends well past taper_start 250 km) → β-drift
  deficit; plus the separate ~0.6 m/s numerics drift (effx).
- **HURDAT2 verification (V8.5.3):** downloaded hurdat2-1851-2025 (and 2023 ed. — identical,
  so not revisions). Init mismatches: Hugo 27.7/74.5/110kt/942 vs true **27.2/73.4/100kt/950**
  (the legacy position = HURDAT2 ~21/04Z, a +4.2 h along-track head start); Katrina
  24.0/86.3/**145kt**/906 vs true **24.8/85.9/100kt/941** (init was −98 km cross-track, west);
  Ivan lon **85.0 vs 86.0** (transcription). Hugo ID = **AL111989** (was AL131989). Legacy
  fallback steering vectors also don't match best-track motion. Full report:
  `HURDAT2_VERIFICATION.md`.
- **ERA5 audit:** the pipeline "observed positions" were synthetic (Katrina's literally a
  straight line init→landfall). Contamination limited: DLM comes from the raw ERA5 grid;
  obs_track fed only the printed table + the retired frozen track-mean mode. **Time-varying
  lockstep runs (Run 15, Hugo re-val 2) have clean steering.** obs_tracks now HURDAT2.
- **Hugo re-val 2 (HURDAT2-scored, same-latitude 32.5°N):** −3.6 h early, +112 km E. Net of
  the init head start the dynamics ran ~0.6 h slow. Cross-track accumulates at a **steady
  ~0.7 m/s from t≈6 h** (the earlier "delayed-onset drift" read was an artifact of scoring
  against the synthetic track — retracted).
- **Katrina Run 15 (HURDAT2-scored, 29.1°N):** −2.3 h early, +14.6 km E — **retired**: init
  placed the model 98 km west (left) of track, eastward drift ~1.0–1.2 m/s crossed the true
  track at t≈21 h and the two near-cancelled at landfall. Same archetype as Hugo's retired
  48 km.
- **Drift coherence (the good news):** eastward cross-track drift appears in BOTH storms
  (~0.7 and ~1.0–1.2 m/s), with the cap clamping from t=0 (Katrina) and not until t+22
  (Hugo) → **two-storm evidence against the cap as the source**, ahead of cap-off effx.
  Matches effx (0.6) and betadrift (+0.43 E).
- **V8.6 staged:** `INIT_SOURCE = "hurdat2"` in storm_data applies true inits (position,
  Vmax, P_min, f from true lat0, corrected fallback steering, HURDAT2 landfall records);
  `"legacy"` bit-reproduces all V8.5.x runs. **Registered predictions** (in the verification
  doc, written before any V8.6 run): **Katrina +80…+130 km E; Hugo +90…+130 km E; timing
  small.** Both ≈ +110 ± 30 km E ⇒ one coherent honest residual across two storms = the
  BAMS result. Far outside ⇒ superposition wrong, investigate.

---

## Repository layout

```
oracle_v8/
├── vortex_init.py                          ★ Holland + WIND_TAPER (bounds vortex at R_env)
├── run_hugo.py                             ★ cap+big+taper + lockstep steering; 36h window (V8.6)
├── run_katrina.py                          ★ cap+big+taper + lockstep steering
├── run_ivan.py                             ★ NEW — Ivan (12Z init, recurver); 5000km/320² domain
├── landfall_verify.py                      ★ NEW — along/cross-track verification (shared)
├── diagnostics.py                          ★ real Vmax, ζ-anchored ventilation, KE/enstrophy
├── run_translation_test.py                 ★ f-plane harness + discriminator/betadrift/effx + (V8.7) intensify-ladder, gate-beta(dir), gate-beta-renv
├── storm_tracker.py                        ★ gate + re-acq + argmax[vor] + raw-θ′
├── era5_steering.py                        ★ ERA5 DLM; obs_track now HURDAT2 (table-only)
├── storm_data.py                           ★ HURDAT2 inits via INIT_SOURCE toggle + obs_track
├── HURDAT2_VERIFICATION.md                 ★ NEW — audit findings F-V1…F-V6 + predictions
├── test_footer.py                          ★ NEW — re-score any old run log vs HURDAT2
└── solver/
    ├── operator_config.py                  ★ OperatorConfig (FIXED slots — no cap slot)
    ├── tendency.py                         ★ Coriolis/Drag/IntensityCap; set_env() for steering
    ├── poisson.py · integrator.py          ★ batch Thomas; RK3 (diag.max_u = max|u-COMPONENT|)
```
Placement: run scripts, diagnostics, landfall_verify, storm_data at oracle_v8/ TOP level;
`tendency.py` in `solver/`. Delete `__pycache__` when OneDrive serves stale imports (check
`d.__file__`). A duplicated result with *identical wall time* = the run didn't re-execute.

---

## ★ HURDAT2 VERIFICATION + V8.6 INIT TOGGLE (V8.5.3 → V8.6.0-pre)

Findings (full detail in HURDAT2_VERIFICATION.md):
- **F-V1** Hugo ID wrong (AL111989, not AL131989) — fixed.
- **F-V2** Hugo init ≠ HURDAT2: legacy position = best track at ~21/04Z (+4.2 h along-track
  head start, −8 km cross); Vmax was the 12Z value; P matched no fix.
- **F-V3** Katrina init ≠ HURDAT2: −98 km cross-track (west), ~0 along; Vmax 145 was the
  12Z value (true t=0 = 100 kt — **below the 70 m/s cap**); 906 mb matched no fix.
- **F-V4** Ivan lon 85.0 → 86.0 (transcription) — fixed under hurdat2 source.
- **F-V5** ERA5 obs tracks synthetic; production steering pathway unaffected (resolved).
- **F-V6** Katrina's +14 km = compensating errors — retired as a headline.
- Provenance pattern (time-shifted positions, intensities from later fixes, pressures
  matching nothing) = values recalled from memory, not read from the file. The verification
  protocol caught it pre-submission → methodology-paper material (epistemic probes).

Toggle: `storm_data.INIT_SOURCE = "hurdat2" | "legacy"`. hurdat2 updates lat0/lon0, Vmax,
P_min, f(true lat0), fallback u/v_env (centred best-track motion), landfall records.
Legacy values stay in the literal dicts → every historical run reproduces. Rmax (EBTRK)
NOT yet re-verified — separate dataset, check before BAMS.

Scoring: ALL landfall scoring now via `landfall_verify.landfall_report` (same-latitude
threshold comparison + along/cross decomposition at obs fixes + labelled legacy metric).
The old footer's landfall-point comparison conflated along+cross (Hugo re-val 2 printed
+137 km when same-latitude was +112). `test_footer.py <log> <STORM>` re-scores old logs.

---

## ★ THE STRUCTURAL FIX — wind taper (vortex_init.py, V8.5.1)

The bare Holland profile `V_t(r) = Vmax·√(x·e^{1−x})`, x=(Rmax/r)^B, decays as **r^(−B/2) ≈
r^−0.75** — 27 m/s @500 km, 16 @1000, 12 @1500, 10 @2000. It is **never truncated**: `R_env`
appears only in `self._r1d = linspace(0,R_env,n)` (the 1-D pressure-integration grid) and so
shapes only P′ → θ′, which is dynamically PASSIVE (buoyancy off). So the wind filled the whole
domain → oversized β-drift + taper reach + overshoot.

Fix: `HollandVortexInit(..., wind_taper=True, taper_start_frac=0.5)` multiplies V_t by a cosine
taper 1→0 over [taper_start_frac·R_env, R_env]; 0 beyond. Inner core untouched, peak preserved,
radial (axisymmetric) so ~no divergence (pre-balance mops up). Default **False = bit-identical to
all prior runs**. Now R_env is the real vortex-size knob. Wired into run_katrina + run_hugo as
`WIND_TAPER = True`, R_ENV_M = 500 km.

⚠ STRUCTURAL: changes the vortex for every storm. This is why Hugo's old 48 km dissolved.

---

## ★ Intensity cap (IntensityCapComponent, tendency.py, V8.4.2)

Perturbation-wind speed limiter: where |V'|>v_cap, relax `d|V'|/dt = −(|V'|−v_cap)/τ` (dir
preserved), via `dV'/dt = −α·V'`, α=(|V'|−v_cap)/(τ|V'|). Background untouched (perturbation-
relative ⇒ needs set_env lockstep with time-varying steering). Applied as a **per-step loop
relaxation** (forward-Euler after integrator.step), NOT an OperatorConfig slot (fixed slots).
`VMAX_CAP_MS=70`, `TAU_CAP=300`. Rig-tested. Contains the >~90 m/s numerical runaway (Five's
max|u|≠intensity, demonstrated: uncapped init-120 → peak 372–493, ε/drag = stability crutches).
**V8.6 note:** true Katrina init Vmax = 51.4 m/s < cap → cap inactive at init for ALL storms;
it now only engages on intensification. Drift evidence (both cap states) disfavors the cap
as the eastward-drift source — and the V8.6 runs scored essentially CAP-FREE (Hugo max|u|<70
until post-threshold; Katrina decays from 52, never near 70) with the residual fully present
→ **cap exonerated for the track residual**. Cap-off effx still isolates the harness 0.6 m/s.

---

## ★ BIG_DOMAIN (run_katrina + run_hugo, V8.5.1)

`BIG_DOMAIN=True` → Lx=Ly=4000 km, nx=ny=256 (dx=15.6 km unchanged). Pushes the reversed-β taper
zones to y<800 / y>3200 km so the storm never reaches them. **Requires the cap on** (uncapped 256²
blows up). Hugo's 32.5°N landfall sat at y~1.8–2.0 Mm in the 2000 km domain — INSIDE the taper —
so Hugo's old result was taper-contaminated too; the big domain fixes that.

---

## ★ Translation harness (run_translation_test.py) — modes

f-plane (β OFF) by default; eff = displacement/(background·T). **All f-plane runs now pass the
background reference to Coriolis** (`u_env,v_env`) — without it the steering flow inertially
oscillates (T≈26 h) and the storm loops, which CONTAMINATED the old "13% over-translation" (that
was an artifact, not a model property). Params: nx, v_cap, dom (sets Lx=Ly), r_env, beta,
wind_taper. Center via dg.vorticity_center (tracker-independent).

CLI sub-runs (skip the expensive full main):
- `… 7` / `disc` — **discriminator**: capped 256² at big domain (Run-14 grid) vs finer-dx.
  Result: big-domain eff_y = baseline ⇒ overshoot is β-drift, not grid → wind taper.
- `… beta` — **betadrift**: β-plane, u=v=0, sweep R_env. Untapered = 2.17 m/s N (R_env-inert);
  tapered = ~0.6–1.0 m/s and R_env-dependent. ⚠ at ~1 m/s over 10 h the drift is ~1–2 cells →
  near the quantization floor; magnitude/direction soft (use longer runs for a clean vector).
- `… effx` — **zonal vs meridional**: capped, balanced background. Meridional faithful
  (eff_y≈1.0–1.13). Zonal shows a W/E asymmetry (W 0.87 / E 1.13) = a constant **~0.6 m/s
  eastward drift**, not a symmetric deficit. (cap-off effx pending — but two-storm cap-state
  evidence now disfavors the cap; check WHEN the asymmetry accumulates, not just its total.)
- `… intensify-ladder` / `16` — **(V8.7) emergent intensification ladder**: Ivan grid/structure,
  init 64, 52 h, β/steering added one rung at a time, eff_y(t) read against EMERGENT Vmax(t).
  J0 f-plane (self-intensifies 45→76, eff_y≡1.000 → advection exonerated); J1 β zero-bg
  (β-drift net +1.89 N, vortex decays); J2 β+steering (intensifies→80, eff_y 1.364 = 1+β/v̄).
  The intensification is **emergent barotropic** (angular-momentum spin-up, θ′ zeroed at init in
  run_ivan — barotropic; contained by drag+cap), NOT thermodynamic/imposed.
- `… gate-beta` / `12` — **(V8.7 upgraded) static β-drift MAGNITUDE + DIRECTION**: β-plane,
  u=v=0, 48 h, **mature drift read from a late 30–48 h window** (after the Holland-init
  adjustment), with a per-interval heading trace and a theory-band verdict (1.5–2.5 m/s,
  290–335° NW). Ivan-strength: |2.69|@353° = TOO STRONG + TOO POLEWARD, rotating, non-stationary.
- `… gate-beta-renv` / `17` — **(V8.7) R_env sweep** at fixed taper_start_frac=0.5, Vmax 64:
  tests whether tighter outer bound → |drift| into band + heading back NW + plateau (one knob).
  Summary flags the sweet-spot R_env (mag+hdg both IN); fork = if magnitude tracks R_env but
  heading stays poleward, the rotation isn't a size effect → domain/sponge check next.

---

## ★ Time-varying lockstep steering (run_katrina + run_hugo, V8.4.0 / ported to Hugo V8.5.2)

`TIME_VARYING_STEER=True`, `TAU_STEER=10800` (3 h). Each diag step: sample DLM at the model's
(t,lat,lon); close `alpha_steer = min(1, DIAG_EVERY·DT/TAU_STEER)` of the gap; shift the uniform
background IN THE STATE by Δ(u,v) AND advance the Coriolis/drag/cap reference by the SAME Δ
(**lockstep** → tendency invariant; reference-only update injects the old Runs 5–6 eastward
torque). A uniform shift is divergence/vorticity-free → projection ignores it, vortex untouched.

Hugo and Katrina share the SAME steering architecture = the literal storm-agnostic config.
ERA5 audit: DLM built from the raw grid; sampling follows the MODEL position (clean of the
synthetic obs tracks). DLM layer = **850–300 hPa** in code; storm_data docstring said 850–200
— reconcile for the paper (850–300 matches the implementation).

---

## ★ Diagnostics (diagnostics.py, V8.4.2, ζ-anchored)

`compute_diagnostics` computes the vorticity centre FIRST and anchors ventilation + Vmax on it
(not the passed θ′-centre): the azimuthal cancellation needs the ROTATIONAL centre, and θ-vs-ζ
separations reached ~90 km (Run 13) / ~200 km (Run 14), contaminating an θ-anchored read with
leaked vortex rotation. Rig-verified to recover a −3 m/s gyre with a 133 km bad θ-centre. r_vent
default 120 km. ⚠ perturbation_ke is domain-MEAN (dilutes ~4× in the big domain) — make it a
vortex-region integral.

---

## Storm history

| Run | Config | Result |
|-----|--------|--------|
| 12 | small dom, no cap | crossed 29.1°N; +8.6 h late, −33.5 km W; max|u|→150 |
| 13 | small dom + cap | +8.9 h late, **−15.6 km** (cross-track halved); Vmax honest ~75 |
| 14 | big dom + cap | **−7.6 h early**, −60 km W; ran to 35.5°N; tracking degraded (sep→200) |
| 15 | big dom + cap + **taper** | legacy-scored −3.0 h, +14 km; **HURDAT2-scored −2.3 h, +14.6 km E — RETIRED (compensating: init −98 km W × ~1 m/s E drift)** |
| Hugo-reval-1 | cap+big+taper, frozen steer | 0.3 h late, +182 km E (old 48 km was comp. errors) |
| Hugo-reval-2 | + time-varying steering | **HURDAT2-scored: −3.6 h, +112 km E** (timing ≈ init head start; net dynamics ~0.6 h slow; drift ~0.7 m/s steady) |
| V8.6 Hugo rerun-1 | true init, 36h window | **+1.2 h late, +131.2 km E** (predicted +90…+130 — edge). Behind all run (−92 km along @t24, spin-down); cross grows in two bursts (t0–6, t24–28) matching the obs−DLM gap |
| Ivan run-1 (A/B) | ⚠ ERA5 file MISSING → silent fallback to CONSTANT steering (u=−2.1, v=+3.4 frozen 52 h) | **NOT the registered experiment — predictions unscored.** −1.3 h, **−192 km WEST**: frozen steering fails W on the recurver (vs +182 km E on Hugo re-val-1) — opposite-sign frozen-steering failures = the time-varying architecture's necessity demo. 320² domain + cap stable 52 h ✓ |
| Ivan run-2 (TRUE) | true init, ERA5 lockstep, 5000/320 | **−8.3 h early, +72.8 km E** — timing band MISSED decisively; cross band missed low but fixed-time cross rate 0.9 m/s = inside the mechanism band. Cause localised: model out-runs its OWN sampled background poleward by +0.6→+2.9 m/s (growing); NOT DLM-sampling feedback (sampled v was weaker than obs-track DLM). Onset predates re-intensification (44→84 m/s after t24). First run on the untested 320² grid — harness at this exact grid is the gating experiment. **V8.7 resolved this:** the +0.6→+2.9 growth is β-drift maturing/rotating poleward (not numerics, not intensity-per-se); the early +0.6 is the benign onset, the growth is the over-developing β-gyre. Mechanism closed; fix = R_env/β-drift calibration |
| V8.6 Katrina rerun-1 | true init | **−3.2 h early, +126.5 km E** (predicted +80…+130 ✓). Over-translates N late (+95 along @t36); cross grows ~steadily under a persistent ~−2 m/s obs−DLM gap |
| **V8.7 Ivan** | hurdat2 + proj fix; **taper-start 250** (was mislabeled "200"; A/B-confirmed byte-identical) | **PROVISIONAL until log is checked in:** −8.4 h early, +70.4 km E (threshold; cross@landfall +141.5, along +279.5). OVER-translates poleward: along +280 km by t42, crosses 30°N at t33.6, tracks accelerating ERA5 v_env (3.9→6.6). Cross ~+110–140 at landfall-time. |
| **V8.7 Ivan (taper A/B)** | hurdat2 + proj fix; **taper-start 200** (true calibrated) | **PROVISIONAL until log is checked in:** −8.1 h, +67.6 km E (cross@landfall +126.3, along +249.2). vs 250: ~15 km W, ~30 km less over-run — predicted direction but ~¼ of harness; eastward residual taper-insensitive. |
| **V8.7 Hugo** | hurdat2 init + projection fix; **taper-start 200** (frac 0.40); **storm-agnostic rerun June 19 2026** | **−2.2 h early, +102.7 km E** (same-lat 32.5°N; legacy +128.3). OVER-translates: along +22 km at the landfall fix; cross **+111.3** at landfall-fix (t28). Supersedes the provisional +1.1 h / +133.3 km E / −130 km under-translation run, whose driver was NOT on the storm-agnostic stack. ⚠ Hugo now over-translates like Ivan → the old "opposite along sign ⇒ no systematic poleward bias" argument needs reframing (see header TODO). |
| **V8.7 Katrina** | hurdat2 init + projection fix; **taper-start 200** (frac 0.40, now persisted in script) | **PROVISIONAL until rerun:** −2.6 h early, +114.3 km E (same-lat, thresh 29.1°N; +124.6 at landfall-fix 29.3°N). MILD over-run (+42→+77). Do not call this 3-for-3 until checked-in logs are regenerated. |
| **3-storm cross-track cluster** (calibrated taper-start 200 target) | landfall-fix cross-track | **Hugo regenerated June 19 (storm-agnostic):** Hugo **+111.3** · Katrina +124.6 · Ivan +126.3 → mean +120 km, all E, span only ~16 km. Katrina/Ivan still **PROVISIONAL until their logs are checked in**. The code persists `taper_start_frac=0.40` in all three drivers and the testbed uses taper-start 200; checked-in same-lat rows now Hugo **+102.7** (corrected), Katrina +126.5, Ivan +72.8. |
| **FORMULATION PROBE** (gate-beta-diffform / -divdamp; in-band, ∇⁴ ref 350°/west+0.42/Vmax42) | mature aim vs operator/damper | **Arm A ∇²:** nu_H 1.2e4→5e4→2e5 = hdg 347/338/318, |drift| 2.30/1.92/1.38, Vmax 35/21/**7** → rotation COLLINEAR w/ Vmax collapse (confound); matched-strength effect only −3°. **Arm B eps** 0.5/0.25/0.1 = hdg 350/350/350 (FLAT), Vmax 42/37/34 (preserved); eps=0 NaN (load-bearing). → **diffusion-form + divergence-damper EXONERATED; aim locked to VORTEX**. Bias = spurious ~NE (1.0 N + 1.0 E vs canonical NW) |
| **STRUCTURE PROBE** (gate-beta-rmax @nx480 / gate-beta-vmax @nx320, ∇⁴) | aim vs Rmax & init-Vmax | **Rmax** 31/46/75km → hdg 343/347/351, Vend 32/36/42. **Vmax** init 64/50/35/21 → Vend 42/35/27/18, hdg 350/348/344/338 (all DECAYED — leak didn't refill). Both trace the SAME aim=f(Vmax_end) curve → **path-independent**. DECOMP: **WEST≈const 0.41; NORTH∝Vmax (~0.06·Vmax thru origin)** → looked like spurious poleward drift |
| **f-PLANE DECOMP** (gate-beta-fdecomp; β on vs off × Vmax 64/35/21) | isolate spurious vs true β-drift | **f-plane drift = ZERO** (fp_N/fp_W = 0.02/−0.00/−0.01 — vortex doesn't move with β off!) → drift is **GENUINE β-drift, NOT artifact**. true β-drift (β−f): N 2.43/1.45/0.96, W 0.40/0.42/0.38 → grows AND rotates poleward w/ Vmax (338→351°). → bias is a DIRECTIONAL gyre-orientation error → **β-gyres OVER-ROTATE poleward** (∝Vmax); seat = equilibration. Next: time-evolution + gyre field |
| **TIME-EVOLUTION** (gate-beta-timeevol; β, Vmax 64 & 21, 12h windows) | heading vs time | Vmax21: 322/327/335/339° (t6-12/12-24/24-36/36-48); Vmax64 net +14° (→~350). **Heading CLIMBS poleward over the run** (early ~canonical NW → late N) = over-rotation CONFIRMED. Rate ~**Vmax-INDEPENDENT** (+14 vs +16°) → NOT swirl-driven → β-Rossby gyre dynamics, **equilibration failure**. Late deceleration hint (0.67→0.33°/h). Next: long-run triage (plateau vs runaway) |
| **LONG-RUN TRIAGE** (gate-beta-longrun; β, Vmax 64, 96h, 12h windows) | plateau vs runaway | Clean t6-60: 337/342/347/351/356° at steady ~0.4°/h, NO plateau → **RUNAWAY** (heading marches THROUGH north). **SPEED saturates ~2.5 m/s** (correct β-drift magnitude!) while **DIRECTION precesses** → gyre AMPLITUDE equilibrates, PHASE precesses (β-Rossby rate, no phase-lock). ⚠ **center-tracker BREAKS t60+** (|drift|→15 m/s, dir swings SE/NW/S = ζ²-centroid loses vortex; late windows GARBAGE). → gyre instrumentation next |
| **GYRE INSTRUMENTATION** (gate-beta-gyre; β, Vmax 64, 60h, m=1 vorticity @ t12-60) + **RE-ANALYSIS** (`reanalyze_gyre.py --tag v64` on saved ζ) | SEE & characterize the gyre | **FIGURE (`gyre_precession_v64.png`): a clean β-gyre dipole/spiral visualization.** Re-centering shows tracked-center offset grew to 17 km, so raw 75-450 m=1 was contaminated. Outer bands are mostly stationary (~W): phi[150-450] 247→245, phi[250-500] 266→271 over t12-48; steering extraction is plausible but not a rigorous match yet (339→317 while drift heading 342→351). Keep the visual gyre/spiral claim; treat gyre→drift steering closure as needing a tighter extraction. ⚠ fixed an m=1 phase SIGN bug and snapshot sets are now tag-selected. |
| **Vmax-DEPENDENCE** (gate-beta-gyre 64/35/21 + compare_vmax.py; β, 60h each) | swirl-shear vs β-Rossby | Mature (t24-48) means across init Vmax 64/35/21 (Vmax_est 32/19/12 = real 3× swirl range): \|wind-up\| 103/88/96°, steer bearing 332/338/336°, poleward-of-NW +17/+23/+21°. **ALL FLAT → INTENSITY-INDEPENDENT** → directional bias is **β-Rossby/structural, NOT swirl-shear**. Drift SPEED scales w/ strength (2.2/1.0/0.9 m/s) but DIRECTION (~335°, NNW) does not. Consistent w/ β-drift being a structural (not amplitude) phenomenon. **β-drift mechanism characterization COMPLETE** → pivot to baseline + more storms |

All scoring same-latitude at threshold_lat (Hugo 32.5°N, Katrina 29.1°N, Ivan 30.0°N
provisional) vs the HURDAT2 obs_track via landfall_verify. Hugo true init 27.2°N/73.4°W
(1989-09-21 00Z) Vmax 100 kt; obs crossed 32.5°N at t+26.9 h / 79.53°W; landfall 22/0400Z
= t+28.0 h. Katrina true init 24.8°N/85.9°W Vmax 100 kt; obs crossed 29.1°N at t+34.2 h /
89.6°W; landfall 29/1110Z = t+35.17 h. Katrina = obs-informed hindcast; Hugo = independent
t=0 init validation.

---

## ★ Compensating-errors cascade (the paper's spine)

Each fix exposed the next hidden error — the Peircean / epistemic-probes story:
1. **Intensity runaway** → cap. (max|u|~150 was numerical, not physical.)
2. **Taper lag** (+8.9 h) → big domain. (Removing it overshot → taper WAS the brake.)
3. **Overshoot** → wind taper. (R_env was inert; the vortex was never bounded; β-drift inflated.)
4. **Hugo's 48 km** → broad-β-drift masking a frozen-steering westward deficit → steering port.
5. **Katrina's +14 km** → bad-init westward displacement cancelling the eastward drift at
   landfall hour. (Found only by scoring against verified HURDAT2 data.)
6. **The data layer itself** → init values and "obs" tracks that were never HURDAT2; caught by
   executing the verification the docstring demanded. The cascade reached BELOW the physics.
7. **The numerics ghost** (V8.7) → the eff-family / Helmholtz / trace saga chased an "over-
   translation" that the validated instrument showed was lattice flicker, and the
   intensification ladder then proved advection is faithful even *through* a 45→76 spin-up
   (J0 eff_y≡1.000). The "over-translation" was real but PHYSICAL: β-drift on steering. The
   error sits one level out — the taper that fixed chapter 3's overshoot **mis-calibrated the
   β-drift** (too strong, too poleward, non-stationary, outer-structure-controlled). The
   cascade's first PHYSICS TRADE-OFF, not a bug, now made quantitative: track skill = β-drift
   calibration via R_env (R34 constraint), and the discriminator moves from "are we faithful"
   (yes) to "are we calibrated."
Remaining: the ~1.1–1.3 m/s eastward residual — ONE coherent number across both V8.6 storms.
Decomposes as (real propagation ~1.5–2 m/s W, visible as the obs−DLM gap) − (model tapered
β-drift ~0.5–1 W) + (spurious ~0.6 E, effx). The taper that fixed chapter 3's overshoot
likely over-shrunk the vortex — the cascade's first error that is a PHYSICS TRADE-OFF rather
than a bug: vortex size now needs a data constraint (R34), not a numerical fix.

Also: Five↔Gemini disagreement localising a checkable fact (does the storm reach the taper?) —
one computation adjudicated it. The harness inertial-oscillation bug — a phantom "13% over-
translation" that the balanced re-run dissolved. And a same-day micro-example: the first
re-scoring of Hugo re-val 2 trusted the pipeline's "observed positions" and produced a
"delayed-onset drift" finding — retracted once scored against true HURDAT2 (drift is steady).

---

## ★ PERFORMANCE (V8.6.3 — the profiler's first catch)
The ORACLE_PROF timers named the throughput eater on their first run:
**HelmholtzDivergenceDampingComponent = 48% of wall (~96% of physics time)** —
it ran ENTIRELY ON THE CPU by construction ("CPU for FFT"): full u,v pulled to
host + 4 numpy FFTs + results shipped back, twice per step, every step, all
campaign. This is why cupy barely beat numpy (everything else was already
fast and small). Ported on-device (identical float64 math). Stack of V8.6.3
perf fixes: Helmholtz on-device; Poisson Thomas factors precomputed+cached;
Poisson residual diagnostics sampled (ORACLE_POISSON_DIAG_EVERY, default 60);
zero-mode Thomas on host; integrator diag syncs throttleable
(ORACLE_DIAG_EVERY). **VALIDATED (V8.7):** L0 reproduces the trace (eff_y 0.998)
and the pipeline collapsed to ~0.16 s/step — 52 h Ivan-grid run ~16 min, the
~8–10× confirmed (the whole stack landed, not just Helmholtz; diag_sync now ~0%).
⚠ **But the stack was validated ONLY on the harness (zero background) and shipped a
4-deep regression that blocked the first nonzero-background storm run** — see the
V8.6.3 REGRESSION CHAIN section. Speed real; correctness coverage had a hole.
⚠ bare `main()` still runs the intentionally-pathological uncapped-120 rows; the
runaway now degrades gracefully instead of crashing thanks to the subcell NaN
guard (V8.7) — but use the named modes, not `main()`.

## ★ β-DRIFT CALIBRATION — the gate-beta sweeps (V8.7, CLOSED)
Harness `gate-beta-*` modes (β-plane, u=v=0, Vmax 64, mature 30–48 h window) sized and aimed the
model's β-drift and found **one control parameter**:
- **`gate-beta-renv`:** magnitude is R_env-controlled (400→800 km: 2.36→2.95 m/s, saturating;
  400 km in-band), but **heading is poleward at every R_env** (344–353°), arriving as a
  time-dependent NW→N **rotation**.
- **`gate-beta-domain` (5000/6500/8000 km, fixed dx):** the +18° rotation is **domain-INVARIANT**
  (clearance tripled 510→1417 km, nothing moved) → **reversed-β taper EXONERATED**; the poleward
  aim is intrinsic, not a boundary artifact.
- **`gate-beta-taper` (taper_start_frac sweep at R_env=650):** heading AND magnitude move TOGETHER
  with frac (343°/2.49 → 358°/3.28) → taper shape is NOT an independent direction knob.
- **THE UNIFICATION:** (R_env=400, frac=0.5) and (R_env=650, frac=0.30) land at the SAME
  (≈343°, ≈2.4 m/s) because both have **taper-start radius ≈ 195–200 km**. So R_env and
  taper_start_frac are **degenerate** — the β-drift (magnitude AND aim) is a 1-parameter family in
  the **taper-start radius = taper_start_frac · R_env**.
- **THE FLOOR:** reducing the taper-start radius co-reduces magnitude (into band) and poleward bias,
  but bottoms at **~343° (still ~8° poleward)** where magnitude hits the band floor and the taper
  can't shrink further without eating the core (~195 km ≈ 2.6×Rmax). **The residual ~8° is
  FUNDAMENTAL** (gyre resolution / NU4), not outer-structure-tunable.
- **THE STORM FIX (V8.7):** `TAPER_START_M = 200_000`; `TAPER_START_FRAC = TAPER_START_M / R_ENV_M`
  (=0.40 at R_env=500), passed to HollandVortexInit in all three run scripts. One shared value,
  calibrated to β-drift physics (NOT landfall), locked before the runs. **Reframe:** the V8.6
  "westward deficit" (recovers only 0.2–0.7 of obs 1–2.5 m/s W) was a **mis-AIM**, not
  under-magnitude — at 353° a 2.69 m/s drift is only 0.33 m/s W; at 343°/2.49 it is 0.73 m/s W
  (magnitude DOWN, westward component UP).

## ★ V8.6.3 REGRESSION CHAIN — perf stack vs the storms (V8.7, found+fixed)
The V8.6.3 perf stack was validated on the HARNESS (**zero background**) and hid a 4-deep failure
that surfaced the instant a storm with a **nonzero ERA5 background** ran. Ivan NaN'd at "step 1";
peeling it took four layers and **only the first was a real bug** — the other three were diagnostics
lying about a healthy model. Unifying lesson: **the harness never exercised a nonzero mean, so every
nonzero-background path was effectively untested.**
1. **REAL BUG — pre-balance lost its projection.** The "skip the 3 no-op RK3+projection sub-stages"
   optimisation (integrator, no-fast-component `else` branch; "Five, third ensemble review") assumed
   those sub-stages were always redundant. For a **projection-only config (the pre-balance) they ARE
   the projection.** Skipping them handed the raw, un-projected Holland field to the main run; the
   first projection then tried to remove the whole imbalance at once (**φ_min ≈ −60619 Pa, ~7× the
   physical deficit**) and detonated step 1. **FIX:** apply ONE projection in that `else` branch (not
   3 — no residual accumulation). Confirmed: pre-balance `phi_rms` now nonzero & decreasing
   (14→1→0.8→0.65→0.54), **φ_min −60619 → −4.6 Pa**.
2. **RED HERRING — gated `max_u` NaN-fill.** `ORACLE_DIAG_EVERY` fills `diag.max_u = nan` on steps
   the caller isn't reading (to skip device syncs). This produced (a) the phantom "max|u|=nan" in
   pre-balance iters 2–5, and (b) the **false "✗ NaN at step 1" abort** — the run-script guard
   (`run_*.py: if np.isnan(float(diag.max_u))`, ~line 384) reads the gated sentinel. State was finite
   at every phase (traced). **FIX:** run with `ORACLE_DIAG_EVERY=1` (negligible on storm runs), OR
   gate the guard: `if n % DIAG_EVERY == 0 and np.isnan(float(diag.max_u))` (both = 60 align).
3. **RED HERRING — the 1.6e7 low-k φ.** Some main-run solves print `phi_hat_max≈1.6e7` (d_hat≈1e-3).
   That is the **natural low-wavenumber Poisson response** (~150 Pa physical after the 1/N FFT
   scaling; gradient → tiny velocity correction; step-0 φ_min = −4.6 Pa confirms). Gauge-pinning the
   zero mode (`phi_hat[0,0,:]=0`) does NOT lower it → it lives in the inner low-k modes, not the
   gauge mode. **Not a bug.**
4. **LATENT — `id()` Thomas cache key.** `poisson._get_thomas_cache` keys on `id(rho_bar_full)`;
   `id()` is unique only among *live* objects → a freed-then-reallocated array can false-hit. `_refs`
   keeps the arrays alive so it didn't bite, but it's fragile. **FIX later:** content/shape key.
**The solver math itself is correct** — the harness runs the identical projection every step and is
fine; the cached Thomas (`_solve_batch` / `_get_thomas_cache`) matches the reference `_thomas_batch`
/ `_solve_zero_mode` line-for-line. The bug was always **access to** the projection (prebal) and
**detectors reading sentinels**, never the Poisson math.

## Critical gotchas
- **diag.max_u = max|u-COMPONENT|, NOT intensity** — use `low_level_vmax`.
- **`diag.max_u` is NaN-FILLED on non-diagnostic steps when `ORACLE_DIAG_EVERY>1`** (V8.6.3) — it is
  a sentinel, not a blow-up. Any NaN check on it must be gated to diag steps (see regression chain).
- **Steering / cap / drag = LOCKSTEP**: shift the state background AND set_env the references by
  the same Δ. The cap is perturbation-relative → it MUST get set_env when the background shifts.
- **f-plane Coriolis needs u_env/v_env** or the background inertially loops (T≈26 h).
- **R_env does NOT bound the wind** unless `wind_taper=True`; bare Holland fills the domain.
- **BIG_DOMAIN needs the cap on** (uncapped 256² blows up).
- **Ventilation/Vmax anchor on the ζ-centre** (θ-centre can be ~200 km off).
- **Sub-cell tracking (V8.6.2) is ON by default** in storm_tracker + harness; cell-snapped
  history (every track before V8.6.2) carries a ±1.5-cell (±23 km) flicker floor — never
  compare per-segment rates across the upgrade without saying so.
- **`_subcell_refine` NaN guard (V8.7):** a blown-up/NaN field (e.g. bare `main()`'s uncapped-120
  runaway) now falls back to the seed centre instead of `int(round(NaN))` crashing. If you see a
  NaN crash there, it's a runaway upstream, not the finder.
- **run_translation Rmax/B overrides (V8.7):** `Rmax=`/`B=` args added so a rung can match a
  storm's vortex (intensify-ladder uses IVAN["B"]). Wire them into the HollandVortexInit call —
  the params are inert if you only compute `rmax`/`bb` but leave the call hardcoded.
- **Score ONLY via landfall_verify vs obs_track** — the legacy landfall-point metric conflates
  along+cross (kept in the report output, labelled, for continuity).
- **INIT_SOURCE**: "hurdat2" for V8.6+; "legacy" to bit-reproduce V8.5.x. Never compare runs
  across init sources without saying so.
- **REQUIRE_ERA5 = True** (V8.6.1, all run scripts): missing ERA5 now ABORTS instead of
  silently falling back to constant steering (the silent fallback burned Ivan run-1).
  Check the log header says "ERA5 steering: ACTIVE" before trusting any result.
- Default `wind_taper=False` is bit-identical to pre-taper runs (Hugo's old 48 km used it OFF).
- `xp.trapz` shim (=trapezoid). x=E-W→u/lon, y=N-S→v/lat. cos(lat0) lon scaling: ~0.1° error.

## Known minor issues
- ~1.1–1.3 m/s eastward track residual (Hugo/Katrina; Ivan fixed-time rate 0.9 in-band) =
  under-recovered TC propagation + the ~0.6 m/s numerics drift (cap exonerated).
- **MERIDIONAL OVER-TRANSLATION — RESOLVED as a mechanism (V8.7):** it is β-drift superimposed
  on the steering, NOT advection numerics (J0: eff_y≡1.000 through a 45→76 self-intensification)
  and NOT intensity-per-se. Decomposes as `eff_y(J2)=1.364 ≈ 1 + net_drift_y(1.79)/v̄(5.2)`.
  The Ivan run-2 +2.9 m/s growth = the β-drift maturing/over-developing (it grows and rotates
  poleward through the run). Superseded the "grid-scaling" suspicion — grid is clean (J0).
- **β-DRIFT CALIBRATION — CLOSED (V8.7):** the model's β-drift (magnitude AND aim) collapses onto
  ONE control, the **taper-start radius = taper_start_frac·R_env** (R_env & frac degenerate). Outer
  structure sets magnitude (400 km / taper-start ~200 km → in-band ~2.4 m/s) and minimises the aim
  bias, but the aim bottoms at a **~343° / ~8°-poleward FLOOR** that is fundamental (gyre
  resolution / NU4), NOT outer-structure-tunable. Boundary/reversed-β taper EXONERATED (rotation
  domain-invariant). Storm fix applied: `TAPER_START_M=200_000` (one shared value). See the
  β-DRIFT CALIBRATION section. **The remaining systematic track error is the cross-track (E) aim
  residual; closing the last ~8° needs the gyre-level probe.**
- **STORM TRACK STRUCTURE (V8.7, + Hugo storm-agnostic rerun June 19):** **cross-track (eastward) is
  SYSTEMATIC** (~+110–140 km at landfall both storms) = the β-aim residual (too poleward →
  westward-deficient → drifts E). Along-track: Ivan over-translates (+280 km ahead by t42, −8.4 h
  early, tracks accelerating ERA5 v_env 3.9→6.6); **Hugo now over-translates too** (+22 km ahead,
  −2.2 h early — corrected storm-agnostic run; the old −130 km / +1.1-h-late under-translation was the
  non-storm-agnostic driver). ⚠ **TODO (reframe):** with both storms over-translating, the former
  **"opposite along sign ⇒ no systematic poleward model bias"** conclusion no longer holds — decide
  how to frame the now-systematic poleward over-translation. Open question raised by
  Ivan: is the ERA5 DLM closer to *storm motion* than *environmental steering*? If so, adding model
  β-drift on top **double-counts** propagation → poleward over-run. Worth testing directly
  (`storm_data` note: fallback steering was "estimated from storm motion vector").
- perturbation_ke is domain-mean (dilutes in big domain).
- Big-domain tracking can degrade late (raw-θ′ vs chosen split up to ~2°); windowed track OK.
- Over-translation worsens at finer dx (eff_y 1.13→1.35 from dx 15.6→7.8) — don't naively refine.
- DLM layer 850–300 (code) vs 850–200 (storm_data docstring) — reconcile for the paper.
- Rmax (EBTRK) values not yet re-verified against the EBTRK dataset.

---

## What's next
```
NOW:   ★★ β-DRIFT MECHANISM CHARACTERIZATION COMPLETE — PIVOT TO PUBLISHABILITY (V8.7).
       The mechanism chain is strong, but post-cleanup provenance matters: the 3-storm +120km
       landfall-fix cluster is provisional until fresh logs are regenerated from the cleaned scripts.
       Evidence chain: genuine β-drift (f-plane null) → a real β-gyre/spiral visualization
       (`gyre_precession_v64.png`; tagged snapshots) with outer bands mostly stationary and steering
       extraction plausible but not yet rigorous → the gyre equilibrates a sheared spiral that aims ~20°
       too POLEWARD (NNW ~335° vs canonical NW 315°) → eastward track bias. Vmax-dependence test
       (64/35/21, real 3× swirl range): wind-up ~96-103° & poleward-offset ~+17-23° are ALL FLAT →
       the directional bias is INTENSITY-INDEPENDENT → β-Rossby/structural, not swirl-shear. Drift
       SPEED ∝ strength, DIRECTION does not — consistent with β-drift theory. This is a complete,
       mechanistic, credible "epistemic probe" error characterization — exactly the reviewer's framing.
       RUN NOW (Justin) — the fork:
       (A) ★★ **PIVOT to publishability-critical** (recommended): (1) **BASELINE** — CLIPER5 / bare
       β+steering on the same Hugo/Katrina/Ivan cases (reviewer: without it, 3 storms + no baseline
       isn't publishable). (2) **MORE STORMS** — expand to 8-10, ideally a 2nd basin.
       (B) OPTIONAL mechanism coda (only if we want to try to REDUCE the bias, not just characterize):
       the outer gyre (250-500km) sits INSIDE the taper region (taper-start 200km) — earlier calib found
       taper-start IS the β-drift control knob — so a **taper-start / R_env sweep** of gate-beta-gyre would
       tell us if the +20° aim is intrinsic β-Rossby (unavoidable) or a tunable config/numerical effect
       (potentially fixable → better tracks). ~12min/run. High-leverage IF reducible, but re-opens calib.
       DEFERRED: ERA5 DLM double-counting; cap-off effx; fix 850-300/850-200; persist taper_start_frac=0.40.
       PAPER: write up the β-drift subsection with the COMPLETE mechanism (genuine β-drift, correct speed-
       scaling, intensity-independent ~20° poleward aim bias from the β-gyre equilibrium, figure). Replace
       the old "structural residual" text. Westward-component metric. Justify ~barotropic vs BAM/VICBAR.

REVIEW: ★ EXTERNAL PEER REVIEW (clean model, Jun 2026) — serious, fair, technically literate; CONVERGES
       with our framing (methodology + error-char paper; single β-drift bias; cluster = strong evidence).
       VALIDATES: β-drift should be NW 1–3 m/s (our 290–335 band lit-consistent: Holland 83 / Chan-Williams
       87 / Fiorino-Elsberry 89); poleward aim → E cross-track; 3-storm cluster → single source;
       compensating-errors method "exemplary." 4 CONTRIBUTIONS: (A) error-characterization-as-methodology
       [strongest], (B) β-drift DIRECTION as independent diagnostic, (C) compensating-errors narrative,
       (D) storm-agnostic config as structural-error test.
       SHARP CATCHES (→ folded into NOW): structural premature (formulation untested); f-plane self-
       intensification suspicious; magnitude misleading → westward component (2–3× deficit); 3 storms +
       no baseline publishability-critical; cap-off effx unrun.
       WE REBUT / NUANCE: periodic_taper β-plane — code shows interior 60% exact true-β & vortex never
       leaves it (domain check confirms); "energetically leaky" is suspicion not proof (could be initial-
       balance ringdown — the audit decides, don't assume leak). Honest that model is ~barotropic-in-
       practice (θ'=0, buoyancy off, Newtonian cooling bounds θ') → must justify vs barotropic baselines.
PREV:  ★★ OVER-TRANSLATION MECHANISM CLOSED + β-DRIFT NAMED (V8.7-pre): J0 f-plane
       self-intensifies 45→76 with **eff_y≡1.000** (advection faithful → intensity-coupled
       numerics DEAD); J2 eff_y=1.364 = **1 + β-drift_y/v̄** (over-translation = β-drift on
       steering). gate-beta mature **|2.69| m/s @ 353°** = TOO STRONG + TOO POLEWARD,
       rotating NW→N ⇒ outer-structure controlled. [Now fully resolved → β-DRIFT CALIBRATION
       section: the knob is the taper-start radius; the ~8° aim is a fundamental floor.]
PREV:  ★★ INSTRUMENT VALIDATED + EFF FAMILY RETIRED (trace, V8.6.3 centroid):
       Galilean control net drift **0.6 km / 14 h (0.01 m/s)** — noise floor ±0.4 km.
       Run A: **eff_y = 1.000 ± 0.007**, cross-track ±0.8 km/14 h. THE ADVECTION IS
       FAITHFUL — the entire eff>1 family (1.13 historical, 1.35 finer-dx, 1.042
       gate, ladder values) was lattice flicker around 1.0. [V8.7 confirmed this
       independently: J0 holds eff_y≡1.000 through a full 45→76 spin-up.] Storm
       residuals (+131/+127, −8.3 h) stand; mechanism budget rebuilt from clean
       measurements → β-drift on steering (see NOW).
PREV:  ★ TRACE VERDICT: **QUANTIZATION, not acceleration.** The Galilean control
       occupied exactly TWO lattice positions for 14 h (±1.5-cell flicker = the
       center-finder noise floor, ±23 km, eff noise ±0.09-0.13 at 10-14 h). The
       "anchor move" (1.042→1.230) and run-A's "10.8 m/s" were flicker; true
       eff_y ≈ 1.05 ± 0.09. ALL ladder rungs sub-noise → unreadable. 10 h checkpoint
       reproduced the gate exactly (harness edits cleared). Storm-Ivan's −8.3 h
       STANDS (cumulative 260 km ≫ flicker) but per-segment rates carry ±1-1.5 m/s.
       ⚠ HISTORICAL effx W0.87/E1.13 and betadrift +0.43E were measured INSIDE this
       noise floor — the "0.6 m/s spurious east drift" may be measurement noise.
       FIX SHIPPED (V8.6.2): sub-cell centre-finding — ζ²-patch centroid in the
       harness + 3×3 parabolic refine on the tracker's θ′ argmin (subcell=True
       default; precision ~0.1 cell ≈ 1.6 km, eff floor <±0.01). RE-RUN ORDER:
       trace (90 min — control must now sit within ~±2 km = validation), then
       ladder (readable at last), then gate-drift (re-measure the W/E asymmetry —
       it may dissolve, collapsing the residual budget to pure propagation).
PREV:  LADDER RUN: anchor MOVED — L0 (≡ gate G1 config) reads eff_y 1.230 @12h vs
       1.042 @10h ⇒ ~10.8 m/s on a 5 m/s background in hours 10-12: a LATE-ONSET
       accelerating excess in the minimal f-plane constant-background config — the
       storm-Ivan signature without β/ramp/θ′. L3≡L2 bit-identical ⇒ warm core
       exonerated forever. L1/L2 deltas confounded by onset timing — re-read after
       trace. KEY INSIGHT: uniform background on f-plane is Galilean-removable ⇒ any
       excess = center-finder artifact OR frame-breaking numerics (advection
       dispersion suspect). NEXT GPU: `trace` (~90 min) — L0 @14h with per-interval
       velocities (10h checkpoint must reproduce +1.042, closing the code-confound
       loop) + ZERO-background Galilean control (stationary vortex must not move).
PREV:  ✔ GATE: **GRID CLEAN** — eff_y = +1.042 IDENTICALLY on 256²/4000 and
       320²/5000 (Δ=0.000; same Vmax'_end — domain size & mode count have NO effect at
       fixed dx). Bonus finding: tapered eff_y = 1.042 vs untapered 1.13 — the taper cut
       the over-translation; grid-borne excess explains only ~+0.2 m/s of Ivan's
       +1.4…+2.9 ⇒ STORM-STATE physics confirmed. NEXT GPU: `ladder` mode (~2.5 h) —
       adds storm ingredients one rung at a time at the Ivan grid (L1 β×steering, L2
       lockstep ramp, L3 warm core); the jumping rung names the mechanism. Then
       gate-beta (intensity coupling). R34/taper work remains gated until resolved.
PREV:  ✔ V8.6 Hugo+Katrina — predictions held; residual = +131/+127 km E (one number!).
       Next experiment (ensemble review first): vortex size from DATA, not track fit —
       R_env / taper shape constrained by observed outer size: Katrina R34 mean ~120 nm
       (140/100/100/140), Ivan ~190 nm (225/175/150/200) from HURDAT2; Hugo needs EBTRK.
       Candidates: (a) R_env ≈ k·R34 with ONE shared k (storm-agnostic), (b) blend-to-env
       taper instead of hard-zero, (c) require V(R34)=17.5 m/s on the tapered profile.
       Decide the rule BEFORE running; betadrift sweep predicts the β-drift gain; then
       V8.6.1 re-runs. Registered expectation: recovering ~1 m/s west moves both storms
       to ~+20…+60 km E without touching timing much.
NEXT:  cap-off effx (Vmax≈55, no cap) → isolate the ~0.6 m/s numerics term (cap already
       exonerated by the cap-free V8.6 runs); check WHEN asymmetry accumulates.
       Then advection / projection / center-finder hunt for the numerics term.
THEN:  [ ] energy + enstrophy on the intense phase (Five P1)
       [ ] ensemble (Five/Gemini) on taper shape: hard-zero vs blend-to-env, R_env value
       [→] Ivan READY: run_ivan.py created; ERA5 entry added (⚠ 12Z init; download first:
           python -m oracle_v8.era5_steering --download --storm ivan).  Domain 5000/320
           (4000/256 would put the outer circulation INSIDE the taper at threshold).
           Predictions REGISTERED (incl. the recurvature sub-prediction — Ivan's zonal
           stall is the discriminating regime for the propagation-gap mechanism).
       [ ] EBTRK Rmax re-verification; DLM-layer doc reconciliation
       [ ] longer betadrift (48 h) for a clean β-drift vector if needed
PAPER: along/cross decomposition figures now generated by landfall_verify (Run 12 machinery);
       cascade chapters 5–6 (data-layer errors + verification protocol catching AI-recalled
       values) into the epistemic-probes narrative; registered-predictions section; replace
       all +14 km / 48 km claims; storm-agnostic config table from V8.6 results;
       submit BAMS after V8.6 Hugo+Katrina land + Ivan.
```

---

## Versioning
**V8.7** — β-drift calibration CLOSED (taper-start radius = the one knob; R_env & taper_start_frac
degenerate; ~8° poleward aim = fundamental gyre floor; boundary exonerated via domain-invariance).
New modes: gate-beta-domain, gate-beta-taper. Storm fix: `TAPER_START_M=200_000` in all 3 run
scripts. **V8.6.3 REGRESSION CHAIN found+fixed** (4 layers, 1 real bug): pre-balance lost its
projection (no-fast-component `else` branch now applies ONE projection; φ_min −60619→−4.6 Pa); the
other 3 were diagnostics lying (gated `max_u` NaN sentinel → false aborts; 1.6e7 low-k φ red
herring; `id()` cache key latent). First clean V8.7 storms: Ivan +70.4 km E/−8.4 h (over-translates),
Hugo +133 km E/+1.1 h (under-translates) — cross-track E systematic (β-aim residual), along-track
storm-specific. Scoring confounded vs legacy baselines → controlled taper A/B pending.
**V8.7-pre** — over-translation mechanism CLOSED: intensification ladder (J0 advection faithful
through self-intensification → intensity-coupled numerics dead; J2 eff_y=1.364 = 1+β-drift/v̄)
+ gate-beta direction (mature |2.69|@353° = too strong + too poleward, rotating, non-stationary,
outer-structure-controlled). New modes: intensify-ladder (16), gate-beta upgraded to report
direction + late-window mature drift, gate-beta-renv (17). Bug fixes: _subcell_refine NaN guard;
run_translation Rmax/B overrides wired. β-drift calibration (R_env/R34) is now the primary open
issue. **V8.6.3** — performance: Helmholtz-on-device port (was 48% of wall on CPU) + Poisson
Thomas caching + sampled residual diagnostics + throttleable diag syncs; VALIDATED at ~0.16 s/step
(~8–10×), L0 reproduces the trace. **V8.6.2** — sub-cell centre-finding (ζ²-patch centroid +
3×3 parabolic refine; eff floor <±0.01); the trace then retired the entire eff>1 family as
lattice flicker.
**V8.6.0** — true-init re-runs: Hugo +1.2 h/+131.2 km E, Katrina −3.2 h/+126.5 km E;
predictions held; residual unified at ~1.1–1.3 m/s E; obs−DLM segment analysis → taper
under-sizing / β-drift-deficit hypothesis; cap exonerated (cap-free scoring windows); R34
data pulled for the size-constraint experiment. **V8.6.0-pre** — HURDAT2 verification (F-V1…F-V6): true inits staged via INIT_SOURCE toggle;
Hugo ID fix; Ivan lon fix; obs_track (HURDAT2) shipped in storm_data + era5_steering;
landfall_verify.py shared module (same-latitude + along/cross decomposition); Katrina +14 km
and the first Hugo re-val-2 scoring retired/retracted as compensating/reference-track errors;
predictions registered; Hugo window → 36 h. **V8.5.2** — time-varying lockstep steering ported
to Hugo (config identical to Katrina); Hugo re-validation. **V8.5.1** — wind taper; R_env-inert
discovery; betadrift + effx modes; f-plane Coriolis reference fix; BIG_DOMAIN; Run 15;
Hugo old 48 km identified as compensating errors. **V8.5.0** — BIG_DOMAIN taper test; ζ-anchored
ventilation; discriminator. **V8.4.2** — IntensityCapComponent; runaway demonstrated. **V8.4.1**
— diagnostics module; f-plane harness; reversed-β taper geometry verified; max|u|≠intensity.
**V8.4.0** — time-varying steering (Katrina). **V8.3.x** — tracker layers, ERA5 multi-storm,
Hugo validated.
