# HURDAT2 Verification Report — storm_data.py audit

**Date:** June 9, 2026
**Method:** Downloaded NOAA NHC HURDAT2 Atlantic database
(`hurdat2-1851-2025-02272026.txt`, the current edition) and the prior
edition (`hurdat2-1851-2023-051124.txt`); extracted best-track blocks for
Hugo, Katrina, Ivan; diffed against `storm_data.py`. The two editions are
identical for all three storms in the simulation windows — these are
stable values, not recent revisions.

This executes the check the storm_data.py docstring required:
"Verify all values below against the actual HURDAT2 file before
submitting to BAMS." **The check fails for all three storms.**

---

## Findings

### F-V1 (HIGH) — Hugo storm ID wrong
`storm_data.py` says `AL131989`. AL131989 is an unnamed October tropical
depression. **Hugo is AL111989.** (Fixed in the patched storm_data.py.)

### F-V2 (HIGH) — Hugo init does not match HURDAT2
| field | storm_data.py | HURDAT2 21/0000Z |
|---|---|---|
| position | 27.7°N 74.5°W | **27.2°N 73.4°W** |
| Vmax | 110 kt | **100 kt** |
| P_min | 942 mb | **950 mb** |

The storm_data position matches HURDAT2 interpolated to ~21/0400Z
(a +4 h offset). 110 kt is the 21/1200Z value. 942 mb matches no fix.
Init offset vs the true t=0 fix: **122 km, almost purely along-track
(+122 along / −8 cross) = a +4.2 h head start** at the observed 8 m/s
translation speed.

### F-V3 (HIGH) — Katrina init does not match HURDAT2
| field | storm_data.py | HURDAT2 28/0000Z |
|---|---|---|
| position | 24.0°N 86.3°W | **24.8°N 85.9°W** |
| Vmax | 145 kt | **100 kt** |
| P_min | 906 mb | **941 mb** |

145 kt is the 28/1200Z value (25.7°N 87.7°W, 909 mb). 906 mb matches no
fix (Katrina's minimum was 902 at 28/1800Z). The init position is ~98 km
from the true fix, roughly *behind* along track. HURDAT2 landfall:
29.3°N 89.6°W at 29/1110Z = **t+35.17 h** (storm_data: 29.1°N, t+35.0).

### F-V4 (MEDIUM) — Ivan init longitude off by 1.0°
storm_data: 22.9°N **85.0°W**, 130 kt, 928 mb. HURDAT2 14/1200Z:
23.0°N **86.0°W**, 125 kt, 930 mb. Looks like an 85.0/86.0 transcription
slip; ~102 km eastward init error if run as-is. HURDAT2 landfall:
30.2°N 87.9°W at 16/0650Z = t+42.83 h (storm_data: 30.3°N 87.7°W, t+42.0).

### F-V5 (RESOLVED → LOW) — ERA5 "observed positions" are synthetic
Audit of `era5_steering.py` (June 9): the obs tracks are hardcoded in
`STORM_CONFIGS`. **Katrina's is a straight line from the (incorrect)
init point to the landfall point — every fix matches linear
interpolation exactly.** Hugo's is a smooth synthetic track between the
same kind of endpoints. These were sanity-check scaffolding that became
mislabeled as observations downstream (including in run logs and in the
first re-scoring of Hugo re-validation 2).

**Contamination scope is limited:** the DLM is built purely from the raw
ERA5 grid; `obs_track` feeds only `print_summary`, the `--levels` CLI,
and `get_track_mean_dlm` — the frozen track-mean steering mode retired
at V8.4.0. The time-varying lockstep pathway samples `get_dlm` at the
MODEL position and never touches `obs_track`. So Run 15 and Hugo
re-validation 2 steering is clean; only pre-V8.4 frozen-mean Katrina
runs ingested the synthetic track. Both `obs_track` dicts now carry
HURDAT2 fixes and the summary header is relabeled.

Residual second-order effect: the t=0 DLM is sampled at the (wrong)
init position; the 3–7° ring average makes this a small perturbation,
and it self-corrects under V8.6 true inits.

Side nit: storm_data's docstring says DLM = 850–200 hPa; the code uses
850–300. Pick one for the paper (850–300 matches the implementation).

### F-V6 (HIGH) — Katrina Run 15's +14 km is a compensating-errors result
Re-scored against HURDAT2 (same-latitude, 29.1°N): **−2.3 h early,
+14.6 km east.** But the decomposition shows the init error placed the
model **−98 km cross-track (west, left of track)** with ~0 along-track
offset, and the model then migrated east relative to the observed track
at ~1.0–1.2 m/s throughout the run: cross-track −115 → −71 → −47 → −15
→ +30 km over t=0→30 h, crossing the true track near t+21 and landing
+15 km east. The celebrated +14 km is thenear-cancellation of (a) the bad
init's westward displacement and (b) the eastward drift — the same
archetype as Hugo's retired 48 km. **The +14 km is retired as a
validation headline** (the timing and the due-north final approach
remain genuinely good).

Cross-storm coherence: Hugo's drift runs ~0.7 m/s with the intensity
cap NOT clamping until t+22; Katrina's ~1.0–1.2 m/s with the cap
clamping from t=0 (max|u| starts above the 70 m/s ceiling). Drift
present in both cap states → further evidence against the cap as the
drift source, ahead of the cap-off effx discriminator.

---

## Consequences for results

1. **Hugo re-validation 2, re-scored against HURDAT2** (same-latitude,
   32.5°N): **−3.6 h early, +112 km east** (vs −0.7 h / +91 km against
   the pipeline track). But the −3.6 h is almost entirely the +4.2 h
   init head start (F-V2): net of init offset, the dynamics ran
   ~0.6 h *slow*. The init offset is along-track only, so the
   **+112 km east cross-track residual is real and survives.**
2. **Retraction:** the "delayed-onset eastward drift" inference from the
   pipeline-track decomposition was an artifact of the wrong reference
   track. Against HURDAT2 the cross-track error accumulates at a steady
   **~0.68 m/s from t≈6 h onward** (+32 → +76 km over t=6→24) —
   consistent with the known ~0.6 m/s spurious drift as a constant-rate
   bias, and matching effx (0.6) and betadrift (+0.43 E) estimates.
3. **Cap evidence:** the drift accumulates from t≈6 h while the
   intensity cap is not clamping (max|u| ≈ 50–55 < 70 until ~t+22).
   Leans against the cap as the drift source, ahead of the cap-off effx
   discriminator.
4. **Katrina Run 15 (−3.0 h / +14 km)** must be re-scored against the
   HURDAT2 obs_track (init was ~98 km off, roughly behind along-track —
   so the −3.0 h early is *more* anomalous than it looked, and the
   +14 km cross-track is luckily clean because Katrina's final approach
   runs due north along 89.6°W). Needs the Run 15 track log.
5. **Paper claim:** "initialized from HURDAT2 best track" is currently
   not accurate as written. Options: (a) re-run the storms with
   HURDAT2-true init values — predicted to mostly remove Hugo's timing
   error and leave the ~0.6 m/s east drift as the single honest
   residual; (b) keep the runs and re-frame init provenance explicitly.
   (a) is strongly preferred for BAMS, and is cheap relative to the
   campaign already run.

## Recommended actions

1. Ensemble review of F-V2/F-V3/F-V4: adopt HURDAT2-true init values
   (position, Vmax, P_min, f from true lat0) as a V8.6 config; re-run
   Hugo + Katrina. This is a data correction, not parameter tuning —
   the no-retuning rule is not violated.
2. Audit `era5_steering.py` + its track inputs (F-V5); re-extract DLM
   anchored to HURDAT2 positions if contaminated.
3. Trace the provenance of the bad init values. The pattern (positions
   matching time-shifted interpolants, Vmax values from later fixes,
   pressures matching no fix) is characteristic of values recalled from
   memory rather than read from the file — worth documenting in the
   paper's methodology narrative as a working example of why the
   verification protocol exists.
4. Ivan: fix the longitude before its first production run (F-V4).
5. All future landfall scoring via `landfall_verify.py` against the
   HURDAT2 `obs_track` now shipped in `storm_data.py`.

## Registered prediction (before the V8.6 re-runs)

With HURDAT2-true initializations (position, Vmax, P_min, f from true
lat0) and no other changes, the init-offset cancellations disappear and
the ~0.7–1.2 m/s eastward drift is exposed in both storms. Predicted
same-latitude results:
  * **Katrina: +80 to +130 km east** (≈ +15 observed + ~98 km of
    un-cancelled init offset), timing within ~±1.5 h of −2 h.
  * **Hugo: +90 to +130 km east**, timing within ~±1.5 h (the +4.2 h
    init head start removed).
Both storms converging on ~+110 ± 30 km east would make the spurious
eastward drift the single, coherent, honest residual of the
storm-agnostic config — and the cap-off effx test then decides whether
it is removable or reported as the limitation. If the re-runs land far
outside these windows, the superposition assumption is wrong and
something else is in play.


## Prediction scorecard (V8.6 rerun-1, June 10 2026)

| storm | predicted | actual | verdict |
|---|---|---|---|
| Katrina | +80…+130 km E, timing ≈ −2 ± 1.5 h | **+126.5 km E, −3.2 h** | inside band; timing at the edge |
| Hugo † | +90…+130 km E, timing small | **+102.7 km E, −2.3 h** | inside the band; timing −2.3 h (early), just beyond the ±1.5 h window |

† **Correction (June 19 2026).** The original Hugo entry (**+131.2 km E, +1.2 h**)
was scored on a `run_hugo.py` that was NOT on the storm-agnostic stack used by
`run_katrina.py` / `run_ivan.py`. With the driver corrected (its cap + wind-taper +
lockstep-steering settings are now bit-identical to run_katrina — see `run_hugo.py`),
the storm-agnostic rerun gives **+102.7 km E, −2.3 h** (same-latitude 32.5°N; obs
crossed t+26.9 h / 79.53°W, V8 t+24.7 h / 78.43°W; legacy landfall-point metric
+128.3 km; landfall-fix cross-track +110.2 km at t+28). This moves Hugo from a
boundary hit to comfortably inside the +90…+130 km band. ⚠ It also **flips the
along-track sign**: the corrected run crosses 2.3 h *early* (along +23 km at the
landfall fix) where the old run was +1.2 h late (along −130 km).

The superposition assumption held: removing the init artifacts exposed the same
eastward residual in both storms (Δ between storms ≈ 5 km). Segment analysis
(obs zonal motion vs ring-averaged DLM at the obs fixes) shows the observed
storms out-run the DLM westward by ~1–2.5 m/s and the model recovers only
~0.2–0.7 m/s — consistent with (tapered β-drift) − (0.6 m/s numerics drift).
Residual is therefore reframed: predominantly a **vortex-size / β-drift
representation question** (the wind taper trade-off), with the numerics drift
as the smaller second term. Next registered expectation: a data-constrained
vortex size (R_env or taper shape from R34) recovering ~1 m/s westward
propagation moves both storms to ~+20…+60 km E.


## Registered prediction — Ivan first run (before any Ivan GPU time)

Config: identical physics (cap+taper+lockstep, INIT_SOURCE="hurdat2"), Ivan
domain 5000 km / 320² (same dx; taper-zone clearance ~220 km at threshold).
Scored at 30.0°N (registered threshold; obs crossing t+42.0 h, 87.9°W).

* **Cross-track: +120…+200 km E** — IF the propagation-gap mechanism
  generalizes at ~1.0–1.3 m/s over the ~42 h scored window. Wider band than
  Hugo/Katrina because Ivan is the slowest mover (β-drift relatively more
  important) and the only recurver.
* **Timing: −3…+2 h.**
* **Sub-prediction (the mechanism test):** Ivan's observed zonal motion runs
  ~1.9–2.4 m/s westward through t+24, stalls (~0) by t+30–36, then turns
  EAST before landfall. If the residual is under-recovered westward
  propagation, the eastward error should accumulate EARLY (t0–24) and
  flatten — or partially reverse — after the zonal stall. If instead the
  error keeps growing at a constant rate through the stall, the
  constant-rate numerics drift is bigger than we think and the taper
  hypothesis is weakened. Ivan is the discriminating regime the first two
  storms couldn't provide.

Failure of the cross-track band with success of the sub-prediction (or vice
versa) is informative — score them separately.


### Ivan run-1 (June 10): NOT scored against the registered predictions
The ERA5 file was missing; the script silently fell back to constant
HURDAT2-motion steering (u=−2.1, v=+3.4 m/s frozen for 52 h) — a different
steering architecture from the registered config. Result (−1.3 h, −192 km
WEST at 30.0°N) is logged as an accidental frozen-steering A/B: the error
accumulates westward from t+24 onward, accelerating exactly as observed
Ivan's zonal motion stalls and reverses (obs lon 88.2→87.9→87.7°W while the
model marches to 90+°W under the frozen westward vector). Paired with Hugo
re-validation-1 (frozen steering → +182 km EAST), frozen steering fails in
OPPOSITE directions depending on how the storm's environment evolves —
demonstrating the time-varying lockstep architecture is necessary, not a
Katrina-specific fix. Infrastructure note: first 320² run; domain, FFT,
cap, and tracker all stable through 52 h (late raw-θ′ separation is
post-threshold). Predictions remain OPEN for the true Ivan run after
`python -m oracle_v8.era5_steering --download --storm ivan`. All run
scripts now carry REQUIRE_ERA5=True so this failure mode aborts loudly.


### Mid-run refinement — Ivan true run (registered while the run is in progress,
### from the ERA5 header table only; no track data seen)

Computed obs−DLM gaps from the live header (HURDAT2 motion vs DLM at obs fixes):
zonal gap is **steady −1.0 to −1.8 m/s (mean −1.33) through t+30 — including
into the stall**. The recurvature is the DLM turning east (+0.60 by t+30)
against *persistent* westward propagation; propagation does NOT die at the
stall. Meridional gap ≈ **0.0 ± 0.5 m/s** — the DLM carries Ivan's northward
motion almost perfectly.

**Correction to the original sub-prediction, registered before results:** the
"error flattens after the stall ⇒ taper" clause assumed propagation dies at
the stall. It doesn't. Under the propagation-deficit hypothesis the error
rate ≈ (real propagation − model propagation), independent of DLM evolution —
so a steady gap predicts near-CONSTANT eastward growth all run, same as the
numerics-drift hypothesis. Ivan's growth-timing discriminator is therefore
weaker than originally claimed. What Ivan still discriminates is MAGNITUDE:

* numerics drift alone (~0.6 m/s): **≈ +90 km E** at threshold
* propagation deficit (gap 1.33 − recovery 0.2…0.7, bracketing Hugo/Katrina):
  **+92…+164 km E**

Refined registration: **cross-track +90…+165 km E** (center ~+125; refines the
pre-run +120…+200 band downward — Ivan's gap is gentler than Katrina's −2.0);
**timing −1.5…+1.0 h** (tightened from −3…+2 via the ≈0 meridional gap);
growth near-constant throughout, including through the stall. Scoring rule:
≲ +90 km ⇒ numerics-drift-only is sufficient; ≳ +120 km ⇒ propagation deficit
dominant; +90…+120 km ⇒ gray zone, decompose per-segment. Score against BOTH
the original and refined bands and report both.


## Ivan true run — scorecard and analysis (Ivan_run_2)

**Result (same-latitude, 30.0°N): −8.3 h early, +72.8 km east.**

Scored against the registrations, no excuses first:
* **Timing: MISS, decisively.** Registered −1.5…+1.0 h (refined) / −3…+2 h
  (original); actual −8.3 h.
* **Cross-track band: MISS LOW.** Registered +90…+165 (refined) /
  +120…+200 (original); actual +72.8. However the fixed-time cross-track
  rate is ~0.9 m/s (−24 → +103 km over t0–30), **inside the mechanism band
  0.63–1.13 m/s** — the same-latitude shortfall comes from the early
  crossing truncating accumulation plus recurvature geometry splitting the
  two metrics. The zonal/propagation-deficit account survives at the
  mechanism level; the band as written did not.
* **Sub-prediction (growth profile):** cross-track growth roughly constant
  through the stall — consistent with the (corrected) registration; weak
  discriminator as anticipated.

**Where the −8.3 h comes from (log-internal analysis):**
1. NOT DLM-sampling feedback. The v sampled at the model's displaced
   position was slightly WEAKER than the obs-track DLM (−0.1…−0.6 m/s) —
   the displacement damped, not amplified. Hypothesis disproven by data.
2. The model out-translated ITS OWN sampled background poleward by a
   growing excess: +0.6 (t12–18) → +1.4 → +2.2 → +2.7 → **+2.9 m/s**
   (t36–42) — far beyond tapered β-drift (0.6–1.0 measured) and beyond the
   known eff_y over-translation at the 256² grid (≤ +13% of background
   ≈ +0.7 m/s here).
3. Onset (t12–24) predates the re-intensification (44 → 84 m/s over
   t24–48), so intensity-coupled β-drift can amplify the late excess but
   cannot explain the onset.
4. Zonal channel behaved: u excess ≈ −0.4…−0.9 m/s (β-drift west minus the
   numerics east drift), and fixed-time cross growth ~0.9 m/s as above.
5. Retrospective: Katrina's −3.2 h early carried the same signature in
   miniature (+95 km along-track by t+36); Hugo's spin-down masked it.
   Ivan amplified it via the longest run, the strongest re-intensification
   — and/or the UNTESTED GRID (first 320²/5000 km production run; the
   harness eff_y numbers exist only for 256²/4000 and finer-dx variants,
   and the known issue says over-translation worsens with grid changes).

**Decisive next experiments (cheap, in order; GATES the R34/taper work
because eff contamination confounds any propagation experiment):**
1. Harness meridional + effx at the EXACT Ivan grid (nx=320, dom=5000,
   v_cap=70, Vmax≈64): measures pure-translation fidelity. eff_y ≈
   1.25–1.45 ⇒ grid-scaling artifact dominant (−8.3 h substantially
   explained; find the dx-independent scaling — spectral-filter indexing
   vs physical wavenumber is a suspect — or fall back to 4000/256 with
   the documented taper graze).
2. Harness betadrift at 320/5000, including a high-Vmax variant ⇒ β-drift
   and its intensity dependence on this grid.
3. Only if both are clean: the excess is storm-state physics and the
   question gets genuinely deep (vortex-flow interaction during
   re-intensification).


### Gate result (run_translation_test gate): GRID CLEAN
eff_y = +1.042 on BOTH 256²/4000 km and 320²/5000 km (Δ = 0.000, identical
Vmax'_end) — domain size and spectral mode count have no measurable effect
at fixed dx. Two corollaries: (1) the tapered over-translation is 1.042 vs
the historical untapered 1.13 — the wind taper CUT the over-translation,
retroactively explaining part of why it was vortex-size-dependent; (2)
grid-borne excess accounts for only ~+0.2 m/s of storm-Ivan's +1.4…+2.9 ⇒
the −8.3 h is storm-state physics. Next discriminator: the `ladder` mode
(β×steering → lockstep ramp → warm core, one rung at a time on the Ivan
grid), then gate-beta for the re-intensification coupling.


### Ladder result: the anchor moved
L0 (identical config to gate G1) reads eff_y = 1.230 at 12 h vs 1.042 at
10 h ⇒ ~78 km covered in hours 10–12 ≈ 10.8 m/s on a 5 m/s background — a
late-onset, accelerating excess in the MINIMAL configuration (f-plane,
constant background; no β, no ramp, no warm core). This reproduces the
storm-Ivan signature in the harness. L3 ≡ L2 bit-identically: θ′ is fully
passive — warm core exonerated. L1/L2 rung deltas are confounded by onset
timing and are unreadable until the anchor is understood. Analytical
lever: the L0 configuration is Galilean-removable, so genuine excess can
only be a frame-dependent discretization error (advection dispersion) —
or the vorticity-center is hopping late in the run. The `trace` mode
(L0 at 14 h with per-interval velocities + a zero-background stationary
control) discriminates; its 10 h checkpoint also replicates the gate on
the edited harness, closing the old-code/new-code loop.


### Trace verdict: lattice quantization, and a measurement-infrastructure fix
The Galilean control occupied exactly two grid positions for 14 h: the
centre-finder returns cell-snapped fixes with ±1.5-cell flicker (±23 km at
dx=15.6), an eff noise floor of ±0.09–0.13 at 10–14 h horizons — not the
±0.04 previously claimed. Re-reads: the ladder "anchor move" and the
apparent late acceleration were flicker (true eff_y ≈ 1.05 ± 0.09); all
ladder rung deltas are sub-noise; the 10 h checkpoint reproduced the gate
exactly, clearing the harness edits. Storm-Ivan's −8.3 h is cumulative
(~260 km ≫ flicker) and stands; its per-segment excess rates carry
±1–1.5 m/s. CONSEQUENCE FOR THE RESIDUAL BUDGET: the historical effx
asymmetry (W 0.87 / E 1.13) and betadrift (+0.43 E) — the entire evidence
base for the "~0.6 m/s spurious eastward drift" — were measured inside
this noise floor and must be re-measured; the drift may dissolve,
collapsing the two-storm +127/+131 km residual to pure propagation
deficit. FIX (V8.6.2, default-on): sub-cell centre-finding — ζ²-patch
centroid refinement in the harness, 3×3 parabolic refinement of the θ′
argmin in storm_tracker (unit-tested to <0.02-cell error on a synthetic
minimum). Re-run order: trace (the control track is now the validation —
expect ≲ ±2 km), then the ladder, then gate-drift.


### Trace re-run (V8.6.3 core-scale centroid): instrument validated, eff family retired
Galilean control: net drift 0.6 km over 14 h (0.01 m/s), positions within
±0.4 km — acceptance passed; the centre-finder noise floor improved ~50×.
Moving frame: **eff_y = 1.000 ± 0.007** with zero lateral drift (±0.8 km /
14 h). A balanced vortex in uniform flow translates EXACTLY faithfully:
the model is frame-invariant to <1%. Every historical eff>1 measurement
(1.13 untapered, 1.35 finer-dx, 1.042 gate, all ladder rungs) is hereby
attributed to lattice quantization of the old centre-finder and retired.
The taper's apparent reduction of over-translation (1.13→1.04) is likewise
retired as artifact-on-artifact. Presumptive (pending gate-drift): the W/E
asymmetry and the 0.6 m/s "spurious eastward drift" dissolve the same way,
collapsing the Hugo/Katrina +131/+127 km residual to pure propagation
deficit. Storm-scale results are unaffected (cumulative displacements far
exceed flicker); per-segment rates from pre-V8.6.3 logs carry ±1–1.5 m/s.
Minor open note: a −0.3 m/s vy decline in hours 11–14 of run A (eff(14h)
= 0.992), near the noise floor; logged, not load-bearing.
