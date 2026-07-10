# Over-rotation candidate hunt — design + registered predictions

**Date:** 2026-07-03 · **Status:** harness built + CPU-smoked; predictions registered BEFORE any
GPU run (this file must not be edited after the runs except to append results).

## Question

Paper §5.1 characterizes the over-rotation as "right amplitude, wrong orientation, and the
orientation never locks," and names three untested candidates for what fails to arrest the
gyre winding. This campaign tests two of them directly and the third by construction:

1. **The outer-wind cutoff** — the compact-support taper lays a negative-vorticity ring across
   the radii where the β-gyres live (gyre annulus ~150–450 km; ramp 200–500 km).
2. **The barotropic configuration** — a vortex with live thermodynamic structure disperses
   β-Rossby energy differently; production runs are barotropic by configuration (§2.1).
3. **The taper shape** (as distinct from its radius, which gate-beta-taper already swept to a
   floor of ~343°).

**Prior evidence in hand:** untapered Holland also drifted poleward (~2.17 m/s N, R_env-inert,
old betadrift mode) — a prior *against* the strong form of candidate 1, but confounded by the
domain-filling circulation. Aim is a known function of Vmax_end (structure probe: W ≈ const
0.41, N ∝ Vmax) — so **every read below is against the westward component first, heading
second, and always at reported Vmax_end** (the formulation-probe lesson: aim rotations
collinear with vortex collapse are confounds).

## Arm A — `gate-beta-shape` (candidates 1 + 3, one mode)

Fixed production geometry (onset 200 km, R_env 500 km, Vmax 64, cap 70, 320²/5000 km, 48 h,
β-plane, quiescent). Five profiles in a ring-sharpness ladder:

| row | profile | ring character |
|---|---|---|
| linear | ramp with slope kinks at 200 & 500 km | sharpest |
| cos | production control | smooth slope, curvature kinks |
| smooth5 | quintic smoothstep | smoothest compact ramp |
| gauss r_d=420 km | Gaussian envelope, **no cutoff** (V(350) matched to control ≈15.8) | no ring |
| gauss r_d=560 km | Gaussian envelope, broader | no ring; size control |

**Registered predictions (Claude, 2026-07-03):**
- P-A1: The three compact shapes (linear/cos/smooth5) read within **±3° heading, ±0.08 m/s
  west** of each other at comparable Vmax_end — i.e. ramp *form* does nothing. Confidence ~70%.
- P-A2: The gauss rows differ from control by ≲5° and the difference tracks the *size* change
  (gauss-560 more poleward than gauss-420, mirroring gate-beta-renv's broader→more-poleward),
  not the absence of the ring. Confidence ~60%.
- P-A3 (the surprise branch): if the ring is the mechanism, gauss rows rotate ≥10° toward NW
  with west rising ≥+0.2 m/s at preserved Vmax_end. I give this ~20%.

**Decision rule (as coded in the mode):** heading span <6° AND west span <0.15 m/s across all
five rows ⇒ **candidates 1+3 exonerated** at fixed size; anything larger ⇒ check monotonicity
in ring sharpness and non-collinearity with Vmax_end before claiming a mechanism.

## Arm C — `gate-beta-baroclinic` (candidate 2)

Same production geometry, cos taper. Ladder in baroclinicity persistence:

| row | θ′ | buoyancy | cooling τ | expected state at mature window (30–48 h) |
|---|---|---|---|---|
| dry control | zeroed | off | 30 min | the anchor (≡ published testbed) |
| passive null | kept | off | 30 min | θ′ passive ⇒ drift MUST equal control (harness integrity) |
| baroclinic 30 min | kept | on | 30 min | warm core dead by ~2.5 h ⇒ near-dry |
| baroclinic 6 h | kept | on | 6 h | baroclinicity persists into the window |
| baroclinic ∞ | kept | on | off | fully persistent; stability frontier |

CPU smoke confirmed the wiring (θ′ ≈ 45 K balanced core retained with buoyancy live, stable at
smoke scale); 48 h × 320² stability of the no-cooling row is itself a datum, not a failure —
the 6-h row still reads if it blows.

**Registered predictions (Claude, 2026-07-03):**
- P-C1: passive null ≡ dry control to tracker precision (else harness bug — stop). ~95%.
- P-C2: baroclinic-30min within 3° / 0.05 m/s west of control (core dead before the window). ~85%.
- P-C3 (the live question): west gain of the 6-h and no-cooling rows over control. If
  candidate 2 is the mechanism, west rises ≥+0.3 m/s toward the canonical ~1.4 and heading
  rotates NW monotonically with τ. I give "implicated" ~40%, "flat ⇒ exonerated" ~50%,
  "unreadable (instability/Vmax confound)" ~10%.

**Decision rule:** monotone west gain with τ at a live vortex ⇒ candidate 2 implicated (a
paper-grade finding: the barotropic reduction causes the characterized bias). Flat ⇒
exonerated within this core.

## If everything is flat

All three §5.1 candidates dead ⇒ the over-rotation is intrinsic to balanced barotropic β-gyre
dynamics in this core at this resolution — which sharpens §5.4's portability question ("do
other cores over-rotate?") into the primary next move, and makes the NU4×resolution interplay
(res sweep at reduced ν₄ at fine dx, where it's stable) the remaining in-house probe.

## Cost & run order (GPU, woe_env, repo root)

```powershell
python -m oracle_v8.run_translation_test gate-beta-shape        # ~5 × 15-20 min
python -m oracle_v8.run_translation_test gate-beta-baroclinic   # ~5 runs; buoyancy rows slower (fast stages live), ~2-3 h total
```

Arm A first (cheaper, and its outcome doesn't gate Arm C — they can also run back-to-back
unattended). Harness changes: `vortex_init.py` (taper_shape, outer_envelope_m — defaults
bit-identical, verified), `run_translation_test.py` (pass-throughs, buoyancy/tau_cool wiring,
max|w|/max|θ′| telemetry, the two modes, ORACLE_REQUIRE_GPU env override). All smoked on CPU.

## Results (appended 2026-07-03, GPU runs by Justin — predictions above untouched)

### Arm A — gate-beta-shape (mature 30–48 h)

| profile | \|drift\| | hdg | west | Vmax_end |
|---|---|---|---|---|
| linear 200→500 km | 2.19 | 348 | +0.47 | 42.3 |
| cos 200→500 km (control) | 2.49 | 350 | +0.42 | 42.2 |
| smooth5 200→500 km | 2.64 | 352 | +0.35 | 43.7 |
| **gauss r_d=420 km (no cutoff)** | **2.30** | **329** | **+1.20** | **39.7** |
| **gauss r_d=560 km (no cutoff)** | **2.85** | **326** | **+1.60** | **41.2** |

**Prediction scorecard:** P-A1 **CONFIRMED** (compact shapes span 4° / 0.12 m/s — ramp form does
nothing). P-A2 **WRONG** (registered ~60%): the gauss rows did not track size — they rotated
21–24° into the canonical NW band. P-A3, the registered ~20% branch, is what happened —
decisively, and cleaner than its own success criterion (≥10°, ≥+0.2 west; observed 21–24°,
+0.8–1.2 west).

**Verdict: candidate 1 (the compact-support cutoff) IMPLICATED; candidate 3 (ramp shape)
EXONERATED.** With a smooth Gaussian outer envelope the model produces β-drift with canonical
magnitude AND canonical aim: gauss-420 reads |2.30| m/s @ 329° with west +1.20 — inside both
theory bands — at preserved vortex (Vmax_end 39.7 vs control 42.2; not the Vmax-collapse
confound). Confound checks all pass: (i) size ruled out — the two families have OPPOSITE size
trends (compact broader→more poleward per gate-beta-renv; envelope broader→more westward); (ii)
mid-radius wind reduction ruled out — gauss-560 has more mid-radius wind than gauss-420 yet is
more canonical, and the onset-radius sweep never broke the 343° floor; (iii) domain/taper-zone
contamination ruled out — envelope wind ≈ 0.05 m/s at 1000 km. Physical reading: consistent
with Fiorino & Elsberry (1989) — β-drift is controlled by the outer wind (300–1000 km); the
compact taper amputates it, removing the ambient flow that phase-locks the gyre pair. The old
"untapered Holland also poleward" prior is superseded (pre-subcell measurement inside the noise
floor on a domain-filling profile).

### Arm C — gate-beta-baroclinic (mature 30–48 h)

| config | \|drift\| | hdg | west | Vmax_end | max\|w\| | θ′_max |
|---|---|---|---|---|---|---|
| dry control | 2.49 | 350 | +0.42 | 42.2 | 0.04 | 0.2 |
| passive null | 2.49 | 350 | +0.42 | 42.2 | 0.04 | 0.2 |
| baroclinic τ=30 min | 2.51 | 350 | +0.42 | 43.0 | 0.04 | 0.2 |
| baroclinic τ=6 h | 2.44 | 354 | +0.24 | 73.3 (cap-pinned) | 0.91 | 95.5 |
| baroclinic no cooling | — | — | — | BLEW UP | — | — |

**Prediction scorecard:** P-C1 **CONFIRMED** (null ≡ control, Δ = 0.0°). P-C2 **CONFIRMED**
(τ=30 min ≈ control). P-C3: **UNREADABLE as designed** (the registered ~10% branch): relaxing
θ′ toward zero cannot hold a *persistent balanced* warm core — at τ=6 h the adiabatic
generation ran to θ′ = 95 K (3× beyond LH82's validated ~10% θ′/θ̄ regime) and the cap pinned
the vortex; the no-cooling row is the expected stability datum. Candidate 2 is neither
implicated nor exonerated. **Redesign (Arm C-v2):** Newtonian relaxation toward the INITIAL
balanced warm-core θ′ field (θ′_ref) instead of toward zero — persistent baroclinicity at
bounded amplitude with no runaway pathway. One small component change; run only if mechanistic
completeness is wanted after the Arm A follow-ups.

### Phase-lock confirmation (appended 2026-07-03, from the Arm A log's per-window traces)

```
gauss r_d=420:  t12: 1.11@328   t24: 1.73@327   t36: 2.20@328   t48: 2.37@329
gauss r_d=560:  t12: 1.50@329   t24: 2.25@325   t36: 2.83@325   t48: 2.81@326
```

Heading is FLAT from t12 to t48 (±2°, tracker precision) while amplitude grows to saturation —
vs the compact family's steady ~0.4°/h precession through north (gate-beta-timeevol/longrun).
The orientation locks *early* (by t12, at half the final amplitude) and the amplitude grows into
a fixed phase: the canonical equilibration, complete. The "over-rotation" was never a rotation
rate error — it was the absence of the lock, and the lock is supplied by the outer vortex flow
the compact taper amputates. Mechanism closed at the testbed level.

### Consequences & next steps

1. **The characterized bias is a cost of the vortex-boundedness choice, and it appears
   removable.** The single structural parameter changes meaning: envelope scale r_d instead of
   taper-onset radius. gauss-420 is already in-band on both magnitude and aim on first guess.
2. **Immediate (free):** read the per-window heading traces (t12/24/36/48) for the gauss rows
   from the existing Arm A log — if the gauss heading HOLDS ~326–330 from t24 onward instead of
   marching through north, the gyre is genuinely phase-locked, which is the mechanistic
   confirmation (amplitude AND phase equilibrated).
3. **Calibration + invariance checks (one GPU session):** r_d sweep (350/420/500/560) for the
   band sweet spot; Vmax 64/35/21 at the chosen r_d (canonical structure should keep west
   ≈ const and now heading ~NW at all intensities); optional f-plane null under gauss.
4. **The big one — six-storm re-run with the envelope profile,** registered predictions first.
   The paper's own attribution (steering controls cross-track; bias lives in along-track on
   poleward movers) makes falsifiable predictions: along-track/timing on poleward movers
   (Michael, Ivan, Katrina) should improve; cross-track should be ~unchanged. Either outcome is
   valuable: improvement quantifies the bias's real-track footprint; no change double-confirms
   subdominance.
5. **Paper decision (Justin's):** §5.1's candidate list is now partially resolved by these
   runs. Options: one-sentence update to §5.1 (calm: "a subsequent profile-family test
   implicated the compact-support cutoff; the recovered drift is canonical in the testbed"),
   or hold the whole finding for the follow-on paper. The §4.1 characterization is untouched
   either way (profile *family* was never one of its three swept levers).

## Stage 2 — `gate-beta-envelope` calibration (registered 2026-07-03, BEFORE first run)

Seven runs: r_d sweep 350/420/500/560 km at Vmax 64; intensity ladder Vmax 64/35/21 at
r_d = 420; f-plane null at r_d = 420. Usage:
`python -m oracle_v8.run_translation_test gate-beta-envelope` (~7 × 15–20 min).

**Registered predictions (Claude, 2026-07-03):**
- P-E1: |drift| increases with r_d (size-controlled magnitude, both families agree on this);
  heading stays in the NW band across ALL r_d — the phase-lock is not r_d-fine-tuned.
  In-band-on-both sweet spot around r_d ≈ 380–450. Confidence ~75%.
- P-E2: at r_d = 420, heading stays ~NW (within ±8°) across Vmax 64/35/21 while speed scales
  with intensity — canonical structure (Fiorino–Elsberry: outer wind controls, direction
  ~intensity-independent). The compact family read 338→351 with Vmax; if the envelope family
  reproduces that poleward rotation, the lock is amplitude-dependent and the mechanism is
  partial. Confidence heading-invariant ~70%.
- P-E3: f-plane null under the envelope: |drift| ≲ 0.05 m/s (the envelope must not manufacture
  drift). Confidence ~90%.

**Decision rule:** all three pass ⇒ pick the in-band r_d and proceed to production wiring +
the six-storm re-run (its own registered predictions to be written before that run). P-E2
fails ⇒ characterize aim = f(Vmax) within the envelope family before any production change.

### Stage 2 results (appended 2026-07-03, GPU runs by Justin)

| config | \|drift\| | hdg | west | Vmax_end |
|---|---|---|---|---|
| r_d=350 km, Vmax 64 | 1.99 | 330 | +0.99 | 38.7 |
| r_d=420 km, Vmax 64 | 2.30 | 329 | +1.20 | 39.7 |
| r_d=500 km, Vmax 64 | 2.62 | 327 | +1.43 | 40.5 |
| r_d=560 km, Vmax 64 | 2.85 | 326 | +1.60 | 41.2 |
| r_d=420 km, Vmax 35 | 1.63 | 326 | +0.92 | 25.7 |
| r_d=420 km, Vmax 21 | 1.23 | 324 | +0.72 | 17.2 |
| r_d=420 km, Vmax 64, f-plane | 0.00 | — | −0.00 | 37.4 |

**Scorecard: P-E1 CONFIRMED** (heading span 4° across r_d 350–560; magnitude size-controlled —
aim and size have decoupled, unlike the compact family where they were entangled to the 343°
floor). **P-E2 CONFIRMED** (heading span 4° across Vmax 64/35/21 vs the compact family's 13°
poleward march; note the *kind* of invariance changed — compact was fixed-west/scaling-north
i.e. direction rotating with intensity; envelope is fixed-direction/scaling-magnitude, west
0.72→1.20 in proportion — the canonical β-drift structure, recovered whole). **P-E3 CONFIRMED**
(f-plane null 0.002 m/s; the envelope manufactures nothing). Footnote, not confound: envelope
equilibrium vortices run slightly weaker (Vmax_end 37–41 vs 42 — envelope trims mid-radius
wind); P-E2 itself shows aim is Vmax_end-independent.

**Mechanism + fix now validated at testbed level: r_d-robust, intensity-invariant,
artifact-free phase-lock.**

### The r_d freeze (decision required BEFORE any production run)

Two candidates pass both bands. The choice must be frozen before the six-storm re-run and
justified by testbed physics only (the no-landfall-tuning discipline):

- **r_d = 350 km** — the mode's auto-pick (band-center magnitude |1.99|); hdg 330, west +0.99.
- **r_d = 420 km** — recommended (Claude): (i) V(350 km) is matched to the production taper
  control *by construction*, so the six-storm A/B is maximally "same vortex, plus the tail" —
  the cleanest attribution; (ii) drift magnitude 2.30 is closest to the compact control's 2.49,
  minimizing the change attributable to overall drift size rather than aim; (iii) west +1.20
  sits nearer the canonical ~1.4.

Production wiring is in place and inert: `production_config.OUTER_ENVELOPE_M = None` (None ⇒
the published cosine-taper profile, bit-identical). Freezing the value + writing the per-storm
registered predictions is Stage 3's opening move.

### Stage 3 — six-storm re-run: REGISTERED PREDICTIONS (2026-07-03, r_d FROZEN at 420 km,
### written BEFORE any storm run — do not edit; append results below)

**Configuration:** `production_config.OUTER_ENVELOPE_M = 420_000.0` (frozen; chosen from
Stage 2 on testbed grounds only — V(350 km) matched to the taper control, |drift| nearest
control). Treatment = all six storms via `run_storm`, identical everything else; control = the
checked-in `*_Agnostic` logs. Score with `landfall_verify` (both metrics, all six).

**The treatment's drift delta (testbed, mature, Vmax 64):** control (2.49 m/s @ 350°: west
0.42, north 2.45) → treatment (2.30 @ 329°: west 1.20, north 1.96), i.e.
**Δ = (−0.78 east, −0.49 north) m/s.** Projecting Δ through each storm's landfall geometry
(same headings/transits as the projection test — the *differential* form of that calculation,
which is better-posed: steering errors cancel in the A/B):

| storm | Δcross (km) | pred cross (obs + Δ) | Δalong (km) | pred along |
|---|---|---|---|---|
| Hugo | −93 | +110 → **+17** | +5 | +23 → +28 |
| Katrina | −108 | +125 → **+17** | −44 | +77 → +33 |
| Ivan | −139 | +126 → **−13** | −30 | +249 → +219 |
| Fran | −81 | +8 → **−73** | −10 | −46 → −56 |
| Michael | −73 | −99 → **−172** | −66 | +124 → +58 |
| Laura | −74 | −32 → **−106** | −30 | +37 → +7 |

**The honest headline of this table: the strong-form prediction is that the fix TRADES the
discovery set's eastward errors for test-set westward errors** (six-storm cross RMS 95 → ~87,
approximately flat) — which is precisely what the paper's steering-dominance verdict implies.
A physics-correct fix to a subdominant term should not buy much cross-track skill. The robust
gains are where the bias demonstrably lived: **along-track on the poleward movers.**

**Registered predictions (Claude, 2026-07-03; bands ±50% on all Δ magnitudes — linear
accumulation and intensity-scaling of Δ are both approximate):**
- P-S1 (robust, hypothesis-independent): along-track improves on Katrina, Michael, and Laura;
  Michael's early-arrival timing shrinks by ≥1.5 h. Confidence ~70%.
- P-S2 (strong form, "bias-additive" H1): per-storm cross-track shifts west by 60–150 km
  (centrals above); discovery storms land +0 ± 60; Fran and Laura go clearly west. ~45%.
- P-S3 ("feedback-compensated" H2): the lockstep steering relaxation partially absorbs
  self-propagation changes (the storm samples the DLM where it actually is), so observed
  |Δcross| come in at LESS THAN HALF the strong-form centrals, roughly uniformly. ~35%.
  (H2 would itself be a finding: it quantifies how much the steering architecture buffers
  self-propagation error — directly relevant to the paper's attribution logic.)
- P-S4 (guards): no storm's Vmax history changes by more than ~10 m/s (envelope trims
  mid-radius wind; testbed Vmax_end dropped ~2.5); ERA5/steering path identical; any timing
  change beyond ±3 h vs control on the direct storms (Hugo, Katrina) is a red flag, not a
  result. Remaining ~20% mass: something outside H1/H2 (e.g., envelope-altered decay
  interacting with steering sampling) — decompose per-segment before interpreting.

**Scoring rules:** (i) score P-S1 first — it is the claim the whole campaign stands on;
(ii) classify H1 vs H2 by the ratio of observed to strong-form Δcross; (iii) report six-storm
cross RMS and the same-latitude metric alongside landfall-fix (Ivan/Michael recurve geometry
inflates landfall-fix); (iv) either H1 or H2 confirms the bias's real-track footprint is
bounded and now characterized from BOTH sides — the follow-on paper's empirical spine.

### Stage 3 results (appended 2026-07-03; six treatment runs by Justin, logs
### `Logs/<Storm>/<storm>_gauss_envelope.txt`; scored against the registrations above)

Landfall-fix decomposition, control → treatment (Δ observed vs strong-form central):

| storm | cross c→t | Δcross (pred) | ratio | along c→t | Δalong (pred) | timing c→t | ΔVmax_end |
|---|---|---|---|---|---|---|---|
| Hugo | +110.2 → +78.7 | −31.5 (−93) | 0.34 | +23.3 → +39.7 | +16.4 (+5) | −2.3 → −2.1 | **−14.2 ⚠** |
| Katrina | +124.6 → +78.2 | −46.4 (−108) | 0.43 | +76.5 → +61.4 | −15.1 (−44) | −2.6 → −2.4 | −2.1 |
| Ivan | +126.3 → +40.2 | −86.1 (−139) | 0.62 | +249.3 → +210.1 | −39.2 (−30) | −8.1 → −7.1 | **−24.6 ⚠** |
| Fran | +7.7 → −29.9 | −37.6 (−81) | 0.46 | −45.5 → −50.1 | −4.6 (−10) | +1.2 → +1.9 | −2.1 |
| Michael | −98.7 → −113.2 | −14.5 (−73) | 0.20 | +123.8 → +93.4 | −30.4 (−66) | −5.2 → −5.1 | −1.7 |
| Laura | −31.9 → −51.8 | −19.9 (−74) | 0.27 | +37.4 → +30.0 | −7.4 (−29) | −1.3 → −1.8 | −1.5 |

Same-latitude: Hugo +102.7→+74.8, Katrina +114.3→+78.2, Ivan +67.6→**+5.7**, Fran +23→−20.7,
Michael −75→−97.0, Laura −39→−42.0. Aggregates: landfall-fix cross RMS **95.3 → 71.1 (−25%)**;
clean-four (guard-passing storms) cross RMS 81.2 → 75.0 (−8%); same-lat RMS 77.3 → 62.4; along
RMS 120.7 → 101.3.

**Prediction scorecard:**
- **P-S1: PASS on its core clause** — along-track improved on Katrina, Michael, AND Laura (all
  guard-clean). The Michael-timing sub-clause FAILED (−5.2 → −5.1 h; the along improvement at
  the landfall fix did not move the threshold-crossing clock).
- **P-S2 (H1, strong form): FAILED** as registered — observed shifts are far below the
  centrals; the discovery storms did not land at +0 ± 60 (Hugo +79, Katrina +78).
- **P-S3 (H2): CONFIRMED** — the sign was right SIX FOR SIX (every storm shifted west), at a
  mean transmission ratio ≈ 0.34 on the guard-clean four (0.20–0.46; Ivan's 0.62 excluded as
  confounded). **The lockstep steering relaxation absorbs roughly two-thirds of a
  self-propagation change over a landfall transit** — a built-in, physical, now-measured
  buffering pathway, and the quantitative reason the aim bias was subdominant at landfall.
- **P-S4: the Vmax guard CAUGHT REAL CONFOUNDS** — Hugo (−14.2) and Ivan (−24.6) ran
  substantially weaker under the envelope (both are the storms with strong barotropic
  re-intensification phases; the envelope trims the mid-radius angular-momentum reservoir that
  feeds dry spin-up). Their large improvements (Ivan same-lat +67.6 → +5.7!) are therefore
  entangled with reduced drift magnitude, not pure re-aim, and are flagged, not headlined.
  Timing guard passed everywhere (direct storms Δ ≤ 0.2 h).

**Honest bottom line:** the testbed fix propagates to real storms with the predicted sign in
all six cases; cross-track skill improves modestly on the clean storms (−8% RMS) and more in
aggregate (−25%) only via the two intensity-confounded storms; along-track improves broadly.
The paper's steering-dominance claim is doubly confirmed — even *removing* the bias buys only
modest cross-track change through the steering buffer. New flag for the roadmap: **the envelope
profile damps barotropic re-intensification** (Hugo/Ivan) — mechanism plausibly the trimmed
mid-radius wind reservoir; needs its own study if intensity ever becomes a quantity.

**Campaign status: COMPLETE at the testbed-to-storms level.** Mechanism (compact-support
cutoff → no gyre phase-lock), fix (Gaussian envelope, r_d = 420 km frozen), calibration
(r_d-robust, intensity-invariant, artifact-free), and real-track validation (sign 6/6,
transmission ratio ~1/3, guards caught the confounds) — every stage under registered
predictions. Follow-on paper skeleton is all here.

## Arm C-v2 — relax-to-θ′_ref (registered 2026-07-07, BEFORE any GPU run)

The corrected baroclinicity test paper-2 §2.3 awaits: `NewtonianCoolingComponent` gains an
optional `theta_ref`; `run_translation(cool_to_init=True)` holds θ′ at the post-prebalance
balanced core, giving persistent BOUNDED baroclinicity with departures damped at τ — the design
Arm C-v1's runaway demanded. Mode: `gate-beta-baroclinic-v2` (CLI 34; 4 rows × 48 h ≈ 1 h GPU).
Context: with the cutoff mechanism established, this is a completeness check — "does live
vertical structure ALSO move the aim," not "which candidate explains the bias."

**Registered predictions (Claude, 2026-07-07):**
- **P-C2v1 (integrity, ~90%):** passive null v2 (held core, buoyancy OFF) reproduces the dry
  control's drift (Δheading < 1°, Δwest < 0.03) — θ′ is passive without buoyancy, whatever the
  cooling target.
- **P-C2v2 (boundedness, ~75%):** the held baroclinic rows stay bounded — max θ′ within ~2× the
  balanced core (≈42 K), max|w| < 5 m/s, no cap-pinning — the runaway pathway is closed by
  construction (departure equilibrium ≈ w·dθ̄/dz·τ ≈ 2 K at τ=30 min).
- **P-C2v3 (the verdict, ~70% exonerate):** at bounded θ′, the compact-taper aim moves < 5° and
  west changes < 0.15 m s⁻¹ vs the dry control ⇒ **candidate 2 EXONERATED** — the cutoff owns
  the whole bias (the envelope already recovered canonical aim with zero baroclinicity, leaving
  no residual to explain). Larger movement toward NW ⇒ vertical structure contributes
  independently ⇒ paper-2 §2.3 and §6 get the richer sentence.
- **Decision:** either branch closes paper-2's §2.3 hook; the τ=6h row checks the verdict is
  not an artifact of hard anchoring.

### Arm C-v2 results

*(append after the GPU run; predictions frozen)*

### Arm C-v2 results (appended 2026-07-07, GPU run by Justin)

| config | |drift| | hdg | west | Vmax_end | max|w| | θ′_max |
|---|---|---|---|---|---|---|
| dry control | 2.49 | 350 | +0.42 | 42.2 | 0.04 | 0.2 |
| passive null v2 (held, buoy OFF) | 2.49 | 350 | +0.42 | 42.2 | 0.04 | 44.7 |
| baroclinic HELD τ=30min | 1.76 | 356 | +0.13 | **79.9 (cap)** | 0.25 | 44.8 |
| baroclinic HELD τ=6h | BLEW UP | | | | | |

**Scorecard:** P-C2v1 **PASS exactly** (null ≡ control at 0.0° carrying the held core
passively). P-C2v2 **PARTIAL** — the runaway is closed as designed (θ′ pinned at 44.8 K, w
0.25) but the no-cap-pinning clause failed for a physical reason: **a maintained warm core is
an energy source** — the relax-to-ref term restores the core against every erosion, a
continuous input, and the vortex intensifies to the cap (42 → 80). The τ=6h blow-up brackets
the other side: the anchoring window between "driven" and "runaway" is narrow. P-C2v3: the
registered thresholds fired the "contributes" branch (Δwest −0.30, rot +5.7°), but the
comparison is intensity-entangled (79.9 vs 42.2) — **and decidable by SIGN regardless: live
maintained baroclinicity moved the aim POLEWARD (west 0.42 → 0.13), the wrong direction to
explain the westward deficit. Candidate 2 is EXONERATED AS CAUSE of the bias**; only the
magnitude of its wrong-way nudge remains intensity-entangled.

**Confidence audit:** P-C2v1 90% ✓; P-C2v2 75% partial; P-C2v3's registered branches did not
include the decisive one (sign-resolved) — lesson, again: register sign-resolved branches for
directional quantities.

**P-C2v4 (matched-intensity cleanup — registered 2026-07-07, optional, ZERO new code):** the
existing `gate-beta` mode's second row (init 120 → cap 70) is a dry cap-pinned control. If its
mature aim reads ≥354° with west ≤0.25, the baroclinic row's shift is intensity alone (the
nudge evaporates); if it holds ~350–352/≈+0.4, the small poleward baroclinic nudge is real.
~60% on the first branch. Either way the paper-level verdict above stands.
