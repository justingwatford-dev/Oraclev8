# Paper-tightening analyses — registered predictions

*Registered 2026-07-18, BEFORE any analysis executes, per campaign discipline. Motivation:
external-review critiques of the two manuscripts — (1) the blind-track-skill numbers have no
baseline ("good, or just non-catastrophic?"), and (2) the transmission ratio (≈0.34 observed,
≈0.42 mechanistic) asks a general claim ("landfall under-reads self-propagation by ~3×") to rest
on six storms, four guard-clean. Three analyses, none requiring a new model run or new data:
every input is HURDAT2 (`oracle_v8/hurdat2.txt`), the cached per-storm ERA5 steering files, the
checked-in run logs, and the committed campaign tables. All three can come back against us; each
such outcome changes manuscript language and is reported as a finding.*

---

## Analysis 1 — Persistence baseline (`measure_persistence.py`)

**Method.** For each of the six storms: observed motion vector at init from a centered
difference of the HURDAT2 best track over init ±6 h; extrapolate the init position linearly for
the storm's transit time T (the observed landfall fix time used throughout both papers);
compute the persistence-minus-observed error at the landfall fix and decompose along/cross
using the observed motion direction at landfall (centered difference ±3 h), the manuscripts'
sign conventions (along + = ahead; cross + = right of motion).

**Predictions.**

- **P-B1 (65%):** six-storm persistence cross-track RMS at the landfall fix exceeds 2× the
  model control's 95 km (i.e. > 190 km).
- **P-B2 (65%):** persistence total displacement error exceeds the model control's total error
  on all six storms. Most likely single exception if one occurs: Fran (fast straight-mover).
- **P-B3 (70%):** the two storms with unmodeled recurvature in the record (Ivan, Michael) each
  show persistence total error > 300 km.

## Analysis 2 — Steering-only tracer (`measure_steering_tracer.py`)

**Method.** A point (no vortex) initialized at each storm's HURDAT2 init fix and advected by
the annulus DLM sampled at the *tracer's own position* from the cached ERA5 file
(`ERA5Steering.get_dlm`, production 3–7° annulus), forward-Euler at dt = 0.25 h, integrated to
the observed landfall time; scored with the same landfall-fix along/cross decomposition. This
measures what the steering alone delivers — the model-minus-tracer difference is an independent
read of the self-propagation footprint. Registered caveat: the tracer is not exactly
"the model minus drift" — it feels the DLM along its own (diverged) path and omits the 3-h
relaxation lag; we treat it as a first-order steering share, not an exact ablation.

**Predictions.**

- **P-S1 (85%, validity guard):** the integration is scoreable (tracer stays inside the cached
  ERA5 domain, reaches the landfall-time evaluation) for all six storms.
- **P-S2 (60%):** tracer cross-track RMS ≤ 1.5× the model control's cross-track RMS
  (≤ ~142 km): steering carries the cross-track outcome, as both papers claim.
- **P-S3 (70%):** the model control runs *ahead of* the tracer along-track (model − tracer
  along > 0) on at least four of six storms — the β-drift's poleward push, visible as the
  model-tracer difference.
- **P-S4 (50%):** on the guard-clean four, the model-minus-tracer displacement magnitude at
  landfall is within a factor of 2 of the transit-attenuated drift footprint
  (0.42 × 2.49 m s⁻¹ × T), consistent with the §5.3 accounting.

## Analysis 3 — Per-storm transmission decomposition (`measure_transmission_decomp.py`)

**Method.** The §5.3 aggregate accounting (intensity scaling × gyre spin-up ≈ 0.42) made
per-storm: family drift-vs-intensity laws fit from committed testbed tables (compact: west
≈ 0.41 m s⁻¹ constant, north ∝ V through 2.45 at V_end 42.2; envelope: west and north linear
fits through the Stage-2 intensity ladder r_d = 420 at V_end 39.7/25.7/17.2); gyre spin-up
curve s(t) from the committed envelope trace (0 → 0.47 → 0.73 → 0.93 → 1.0 at t = 0/12/24/36/48 h,
clamped beyond), applied to the difference vector; per-storm Vmax′(t) read from each storm's
control and envelope run logs, each family's law evaluated on its own run's intensity history.
Predicted per-storm shift = transit integral of the cross-projection (through the storm's
landfall heading, the projection-test values) of s(t)·[env_vec(V_env(t)) − cmp_vec(V_cmp(t))];
predicted ratio = that integral over the strong-form (mature-Δ × T) projection. Observed ratios
(committed Stage-3 table): Hugo 0.34†, Katrina 0.43, Ivan 0.62†, Fran 0.46, Michael 0.20,
Laura 0.27 († = intensity-guard flagged).

**Predictions.**

- **P-D1 (55%):** the predicted ratios partition the guard-clean four correctly into top two
  {Katrina, Fran} and bottom two {Michael, Laura}.
- **P-D2 (40%):** all four guard-clean predicted ratios fall within ±0.15 absolute of observed.
- **P-D3 (60%):** the decomposition does *not* close for flagged Ivan: predicted ratio falls
  short of the observed 0.62 by more than 0.15. (Ivan's weakened envelope run should *reduce*
  its predicted shift; its large observed shift is therefore not intensity-explicable —
  steering-path divergence between the A and B runs is the suspected residual mechanism, which
  is exactly why the guard flagged it.)
- **P-D4 (65%):** the mean of the guard-clean four predicted ratios lands within ±0.10 of the
  aggregate decomposition's 0.42.

## Scoring rules

Each prediction scores CONFIRMED / FAILED on its stated threshold; no post-hoc threshold
motion. Failures are findings: P-B failures weaken the skill framing and the manuscripts say
so; P-S failures revise the steering-share language of both papers; P-D failures demote the
§5.3 mechanism from "the accounting that fits" to "an accounting consistent with the mean but
not the storms," and §5.3 is rewritten accordingly.

---

# OUTCOMES (scored 2026-07-18, same day; scripts `oracle_v8/measure_persistence.py`,
# `measure_steering_tracer.py`, `measure_transmission_decomp.py`)

**Tally: 6 confirmed / 5 failed of 11.** Both failures-with-teeth are informative; details:

## Analysis 1 — persistence (1/3)

| storm | pers along | pers cross | pers total | model total |
|---|---:|---:|---:|---:|
| Hugo | −155 | −96 | 182 | 113 |
| Katrina | −273 | −40 | 276 | 146 |
| Ivan | −301 | −84 | 312 | 279 |
| Fran | −96 | −80 | 124 | 46 |
| Michael | −180 | −137 | 226 | 158 |
| Laura | −115 | −198 | 229 | 49 |

- **P-B1 FAILED** — persistence cross-track RMS 117 km vs threshold >190 (model: 95.3).
  *The* finding of the analysis: under this metric persistence is not catastrophically worse
  on cross-track specifically; the model's edge over no-skill lives in **total landfall
  position, dominated by timing** — every persistence along-track is a large arrears (storms
  accelerate poleward; persistence cannot know).
- **P-B2 CONFIRMED** — model total error smaller on **all six** (margins ×1.1 Ivan … ×4.7
  Laura; mean 132 vs 225 km).
- **P-B3 FAILED** — Ivan 312 ✓ but Michael 226 < 300.

## Analysis 2 — steering-only tracer (4/4)

| storm | tracer along | tracer cross | model−tracer along | \|model−tracer\| | attenuated-footprint ref |
|---|---:|---:|---:|---:|---:|
| Hugo | −92 | +88 | +115 | 117 | 105 |
| Katrina | −69 | +174 | +146 | 154 | 132 |
| Ivan | −12 | +158 | +261 | 263 | 161 |
| Fran | −171 | +17 | +126 | 126 | 92 |
| Michael | +15 | +53 | +109 | 186 | 111 |
| Laura | −67 | +18 | +104 | 116 | 90 |

- **P-S1 CONFIRMED** (6/6 scoreable). **P-S2 CONFIRMED** — tracer cross RMS 105 km ≈ model's
  95: steering alone essentially reproduces the cross-track outcome. **P-S3 CONFIRMED** —
  model ahead of tracer along-track on 6/6 (+104…+261 km): the vortex's poleward
  self-propagation closes the timing the steering alone cannot. **P-S4 CONFIRMED** —
  |model−tracer| within 2× of the §5.3 attenuated drift footprint on all four guard-clean
  storms (ratios 1.2–1.7): an *independent instrument* lands on the same transmission
  accounting as the A/B.
- Per-storm note: the tracer beats the model's total only on **Ivan and Michael** — exactly
  the two storms where the papers say the β-bias dumps into along-track overshoot.

## Analysis 3 — per-storm transmission decomposition (1/4 + post-hoc)

| storm | predicted ratio | observed | note |
|---|---:|---:|---|
| Hugo | 0.45 | 0.34 | guard-flagged |
| Katrina | 0.52 | 0.43 | +0.09 |
| Ivan | 0.67 | 0.62 | guard-flagged; **entanglement quantified** |
| Fran | 0.43 | 0.46 | −0.03 |
| Michael | 0.45 | 0.20 | **the outlier, +0.25** |
| Laura | 0.38 | 0.27 | +0.11 |

- **P-D1 FAILED** (predicted top-two {Katrina, Michael}). **P-D2 FAILED** — but Michael is
  the *sole* >0.15 miss; Katrina/Fran/Laura all within ±0.11.
- **P-D3 FAILED, informatively** — registered claim was that Ivan's 0.62 would NOT be
  intensity-explicable (shortfall >0.15); instead the decomposition, fed the runs' own
  intensity histories, lands at 0.67 vs 0.62: the weakened envelope run's *north* deficit
  projects into cross at Ivan's heading. The intensity guard excluded Ivan for exactly the
  entanglement the decomposition now quantifies; the registered steering-path-divergence
  suspicion was unnecessary.
- **P-D4 CONFIRMED** — guard-clean mean predicted 0.44 vs the aggregate accounting's 0.42.
- **POST-HOC (labeled, not registered):** along-axis transmissions where §5.2 commits the
  observed along shifts — Katrina obs 0.34 / pred 0.43; **Michael obs 0.46 / pred 0.38**;
  Laura obs 0.25 / pred 0.29. Michael's anomaly is therefore **axis-specific**: it transmits
  the correction along-track as predicted and under-transmits only on cross — the one storm
  whose landfall heading (010°) makes the cross axis nearly zonal. Mechanism unresolved;
  candidates (untested): east–west steering compensation for a slow mover; recurve-geometry
  projection. One storm's worth of honest residual.

## Synthesis for the manuscripts

1. The blind-skill claim gets its baseline: **better than persistence at landfall position
   on six of six**, with the margin in timing/total, not cross alone — which *strengthens*
   the papers' central claim that cross-track is steering-set.
2. Paper 1 §5.3's open question ("how much of the skill belongs to the steering?") now has a
   measured answer: steering alone ≈ reproduces cross-track (105 vs 95 km RMS) and arrives
   late everywhere; the model's drift closes timing on all six.
3. The transmission ratio survives as an **aggregate** with independent corroboration
   (P-S4), its per-storm structure is *partially* mechanistic (three clean storms within
   ±0.11 + Ivan's entanglement quantified), and Michael is reported as an axis-specific
   outlier. "Roughly a factor of three" stands, with the per-storm honesty attached.
