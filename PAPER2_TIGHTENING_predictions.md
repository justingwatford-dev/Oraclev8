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
