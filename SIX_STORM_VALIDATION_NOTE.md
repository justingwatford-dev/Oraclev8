# Six-storm validation — what it does to the track-error claim

*Briefing note for external review (Woodruff + ensemble). Status: the out-of-sample runs are
complete and the negative result is solid; the positive reframing is open and is the reason this is
circulating. Nothing in the manuscript has been changed. This note is written to be argued with.*

---

## Bottom line

We added three storms the storm-agnostic configuration had **never been fit to** (Fran 1996,
Michael 2018, Laura 2020), specifically to test whether the paper's eastward cross-track "cluster"
generalizes beyond the storms used to set the one global tuning parameter.

It does not. The large eastward cross-track is a property of the **calibration set only**. All three
out-of-sample storms came in small, and two of the three **flipped sign**. The paper's empirical
centerpiece — "a consistent eastward cross-track displacement (~+120 km) at every simulated
landfall" — does not survive out-of-sample validation.

Two things are *not* damaged by this, and one of them is arguably a better headline than the claim
we lost:
1. The idealized **testbed β-gyre over-rotation** is untouched (it is a quiescent-environment
   characterization of the model's β-drift vector; it does not depend on any storm).
2. **Out-of-sample track accuracy is good** — two of three unseen storms landed within ~40 km
   cross-track; the third's error is mostly timing, not position.

---

## The six storms (same-latitude cross-track, the clean metric)

| Storm   | Set         | Cross-track | Timing  | Track character                         |
| ------- | ----------- | ----------- | ------- | --------------------------------------- |
| Hugo    | in-sample   | **+103 E**  | −2.3 h  | NNW→NW, SC coast                        |
| Katrina | in-sample   | **+114 E**  | −2.6 h  | N, Gulf → LA                            |
| Ivan    | in-sample   | **+68 E**   | −8.1 h  | NNW recurver (missed recurve)           |
| Fran    | out-of-sample | **+23 E**  | +1.2 h  | NNW, fast (~15 kt)                      |
| Michael | out-of-sample | **−75 W**  | −5.2 h  | near-pure-N, slow (~11 kt); obs recurves NE at landfall |
| Laura   | out-of-sample | **−39 W**  | −1.3 h  | WNW→NW→N (gentle recurve, ~13–15 kt)    |

In-sample: all three east, mean ≈ **+95 km**, spanning 46 km.
Out-of-sample: mixed sign, mean ≈ **−30 km**, spanning 98 km.

The three storms used to calibrate the taper parameter are the *only* three with a large eastward
cross-track. That is the whole result in one sentence.

**Linchpin to confirm (Justin):** this entire reading rests on Hugo / Katrina / Ivan being the
taper-calibration ("in-sample") set — the `taper_start_frac = 0.40` comment calls 0.40 the
calibrated landfall value. If the calibration history is different from what we've assumed, the
in/out-of-sample framing has to be redrawn. Please confirm before this goes further.

---

## What survives, what breaks

**Breaks:**
- "Consistent eastward cross-track cluster (~+120 km, spanning 16 km)" as a *general* property. Six
  storms read +103 / +114 / +68 / +23 / −75 / −39. That is neither a cluster nor one-signed.
- The **bridge claim**: "the idealized β-gyre bias manifests as a systematic eastward cross-track
  displacement at landfall." This specific attribution is exactly what the out-of-sample storms
  falsified.

**Survives, unaffected:**
- The testbed β-gyre over-rotation (f-plane control, intensity-independence, grid-convergence,
  ν₄/taper/resolution levers). It characterizes the model's self-propagation *vector* in a quiescent
  environment and is independent of any storm.
- Six-storm tracking skill, in aggregate. The simulator puts real storms near their landfalls.
- Out-of-sample accuracy specifically (Fran +23, Laura −39 cross-track; Michael position dominated
  by timing).

---

## How to read the per-storm errors — and where this gets dangerous

The tempting story: there is **one** characterized too-poleward β-drift, and it simply *manifests
differently* depending on storm speed and track curvature —
- moderate-speed storms whose observed track curves NW → eastward cross-track (Hugo, Katrina);
- a fast storm → small eastward cross-track, the drift diluted by strong steering (Fran);
- a slow storm whose observed track recurves NE → the excess poleward drift dumps into **along-track
  (early arrival)**, and the recurve geometry then projects part of that overshoot onto the
  cross-track axis as a spurious *west* value (Michael);
- a moderate recurver → modest, fairly clean west residual (Laura).

**This is seductive and I want to flag it as a hazard, not hand it over as a finding.** A bias that
can produce east, west, or pure-timing errors "depending on circumstances" explains everything and
therefore predicts nothing. As written it is an unfalsifiable just-so story, and a good reviewer
(or Woodruff) should attack it on exactly that ground. The per-storm explanations above are
**post-hoc rationalizations**, not predictions the model made in advance.

There is a deeper and more honest point underneath it. A storm's landfall track error is a
**convolution** of at least three sources:
1. the characterized β-drift bias (real, isolated in the testbed),
2. errors in the imposed ERA5 steering, and
3. missed track features (e.g., Michael's and Ivan's unmodeled NE/NW recurves).

**Landfall data alone cannot separate these.** The testbed isolates source (1); the six landfall
errors are dominated by sources (2) and (3) in ways that differ per storm. The cluster claim
implicitly attributed the landfall cross-track to source (1) — and that attribution is precisely
what failed out-of-sample. The defensible statement is narrower: source (1) is characterized; the
landfall errors are modest and storm-dependent and are *consistent with but not attributable to* it.

**What would turn the just-so story into a real claim** (open methodological question): take the
testbed β-drift vector and each storm's *actual* observed motion vector, predict the along/cross
split each storm *should* show if (1) were the dominant error, and test whether the six observed
splits match. If they do, "one bias, geometry-dependent projection" becomes a quantitative claim. If
they don't, (1) is simply a small term swamped by steering error, and the honest paper says so. We
have not done this; it may be the single most valuable next analysis.

**Metric note (Michael):** Michael's headline −75 km is partly an artifact of the same-latitude
metric. Its real error is along-track (+124 km *ahead* of obs at the landfall time); because obs
Michael is moving NNE there, the model's northward overshoot projects onto the cross-track axis as a
spurious west value. The run's own along/cross decomposition already shows this. Recommendation:
report the **along/cross decomposition for all six storms**, not a single cross-track scalar — it
makes Michael's error legible as timing and stops the metric from manufacturing a "west bias."

---

## The open question for review

**What is the paper's empirical claim now?** The working proposal, for the group to attack:

> Retire the eastward-cluster claim. Make the empirical headline **out-of-sample track accuracy** —
> a storm-agnostic, first-principles configuration with no per-storm tuning reproduces six historical
> landfalls (three of them never used in calibration) with cross-track errors of tens of km and
> characterized timing errors. Keep the testbed β-gyre over-rotation as a clean *characterized model
> property*, and describe its real-storm footprint honestly as modest, variable, and entangled with
> steering — explicitly **not** as a systematic eastward displacement.

Specific questions:
1. Does six storms (3 out-of-sample) suffice to **retire** the cluster claim with confidence? (We
   think yes.)
2. To **positively assert** out-of-sample accuracy as the new headline, is n = 3 out-of-sample
   enough, or do we need ~2–4 more (ideally a second basin, and at least one more slow storm and one
   more fast straight-mover)?
3. Is "one characterized bias, geometry-dependent landfall projection" worth pursuing **only if** the
   predicted-vs-observed along/cross test (above) passes — or should the paper not make a
   manifestation claim at all and simply report the testbed bias and the landfall errors as separate
   facts?
4. Does the testbed mechanism section stand fully on its own as a characterized-error contribution,
   independent of any landfall-attribution? (We think yes — confirm.)
5. Should cross-track be reported as a track-perpendicular along/cross decomposition for every storm
   (closing the Michael metric artifact, and incidentally enabling zonal movers like Andrew later)?

---

## Caveats / to-verify before acting

- **Michael's track plot** — confirm the −5.2 h / +124 km along-track is a clean "fails to follow the
  observed NE recurve and runs north early," not an uglier ERA5 steering problem. The signature is
  consistent with excess poleward drift on a slow storm, but eyeball it.
- The per-storm explanations in this note are **post-hoc**. Treat them as hypotheses.
- n = 6 is enough for the **negative** result (cluster isn't universal). A **positive** out-of-sample
  accuracy headline would be firmer with a few more storms.
- Numbers here are the same-latitude metric from the unified `run_storm` runs (logs checked in);
  along/cross figures are from each run's decomposition table as reported.

---

*Recommendation: do not rewrite the manuscript off this note. Synthesize external feedback first,
then reconvene to decide the reframe. The validation did exactly what it was for — three unseen
storms were run to test the cluster, and they answered. Better to learn it here than in review.*
