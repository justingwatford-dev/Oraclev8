# Structural initialization — frozen recipe + pre-registration

*Specification for replacing the fixed-vortex constants (R_env = 500 km, taper_start = 200 km,
B = 1.5) with storm-specific **observed outer-wind structure**. Status after red-team review:
**recipe/data sheet only; directional predictions are not frozen and treatment runs are blocked**
until the hold items below are resolved.*

---

## Red-team hold — do not run treatment yet

The red-team review found real blockers in the interpretation layer. The structural numbers below
are useful, but the prediction sheet is suspended until these checks are complete:

1. **Component-decomposed size sweep:** run a fixed-Vmax gate-beta sweep over the treatment-relevant
   `R_env / taper_start / B` combinations and report north/west components, not just speed/heading.
   The current "broader → more east/right" predictions are not safe until we know whether size adds
   westward drift, poleward drift, or both at fixed Vmax.
2. **ERA5 DLM sensitivity:** production steering already uses a 3–7° annulus, not a point sample at
   the storm center, but that is not a storm-removed field. Before structural-init attribution,
   bracket Ivan at minimum with a wider annulus or storm-removed/lagged environmental steering test.
3. **Pipeline regression:** implement structural-init behind a control mode that forces
   `R_env=500 km`, `taper_start=200 km`, `B=1.5`, `Rmax=75 km`, then verify it reproduces the current
   control pathway before any treatment result is interpreted.
4. **Fixed domains:** hold each storm's control and treatment domain/grid fixed. Structural `R_env`
   must not silently change `choose_domain()` output in the primary comparison.
5. **Falsifiable predictions:** after items 1–2, replace the rough bands with one-sided pass/fail
   bounds, e.g. "cross changes by at least +20 km" or "absolute cross decreases by at least 20 km."
   If all treatment changes are within ±20 km, do not claim a landfall-error β-drift mechanism.

---

## Principle (the one line that keeps this honest)

> Fixed equations, fixed numerics, fixed data-ingestion recipe, **storm-specific observed outer
> structure.** Outer wind radii come from data observed **at initialization** — never from the
> landfall being predicted, and never hand-adjusted per storm after seeing a result. Rmax is the
> documented 75 km resolution floor unless an observed Rmax exceeds it.

This is not a retreat from "let physics do its thing." A fixed 500 km vortex is not storm-agnostic,
it is **storm-indifferent** — a wrong size applied uniformly. β-drift scales with vortex size, so a
constant size guarantees a wrong β-drift magnitude. Giving the physics a correctly-sized vortex is
*more* principled than the current bundle, not less — especially since `taper_start_frac = 0.40`
carries calibration-set fingerprints. This recipe **retires the 0.40 entirely**, for all six storms
including the in-sample three, so that after this there is no tuned knob and no in-sample/out-of-
sample asymmetry in *how the vortex is built*.

---

## The correctness guard — read this first

**B is defined against Rmax. They must be fit as a pair, at the resolved scale you actually
integrate.** The failure mode that would silently corrupt this whole exercise: fit B to a storm's
profile at its *observed* Rmax (e.g., Andrew's ~18 km), then integrate the vortex at the resolved
Rmax (75 km). That applies a shape parameter to a storm four times too wide, and the run completes
and returns a plausible-looking number. Never do this.

The recipe avoids it structurally: **Rmax is fixed first** (resolution floor), and **B is fit at
that fixed Rmax.** B and Rmax are consistent by construction because B is always fit at the Rmax used.

---

## Data sources (and the parity rule)

Structure comes from observed wind radii and the wind profile Oracle actually integrates:

- **2004+ storms** (Ivan, Katrina, Michael, Laura): HURDAT2 carries R34 / R50 / R64 in four
  quadrants — use directly. In the current Atlantic file, the trailing Rmax field is missing at
  these initialization fixes, so the 75 km resolved-Rmax floor applies unless EBTRK later proves a
  storm exceeds it.
- **Pre-2004 storms** (Hugo 1989, Fran 1996; Andrew 1992 if ever run): HURDAT2 has no radii. Use the
  **Extended Best Track (EBTRK, DeMaria; 1988–present)**, which carries Rmax, R34 quadrants, ROCI,
  and Penv for these storms.
  Current source file for this recipe pass:
  `EBTRK_AL_final_1851-2021_new_format_02-Sep-2022-1.txt` from CIRA/RAMMB. The wind-radii groups
  must be parsed as fixed-width four-by-three-character quadrant fields; whitespace splitting is not
  reliable for legacy records. Include only positive finite quadrant radii in the azimuthal mean;
  treat `0`, `-99`, and `-999` as unavailable/no-radius values and log the raw quadrants.
  EBTRK's four-digit source code applies to RMW, eye diameter, POCI, and ROCI, not directly to the
  34/50/64 kt wind-radii groups; still carry a per-record confidence flag for legacy storms.

**Parity rule (non-negotiable):** the pre-2004 storms get *real* EBTRK radii, not a climatological
fallback, specifically so that data quality does not become confounded with the in-sample/out-of-
sample split. Hugo and Fran must get measured structure like everyone else. A climatological fallback
exists (below) but is a last resort, flagged per storm, **not** the default for old storms.

---

## The recipe (deterministic; freeze before any rerun)

All winds are handled in the same wind space as the current Oracle vortex initializer. This is an
implementation guard: `HollandVortexInit` integrates a Vmax-normalized Holland profile using
`Vmax`, `Rmax`, and `B`; it does not currently integrate a pressure-Holland profile from Penv/Δp.
Therefore the primary structural fit below uses the current Oracle profile directly. Penv/Δp is
recorded as an audit field and possible future pressure-wind upgrade, not as a hidden input to this
treatment arm.

**Per storm, at the initialization fix:**

1. **Observed radii** = arithmetic mean of available positive quadrant radii at the initialization
   fix: R34, R50, and R64 where available. Record the raw quadrant values, missing-value handling,
   and asymmetry separately; the model vortex is axisymmetric, but the data are not.

2. **Rmax** = max(observed Rmax, **75 km** resolution floor). Fixed rule; for all six this is
   expected to be 75 km, but EBTRK must confirm this where HURDAT2 Rmax is missing. Document any
   storm where observed Rmax exceeds the floor.

3. **B** = fit the current Oracle Holland profile,
   `V(r) = Vmax * sqrt(x * exp(1 - x))`, `x = (Rmax/r)^B`, at the Rmax from step 2. The objective is
   nonlinear least squares in log-radius space:
   `min_B Σ_i [log(r_model(V_i; B)) - log(r_obs,i)]²`, where `r_model` is the outer-branch radius at
   threshold `V_i` for R34/R50/R64. Use only thresholds whose observed radius lies outside the
   resolved Rmax floor; an observed R64 smaller than 75 km is an inner-core/resolution conflict, not
   a usable outer-radius constraint. Clamp **B ∈ [1.0, 2.5]**. If fewer than two radii are usable,
   or if the fitted B hits either clamp, tag the storm **low-shape-confidence** and do not
   over-interpret B differences.

4. **Complete the profile** V(r; Vmax, Rmax, B) — this is the exact profile the vortex will integrate.

5. **taper_start** = the radius where V(r) = **17.5 m/s (34 kt)** — i.e., the model's own gale
   radius (≈ observed R34 by construction). This is where blending toward the environment begins.

6. **R_env** = the radius where V(r) = **5 m/s (10 kt)** — where the storm's wind is negligible
   against the environment. Clamp **R_env ∈ [400, 800] km** for domain/numerical sanity. This is
   where blending to the environment completes.

7. Done. R_env, taper_start, and B are all derived from **one** fitted profile, so they are mutually
   consistent and there are **no free multipliers**. Implementation still passes a
   `taper_start_frac = taper_start / R_env` to `HollandVortexInit`, but that fraction is now an
   output of the recipe, not the calibrated constant 0.40.

**Pressure-wind audit (not part of this treatment arm):** record Pc, Penv if available, and the
pressure-implied B/Vmax consistency check. A future pressure-Holland arm would require a separate
recipe and code change, because the current initializer does not use Δp to set the wind field.

**Known limitation, applied uniformly (document, don't hide):** the 75 km resolved Rmax means the
*inner core* (r < R34) is resolved-scale, not observed-scale; the fit prioritizes the *outer*
circulation that governs β-drift. For storms whose observed R34 approaches 75 km (very compact
storms — Andrew is the extreme case), the resolution floor prevents faithful representation and the
model vortex is necessarily broader than reality. This bites Andrew hardest (deferred anyway) and is
expected to be minor for the five poleward storms.

---

## Experiment design (planned two-arm after red-team hold)

- **Control:** current constant config (R_env 500 / taper 200 / B 1.5). **Already run** — the
  six-storm results we have.
- **Treatment:** full observed-outer-structure init per the recipe, **0.40 retired**, applied to
  **all six** only after the hold items are complete.
- **Domain/grid:** fixed per storm across control and treatment. Do not let treatment `R_env`
  silently choose a larger domain in the primary comparison.
- **Comparison metric:** the **along/cross decomposition for every storm**, not the same-latitude
  cross-track scalar. (This also closes the Michael metric artifact, where a northward overshoot
  projected onto the cross-track axis as a spurious west value.)

This is not "true observed vortex vs. wrong vortex." With a 75 km Rmax floor and frequent B-clamp
saturation, it is a **data-informed outer-size vortex vs. generic constant-size vortex**. That is a
legitimate experiment, but it must be described that way.

**Why no primary "size-only, keep B = 1.5" rung:** a naive version would staple an observed-size
outer profile onto an unfitted inner shape. A self-consistent size-only arm could be designed later
(for example, fit the one free radius at fixed B to match R34), but it is not the primary treatment.

*(Optional, only if size-vs-shape isolation is wanted later: a pure-size-scaling arm that scales the
**whole** consistent profile by observed/default size at fixed B — self-consistent, unlike step 2.
Not primary.)*

**This experiment also tests the convolution hypothesis directly.** The briefing note argued the
landfall error is β-drift bias + steering error + missed recurve, and that steering may dominate.
Structural init changes the **β-drift** component (size-dependent) while leaving steering untouched
(size-independent). So:
- if the landfall errors **move** with structural init → β-drift was a significant component, and
  the direction/magnitude should match the predictions below;
- if they **barely move** → steering/recurve dominates, β-drift is a small term, and the testbed
  β-gyre stays a clean characterized property with limited real-storm footprint.

Falsifiability commitment: if the treatment changes all storm along/cross errors by less than
±20 km, or if the pre-committed signs fail, do **not** claim a β-drift landfall-error mechanism.
Frame the result as a testbed β-gyre characterization plus six-storm track-error audit.

---

## Pre-registration: structural values only; predictions suspended

**Do not treat the direction/magnitude predictions as frozen.** Sequencing matters and is itself a
guard:

1. Apply the recipe → get (B, taper_start, R_env) per storm. *Deterministic from EBTRK/HURDAT2; not
   a model run.*
2. For each storm, compare to the constant (75 / 200 / 500 / 1.5): is the storm **broader** or more
   **compact** than the default vortex?
3. Run the hold-item size sweep and DLM sensitivity test. Then predict the **direction** and
   one-sided bound of the landfall-error change from the component-decomposed β-drift response.
   Convert that expected motion change into along/cross using each storm's observed landfall motion;
   do **not** assume that "stronger β-drift" always means "more east."
4. **Then** run the model. Compare actual to predicted.

The frozen recipe protects against tuning; the **written-down predictions** protect against post-hoc
rationalization. Two different traps, two different guards.

| Storm   | obs Rmax / R34 / ROCI | recipe R_env / taper / B | broader or compact vs 500/200/1.5 | predicted β-drift change | predicted landfall-error change (dir + rough mag) | — actual change — | match? |
| ------- | --------------------- | ------------------------ | --------------------------------- | ------------------------ | ------------------------------------------------- | ----------------- | ------ |
| Hugo    | EBTRK Rmax **30 nm** (**56 km**, floored to 75 km); R34/R50/R64 mean = **138 / 88 / 50 nm** (**255 / 162 / 93 km**). EBTRK POCI/ROCI fields are **15 / 15** with source `4444`; treat as legacy-unreliable audit fields, not structural inputs. | **R_env 721 km / taper 260 km / B 2.50 clamp**; low-shape-confidence because B saturates. | Broader outer vortex; taper +60 km, R_env +221 km. | **Suspended.** Need fixed-Vmax size-sweep N/W components. | **No committed prediction.** Existing +23 along / +110 cross is the control baseline only. | not run | TBD |
| Katrina | HURDAT2 Rmax missing; R34/R50/R64 mean = **120 / 76 / 50 nm** (**222 / 141 / 93 km**); ROCI not in HURDAT2. | **R_env 721 km / taper 261 km / B 2.50 clamp**; low-shape-confidence because B saturates. | Broader outer vortex; taper +61 km, R_env +221 km. | **Suspended.** Need fixed-Vmax size-sweep N/W components. | **No committed prediction.** Existing +77 along / +125 cross is the control baseline only. | not run | TBD |
| Ivan    | HURDAT2 Rmax missing; R34/R50/R64 mean = **188 / 101 / 74 nm** (**347 / 188 / 137 km**); ROCI not in HURDAT2. | **R_env 800 km clamp / taper 314 km / B 2.50 clamp**; low-shape-confidence and R_env saturated. | Much broader, but R_env-clamped. | **Suspended.** Need size-sweep N/W components and DLM sensitivity first. | **No committed prediction.** Existing +249 along / +126 cross is the control baseline only. | not run | TBD |
| Fran    | EBTRK Rmax **20 nm** (**37 km**, floored to 75 km); R34/R50/R64 mean = **181 / 124 / 73 nm** (**336 / 229 / 134 km**); POCI **1009 hPa**, ROCI **250 nm**; source `1411`. | **R_env 800 km clamp / taper 309 km / B 2.27**. | Much broader; taper +109 km, R_env +300 km, but R_env-clamped. | **Suspended.** Fran's fast motion makes along/cross projection especially uncertain. | **No committed prediction.** Existing about -45 along / +8 cross is the control baseline only. | not run | TBD |
| Michael | HURDAT2 Rmax missing; R34/R50/R64 mean = **128 / 50 / 28 nm** (**236 / 93 / 51 km**); R64 is inside the 75 km Rmax floor and excluded from fit; ROCI not in HURDAT2. | **R_env 663 km / taper 238 km / B 2.50 clamp**; low-shape-confidence due compact-core conflict. | Broader outer vortex, but inner-core structure unresolved. | **Suspended.** Projection is geometry-sensitive on Michael's NNE motion. | **No committed prediction.** Existing +124 along / -99 cross is the control baseline only. | not run | TBD |
| Laura   | HURDAT2 Rmax missing; R34/R50/R64 mean = **110 / 70 / 43 nm** (**204 / 130 / 79 km**); R64 barely clears the Rmax floor; ROCI not in HURDAT2. | **R_env 663 km / taper 238 km / B 2.50 clamp**; low-shape-confidence because B saturates. | Moderately broader. | **Suspended.** Need component response before claiming sign. | **No committed prediction.** Existing +37 along / -32 cross is the control baseline only. | not run | TBD |

*Computation note for the pre-run rows:* modern-storm radii are arithmetic means of nonzero HURDAT2
quadrants at the initialization fix; Hugo/Fran radii are arithmetic means of nonzero EBTRK quadrant
radii at the same fixes. Recipe values use the current Vmax-normalized Oracle Holland profile,
Rmax = 75 km, log-radius B fit over usable outer radii, and R_env clamped to [400, 800] km. Hugo's
POCI/ROCI fields are retained only as an audit warning because the EBTRK record gives an
unphysical `15 / 15` pair.

*Rule for the sheet: predictions are not committed yet. The structural values are committed; the
landfall-error signs/bounds become valid only after the red-team hold tests are documented.*

---

## Sanity-check results before implementation

- **EBTRK is now available for the old storms, but parser/source logging is a hard prerequisite
  before any treatment run.** Hugo and Fran have no HURDAT2 radii, so the implementation must ingest
  the fixed-width EBTRK quadrant groups and print the source code per field. Hugo's wind radii and
  RMW are usable; its POCI/ROCI pair is legacy-unphysical and must not silently feed a pressure or
  domain-size calculation. The modern storms' HURDAT2 Rmax fields are still missing at the chosen
  initialization fixes, so Rmax/ROCI/Penv parity still requires EBTRK ingestion or an explicitly
  documented "not available" tag.
- **The first pass is mostly a size experiment, not a clean shape experiment.** With the 75 km
  resolved-Rmax floor and the current Vmax-normalized Oracle profile, Hugo, Katrina, Ivan, Michael,
  and Laura fit to **B = 2.50**, the upper clamp. Fran is the exception at **B = 2.27**, which is
  useful but should still not be over-interpreted until the structural treatment is run. Clamp
  saturation should be reported as a resolution-floor consequence.
- **R_env clamping matters.** Ivan and Fran reach the 800 km R_env ceiling in the pre-run
  calculation. Their responses are therefore lower-bound responses to observed outer size, not fully
  free observed-size initializations.
- **The directional prediction is not yet justified.** Existing gate-beta notes show magnitude and
  heading move with taper-start radius, while the f-plane decomposition shows the westward component
  is small and nearly flat across Vmax. A storm-relevant size sweep must decompose N/W components
  before the recipe predicts along/cross signs.
- **ERA5 steering is annular but not storm-removed.** `ERA5Steering.get_dlm()` uses a 3–7° ring,
  which already avoids point-sampling the vortex core. That is better than the red-team's worst-case
  reading, but still leaves a real environmental-steering sensitivity to test.
- **The pressure-Holland version is a separate experiment.** If the team wants Penv/Δp to set the
  wind field, `HollandVortexInit` must be changed and a new prediction sheet frozen. Mixing
  pressure-fitted B with the current Vmax-normalized initializer would be an implementation mismatch.

---

## What success and failure look like (both reportable)

- **Errors tighten out-of-sample AND pre-committed signs/bounds match** → the strongest possible result: the
  "cluster broke" finding becomes *"the model exposed a missing structural initialization; observed
  structure recovers accuracy across all six storms with no per-storm tuning."* The reframe flips
  from a negative to a positive headline.
- **Errors don't move much** → steering/recurve dominates real-storm landfall error. Do not write a
  β-drift-as-landfall-error paper; write a testbed β-gyre characterization plus track-error audit.
- **Errors move but predictions fail** → the simple "β-drift ∝ size" picture is incomplete; report
  it as such (do not retro-fit an explanation).

---

## Honesty ledger (the guards, made explicit)

1. **Recipe frozen after hold tests**, before any treatment rerun. No "Michael should have a smaller
   taper" after seeing Michael.
2. **Predictions committed** before the treatment runs, but only after component-size and DLM
   sensitivity checks.
3. **Per-storm data source + confidence** recorded; no climatological estimate masquerading as
   measurement. The 75 km Rmax is a resolution floor, not an observation; EBTRK/HURDAT2 radii are
   the observed outer-structure inputs.
4. **Inner-core resolution limit** documented and applied uniformly (75 km floor), not silently.
5. The recipe **removes** the one calibrated knob (0.40); it does not add new ones (no free
   multipliers — R_env/taper fall out of the fitted profile).

---

## Code scope (what this actually requires)

This is a real code change, not a config tweak — flagging the surface so it can be planned:

- **Extend `hurdat2.py`** to retain R34/R50/R64 quadrant radii, not only Rmax. The data are present
  for the modern storms, but the current loader does not expose them.
- **EBTRK parser** for Rmax/ROCI/Penv and for the pre-2004 storms (Hugo, Fran): EBTRK format →
  Rmax, R34 quadrants, ROCI, P_env per fix. Modern-storm Rmax is also missing in the current HURDAT2
  file at these initialization times, so EBTRK is still useful beyond Hugo/Fran.
- **`structural_init()`**: implements steps 1–7; returns (Rmax, B, R_env, taper_start) per storm;
  emits a confidence tag and the source per field.
- **Retire `taper_start_frac = 0.40`** from `production_config.py`; route the vortex builder to the
  structural values. Keep the constant config available behind a flag as the **control** arm.
- Keep **domain/grid fixed per storm** across control and treatment. If a structural `R_env` would
  require a larger domain for numerical safety, that is a separate domain-sensitivity arm, not the
  primary comparison.
- Add a **pipeline regression mode** that routes the control constants through the new structural
  path and checks that the control result is reproduced before treatment runs are interpreted.
- **Log the structural values** per storm (checked in), so the treatment runs have the same
  provenance discipline as everything else.
- Consider a one-shot **structure audit** (`python structural_init.py <Storm> <ebtrk_or_hurdat2>`)
  that prints the observed radii and the derived (Rmax, B, R_env, taper_start) — the analogue of
  `python hurdat2.py <Storm>` — so the recipe output can be eyeballed before any model run.

---

*Recommendation: freeze this recipe and the prediction sheet, take both into the feedback round
alongside the six-storm note, and only then write the EBTRK parser and `structural_init()`. The
order is deliberate: recipe and predictions committed first, code second, runs third, comparison
last. Same discipline that just saved the ν₄ claim and retired the cluster — applied before the work
instead of after.*
