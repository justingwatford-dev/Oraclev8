# Paper outline — Oracle V8 track-error / clean-initialization paper

*Target: **BAMS** (Bulletin of the American Meteorological Society). Status: skeleton, 2026-06-24.
Spine = model-validation paper with two clean contributions (characterized β-gyre + blind track
skill), with the clean-initialization / compensating-errors methodology as a first-class section
(§3) and echoed in the intro and discussion. Section numbers here are the ones the finished
subsection's cross-refs resolve to (see the map at the bottom).*

BAMS notes: written for a broad audience — narrative register, significance-forward, lighter on
dense math than JAS/MWR. Plan for a **Capsule** (~30 words, under the title), a **Significance
statement**, and an optional **sidebar** for the dynamical-core internals so the main text stays
accessible. Rough length target ~6000–8000 words + figures.

---

## Working title (pick one — #1 recommended)

1. *Catching Our Own Errors: Clean Initialization Reveals an Over-Rotating β-Gyre and Blind Track
   Skill in a Tropical-Cyclone Model*
2. *Letting the Physics Miss: Clean Initialization and Compensating Errors in Tropical-Cyclone
   Track Prediction*
3. *What a Hurricane Model Gets Wrong on Purpose: Characterizing an Over-Rotating β-Gyre*

## Capsule (~30 words, BAMS)

> A tropical-cyclone model initialized from observations with no track-fitted parameter reproduces
> six landfalls and exposes an intrinsic, over-rotating β-gyre — a model error characterized in the
> open rather than hidden behind compensating fudge factors.

## Abstract (beats to hit)

- Track error = self-propagation (β-drift) + environmental steering + unmodeled features; hard to
  separate from landfall data alone.
- We adopt **clean initialization**: initialize from observations, never tune to the hindcast, and
  let residuals show. One outer-structure parameter, calibrated to canonical β-drift in a testbed —
  **not** to any landfall.
- **Result 1 (model property):** the model's β-gyre over-rotates poleward — a bounded, characterized
  self-propagation aim bias, robust to structure, diffusion, and resolution.
- **Result 2 (blind skill):** six historical landfalls (three never previously run) reproduced to
  tens of km cross-track with no track-fitted parameter; landfall cross-track is steering-controlled
  and not robust to storm geometry; the β-gyre bias is a detectable but subdominant along-track term.
- **Methodological lesson:** a compensating-errors cascade — each fix exposed the next hidden error —
  and the epistemic probes that retired two wrong claims before submission.

---

## §1 Introduction
*Purpose:* motivate clean initialization; state the two contributions + the methodological lesson.
- TC track forecasting and the along/cross error decomposition; β-drift vs steering.
- The hazard of compensating errors / tuned hindcasts; the clean-initialization alternative.
- Roadmap: §2 model, §3 methodology, §4 the two results, §5–6 discussion/conclusions.
- *Feeds from:* CHEATSHEET intro; SIX_STORM_VALIDATION_NOTE framing.
- *~800–1000 words.*

## §2 The Oracle V8 model and experimental design
*Purpose:* the methods, written for a broad reader (push core internals to a sidebar).
- **§2.1 Dynamical core.** ⚠ *Needs your accurate model-class description* (equations, vertical
  structure). Known components to document: RK3 integrator; projection/Poisson solver; Helmholtz
  divergence damping; ν₄ (∇⁴) hyperdiffusion; intensity cap; surface drag; Newtonian cooling;
  Δx = 15.6 km, Δt = 30 s. → candidate **sidebar**.
- **§2.2 Vortex initialization and the single outer-structure parameter.** HURDAT2 initialization;
  Holland vortex + wind taper; the taper-start radius **calibrated to testbed β-drift, not landfall**
  (the no-landfall-tuned-parameter discipline; cite the gate-beta calibration). → **cross-ref #3**
- **§2.3 Environmental steering.** ERA5 DLM, 3–7° annulus (storm-core excluded), time-varying
  lockstep update.
- **§2.4 Track verification.** Along/cross decomposition at observed fixes; why the same-latitude
  scalar is retired (Michael artifact); `landfall_verify`.
- *Feeds from:* CHEATSHEET (core/cap/taper/steering/scoring); production_config.py; era5_steering.py.
- *~1200–1500 words (excl. sidebar).*

## §3 Clean initialization and the compensating-errors cascade  *(the epistemic spine)*
*Purpose:* how we got here — the methodological contribution the broad BAMS audience values.
- The principle: clean init exposes true residuals instead of concealing them behind compensating
  errors.
- **The cascade** (each fix exposed the next error): intensity runaway → cap; taper lag → big
  domain; overshoot → wind taper; Hugo's 48 km → frozen-steering deficit → steering port; Katrina's
  +14 km → bad-init westward displacement; the data layer itself → HURDAT2 verification; the numerics
  "ghost" → the eff-family/Helmholtz chase that was lattice flicker.
- **The epistemic probes** (skepticism that retired wrong claims *before* review): the f-plane null;
  the inertial-oscillation phantom "13% over-translation"; the angle-wraparound auto-read; the
  retracted "delayed-onset drift"; the **retired eastward cluster**; the **falsified structural-init
  bridge** (β-drift projection test). → **cross-ref #2**
- *Feeds from:* CHEATSHEET "Compensating-errors cascade"; HURDAT2_VERIFICATION.md; RED_TEAM_AUDIT.md;
  BETA_DRIFT_PROJECTION_TEST.md.
- *~1200–1500 words.*

## §4 Results  *(= the finished subsection — drops in whole)*
- **§4.1 An over-rotating β-gyre.** Testbed characterization: f-plane null; three levers (structure,
  diffusion, resolution); intensity-independence + poleward precession; Fig. (gyre_precession),
  Table 1. → **cross-ref #1** (the canonical magnitude band lives here)
- **§4.2 Six-storm landfall track errors.** Table 2 (along/cross, log-verified); cluster doesn't
  generalize; Katrina-vs-Laura falsifier; projection test; blind track skill; β-gyre bias
  subdominant.
- *Feeds from:* PAPER_track_error_characterization.md (complete).

## §5 Discussion
- What the over-rotation is (β-Rossby equilibration failure); relation to canonical β-drift and to
  BAM/VICBAR-style steering models.
- Limits of landfall attribution (the convolution; steering dominance).
- Generality: Oracle-specific, or a class of TC models? (open)
- What blind skill with zero track-fitted parameters implies for clean-init as a method.
- Future: more storms / 2nd basin; deferred pressure-Holland arm.
- *Feeds from:* SIX_STORM open questions; CHEATSHEET "what's next."
- *~1000–1200 words.*

## §6 Conclusions
- Two separate facts + no bridge; the methodological lesson. *~400 words.*

## Significance statement (BAMS, ~120 words)
- Plain-language: why characterizing a model's *own* error honestly, with no hindcast tuning, matters
  for trustworthy TC guidance.

---

## Cross-ref resolution map (closes the PAPER_track_error_characterization.md placeholders)

| Placeholder in subsection | Points to | Resolves as |
|---|---|---|
| `[:139]` "magnitude band (Section X, above)" | β-gyre testbed (same section) | §4.1 (or "above") |
| `[:207]` "compensating errors (Section X)" | methodology spine | §3 |
| `[:134]` "storm-agnostic configuration (Section X)" *(to add)* | init / methods | §2.2 |
| author-note "six-storm landfall table" | self-reference | Table 2 (moot) |

## Figures & tables inventory
- **Fig. 1** gyre_precession.png (β-gyre m=1 asymmetry) — have it.
- **Table 1** resolution convergence — have it.
- **Table 2** six-storm along/cross — have it (log-verified).
- *Candidate new:* schematic of the along/cross decomposition; a six-storm track-map panel;
  projection-test predicted-vs-observed scatter (from BETA_DRIFT_PROJECTION_TEST.md).

## Open inputs needed from author
- Accurate **dynamical-core description** for §2.1 (model class / equation set / vertical structure).
- Confirm **title** choice and whether to use a **sidebar** for the core.
- Reference details to verify (Holland 1983; Chan & Williams 1987; Fiorino & Elsberry 1989;
  Chan 2005) — vol/pages.
