# Draft section — §3 Clean initialization and the compensating-errors cascade

*Drop-in for §3 per [PAPER_outline.md]. Manuscript register, BAMS narrative voice; the bracketed
author notes at the end are for you, not the paper. Drafted 2026-07-01 from CHEATSHEET
("Compensating-errors cascade"), HURDAT2_VERIFICATION.md, SIX_STORM_VALIDATION_NOTE.md,
BETA_DRIFT_PROJECTION_TEST.md, and RED_TEAM_AUDIT.md. Cross-refs use the outline's numbering
(§2.2, §4.1, §4.2); the Results subsection's "compensating errors (Section X)" placeholder
resolves to this section.*

---

## 3. Clean initialization and the compensating-errors cascade

### 3.1 The principle

A landfall hindcast can be right for the wrong reasons. A storm's track error at landfall is a
convolution of at least three sources — the storm's self-propagation, the imposed environmental
steering, and track features the model cannot represent — and a pair of compensating errors in any
two of them will produce a small landfall miss for free. Worse, when a model is developed against
the same landfalls it is judged on, such cancellations are not merely possible but *selected for*:
every tuning step migrates the configuration toward settings where errors offset, and the model's
failure is silently deferred to the first storm where the cancellation does not hold. Nothing in
the final error statistics distinguishes a model that is right from one whose errors are balanced.

We therefore adopted a discipline we call *clean initialization*: (i) every storm is initialized
directly from the HURDAT2 best track, with the initialization values themselves verified against
the primary source file; (ii) all storms run under one storm-agnostic configuration, with no
per-storm adjustments of any kind; (iii) the configuration's single structural free parameter —
the radius at which the outer wind profile begins to taper — is calibrated in a
quiescent-environment testbed against the canonical β-drift range from the published literature
(Section 2.2), never against any landfall; and (iv) residuals are reported as findings, not
absorbed into the configuration. Under this discipline the model has nowhere to hide an error.
The corollary, which this section illustrates, is that development proceeds by *exposing* errors —
including several that flattering early results had concealed, and several introduced by our own
tooling and data handling.

### 3.2 The cascade

It is one thing to state the hazard of compensating errors and another to watch it operate. The
model's development ran as a cascade: seven times, repairing one error exposed the next error it
had been hiding. We recount the sequence because the sequence is the point — no single audit,
however careful, would have found the seventh error while the first six stood in front of it.

1. **A numerical intensity runaway.** Early integrations developed grid-scale perturbation winds
   exceeding 150 m s⁻¹ — numerical growth, not intensification. A bounded relaxation toward a
   70 m s⁻¹ ceiling contained it, and simulated intensities became honest. Only then was the track
   readable at all.

2. **A domain boundary acting as a brake.** With intensity honest, the storm arrived hours late
   (+8.9 h). The cause was geometric: the vortex was feeling the domain-edge zones where the
   β-plane approximation is tapered off. Enlarging the domain removed the brake — and the storm
   now arrived hours *early*. The lateness had been masking an overshoot.

3. **A vortex with no outer edge.** The overshoot traced to vortex structure. The analytic wind
   profile used at initialization decays so slowly with radius that, untruncated, it filled the
   entire domain; the nominal environmental radius shaped only a dynamically passive pressure
   integral. The oversized circulation inflated the storm's β-drift. A cosine taper on the outer
   wind fixed this — and, in doing so, created the model's one genuine structural parameter, the
   taper-onset radius (Section 2.2).

4. **Our best result dissolves.** Hugo's celebrated 48-km landfall error — the strongest early
   validation — did not survive the honest vortex: re-run with the structural fixes, Hugo missed
   by +182 km east under the then-static steering. The 48 km had been a cancellation between an
   oversized vortex's inflated westward β-drift and a westward deficit in steering held frozen in
   time. The repair was time-varying environmental steering, sampled from ERA5 along the model's
   own track.

5. **A second flattering number dissolves.** Katrina's +14-km landfall error repeated the
   archetype. Scored against verified best-track data, the initialization had placed the storm
   98 km west of its true track; an eastward drift of ~1 m s⁻¹ then crossed the observed track at
   t ≈ 21 h, and the two errors near-cancelled at the landfall hour. Two errors, one flattering
   number — found only because the reference data themselves were re-verified (chapter 6).

6. **The data layer itself.** Executing the verification that the initialization module's own
   documentation demanded revealed that the stored initialization values matched no fix in the
   HURDAT2 file: positions matching time-shifted interpolants (one a +4.2-h along-track head
   start), intensities taken from fixes twelve hours later, central pressures matching nothing,
   one storm filed under the wrong cyclone identifier — and the pipeline's "observed" reference
   tracks were synthetic, one of them literally a straight line from the (incorrect) initialization
   point to landfall. The provenance pattern — plausible, internally consistent, and wrong — is
   characteristic of values recalled from memory rather than read from the source. Every
   downstream score had inherited these references. All initializations and reference tracks were
   re-derived from the primary file, and the affected results re-scored or retracted (§3.3).

7. **The instruments.** An apparent 13% over-translation — the model seemingly outrunning its own
   imposed steering — launched the longest hunt of the project, through the advection scheme, the
   pressure solver, and the damping operators. A Galilean control finally ended it: a balanced
   vortex in a uniform background flow, which must translate exactly with the flow, instead
   *appeared* to jitter — because the vortex center-finder reported grid-snapped positions with a
   ±1.5-cell (±23 km) flicker, a noise floor beneath which the entire family of over-translation
   measurements had been made. With sub-cell center-finding the control translates faithfully to
   within 1%, even while self-intensifying, and every over-translation number was retired as
   instrument artifact. The model's advection had been faithful all along; the residual that
   remained, after the instrument was fixed, was the genuine β-drift signal of Section 4.1.

The cascade descended through four strata: model physics (1), configuration geometry (2, 3),
the data layer (5, 6), and finally the measurement instruments themselves (7). Each error was
invisible while the ones above it stood. One chapter also ended differently from the others:
the outer-wind taper that repaired chapter 3 turned out to *mis-calibrate* the β-drift — too
strong and aimed too far poleward — which is not a bug but a physics trade-off, and it was
resolved not by tuning to landfall but by calibrating the taper-onset radius against the canonical
β-drift band in the testbed (Section 2.2). What remained after that calibration is the bounded,
characterized aim residual reported as Result 1 (Section 4.1).

### 3.3 The probes: retiring our own claims

Clean initialization removes the temptation to tune; it does not by itself protect against
over-interpretation. For that we relied on a second habit: subjecting each claim to the cheapest
test that could kill it, and registering predictions before runs rather than after. Several
claims did not survive, and their retirements shaped the paper more than the results that stood.

**Registered predictions.** Before re-running the storms with corrected initializations
(chapter 6 above), we wrote down predicted landfall bands — Katrina +80 to +130 km east, Hugo
+90 to +130 km east — derived from the decomposed error budget, so that the re-runs could
confirm or refute the budget rather than merely produce new numbers. Both re-runs landed inside
their bands. The residual they exposed was real, coherent — and, as it later proved, still
misattributed (below).

**Small retractions.** An early "13% over-translation" was contaminated by a rotating-frame
subtlety in the test harness (an unreferenced background flow undergoes inertial oscillation)
before being retired altogether by the Galilean control of chapter 7. A "delayed-onset drift"
finding was retracted when re-scored against verified best-track references: the onset had been
an artifact of the synthetic reference track, and the drift is in fact steady. An automated
sweep-reader once reported a heading change from 350° to 6° as a rotation *toward the northwest*;
6° is east of north. Each of these is small; each would have survived into review unkilled if the
claim had not been re-derived from primary data.

**The f-plane null.** The characterized self-propagation bias (Section 4.1) rests on a control:
with the planetary-vorticity gradient set to zero and all else identical, the vortex does not
translate at all (residual drift ≲ 0.02 m s⁻¹). The drift is genuine β-gyre dynamics, not an
advection or tracker bias — the instrument fix of chapter 7 made this null test meaningful.

**The two large retirements.** The first concerned our own headline. The three discovery storms
(Hugo, Katrina, Ivan) showed a tight, same-signed eastward cross-track cluster at landfall, and
the manuscript draft claimed it as systematic. Three test storms the configuration had never
run — Fran, Michael, and Laura — answered: +23, −75, and −39 km. The cluster was a property of
the discovery set, and the claim was retired (Section 4.2). The second concerned the bridge
between our two results. If the testbed β-gyre bias were the dominant source of landfall
cross-track error, a single bias vector projected through each storm's landfall geometry should
reproduce all six observed errors — a one-afternoon calculation. It fails decisively: it predicts
eastward cross-track for every poleward-moving storm, whereas Katrina and Laura approach on the
same heading with opposite observed signs (+125 vs −32 km) — no geometry-projected vector can do
that. The same test shows where the bias *does* live: it predicts Michael's along-track error
almost exactly (+123 predicted, +124 observed). Steering controls landfall cross-track; the
β-gyre bias is real but subdominant, expressed mainly as along-track error on poleward movers
(Section 4.2). This cheap test also cancelled an expensive plan: a data-informed vortex-structure
initialization scheme, motivated by the bridge claim, was abandoned when the model's own
structure sweeps showed it would push every discovery storm *further* east — the wrong direction.

What this discipline buys is stated most honestly in the negative. The results of Section 4 are
not the claims we set out to make; they are the claims that survived. A characterized
self-propagation bias whose landfall footprint is bounded and subdominant, and blind track skill
under a configuration with no landfall-tuned parameter, are what remained after the cluster
claim, the bridge claim, the over-translation family, and two flattering landfall errors were
retired by our own tests. We offer the cascade and the probes as the methodological content of
this paper: not that the model is right, but that its errors are where we say they are.

---

### Author notes (delete before submission)

- **Word count:** main text ≈ 1500 words — at the top of the outline's 1200–1500 target. If it
  must shrink, the candidates are the "small retractions" paragraph (compress to two sentences)
  and cascade items 1–2 (merge).
- **Chapter 6 phrasing ("values recalled from memory").** The provenance pattern language is
  deliberately agent-neutral. If you want to state explicitly that parts of the pipeline were
  assembled with AI assistance and that this is a working example of why verification-against-
  primary-source protocols are non-optional in AI-assisted research workflows, this is the place —
  it is arguably a BAMS-worthy point on its own. Your call on how far to go; the section stands
  either way.
- **Numbers to verify against logs before submission:**
  - +8.9 h (chapter 2) = the capped small-domain Katrina run; −7.6 h early = the big-domain run.
  - Hugo +182 km east (chapter 4) = the frozen-steering re-validation run.
  - Katrina −98 km init offset / ~1 m s⁻¹ drift / crossing at t ≈ 21 h (chapter 5) =
    HURDAT2_VERIFICATION.md F-V6 (init cross-track −98 km; drift 1.0–1.2 m s⁻¹).
  - ±1.5-cell = ±23 km at Δx = 15.6 km; control faithful to <1% (measured 1.000 ± 0.007, net
    drift 0.6 km over 14 h) — trace re-run, V8.6.3.
  - Registered bands and outcomes: HURDAT2_VERIFICATION.md "Registered prediction" + scorecard
    (Katrina +126.5 in-band; Hugo corrected run +102.7 in-band — cite the corrected number).
  - f-plane null ≲ 0.02 m s⁻¹ = gate-beta-fdecomp.
  - Projection test values (+125/−32; +123 vs +124) = BETA_DRIFT_PROJECTION_TEST.md.
  - Six-storm same-latitude values (+103/+114/+68; +23/−75/−39) = SIX_STORM_VALIDATION_NOTE.md.
- **Terminology guard:** this draft deliberately avoids internal jargon (no version numbers, no
  harness mode names, no "eff_y"). "Discovery/test" labels used per the calibration-linchpin
  resolution (taper-start calibrated to testbed β-drift, not landfall) — do NOT reintroduce
  "in-sample/out-of-sample."
- **Optional addition:** the LH82 small-perturbation study (LH82_SMALL_PERTURBATION_FINDINGS.md)
  is a solver-level instance of the same discipline — a near-miss false claim ("the equation set
  breaks at 5%") caught by varying damping and time step independently. One sentence could go at
  the end of "small retractions" if you want a physics-level example alongside the track-level
  ones; alternatively it belongs with the §2.1 sidebar. Left out here to keep §3 track-focused.
- **Cross-ref wiring:** the finished Results subsection's "(Section X)" compensating-errors
  placeholder → this section (§3). This draft's forward refs: §2.2 (taper calibration ×3),
  §4.1 (β-gyre result ×3), §4.2 (six-storm + projection ×3). Renumber if the outline shifts.
- **Fran's +23** is quoted from the six-storm note's same-latitude table; the projection-test doc
  lists Fran's landfall-fix cross as +7.7 km. Consistent (different metrics) — but pick ONE metric
  for §3's sentence and match whatever §4.2's table reports at that point in the text.
