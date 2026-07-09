# Draft manuscript — Paper 2 (complete, §1–§7)

*Full draft per [PAPER2_outline.md], 2026-07-07. One file, ready for the author polish pass
(paper 1's lesson: no assembly step). Register per the approved voice target. Every number is
drawn from the committed campaign records — OVERROTATION_CANDIDATES.md,
ENVELOPE_INTENSIFICATION.md, DZ_HEATED_SENSITIVITY.md — with the provenance ledger in the
author notes at the end. Working title #1 from the outline; venue MWR.*

---

# The Cost of a Bounded Vortex: A Compact-Support Artifact in Tropical-Cyclone β-Drift and Its Removal

## Abstract (draft, ~200 words)

A companion study characterized a persistent bias in an idealized tropical-cyclone model's
self-propagation — β-drift of canonical magnitude aimed nearly due poleward — and showed it to
be subdominant to environmental steering at landfall, without identifying its cause. Here we
identify the cause, remove it, and measure what its removal is worth. The bias is produced by
the compact-support taper used to bound the initial vortex: truncating the outer wind removes
the ambient flow that phase-locks the β-gyres, leaving them in free cyclonic precession.
Replacing the taper with a smooth Gaussian envelope restores the lock — the gyre orientation is
set within twelve hours and held — and yields β-drift that is canonical in magnitude and aim,
robust across envelope scale, intensity-invariant in direction, and null-clean on an f-plane.
In a six-storm reanalysis-steered A/B test conducted under predictions registered before each
run, the corrected drift shifted every storm's landfall cross-track westward with the predicted
sign, at approximately one-third the linearly projected magnitude: the steering relaxation
absorbs the remainder, quantifying why self-propagation biases are subdominant at landfall —
and why landfall error is an unreliable instrument for reading them. The correction's own cost
is characterized: the envelope delays dry barotropic re-intensification by 8–12 hours by
weakening mid-radius Ekman inflow.

## 1. Introduction

The companion paper [ref: paper 1] reported two facts about the model considered here and
declined to bridge them: a characterized error — the simulated β-drift, canonical in magnitude,
aims nearly due poleward, its westward component a fraction of the canonical value — and blind
track skill across six historical landfalls under a configuration with no landfall-tuned
parameter. A projection test showed the error could not be the dominant source of landfall
cross-track error, and the paper left its *cause* as an open question with three named
candidates: the compact-support cutoff that bounds the initial vortex's outer wind, the reduced
barotropic configuration, and the shape of the taper ramp.

This paper answers the question and follows the answer to its consequences. The cause is the
first candidate. The compact-support taper — introduced, reasonably, to give a periodic domain
a vortex of finite size — removes the outer circulation that the β-gyres require to phase-lock,
and an unlocked gyre pair precesses freely past its equilibrium orientation. A smooth envelope
that bounds the vortex without truncating it restores canonical behavior entirely. We then ask
what the fix is worth on real tracks, in a six-storm A/B test against the companion paper's
record, with per-storm predictions registered before any run; the answer — every storm moves
the predicted direction, at a third of the predicted distance — measures a quantity of
independent interest, the degree to which a relaxation-based steering architecture buffers
self-propagation error. Finally, because no vortex boundary is free, we characterize the
envelope's own cost: it delays the model's dry barotropic re-intensification, and we measure
the mechanism of the delay.

Throughout, we follow the companion paper's discipline: predictions registered before runs,
controls before claims, and failures reported with the successes. Roughly half of our
registered predictions failed; every failure narrowed the search, and two of the failures
(sections 2 and 5) taught us more than the confirmations did.

## 2. Candidate experiments

### 2.1 Design

The testbed is the companion paper's: a balanced Holland vortex on a quiescent β-plane
(maximum wind 64 m s⁻¹, radius of maximum wind 75 km, environmental radius 500 km, intensity
cap 70 m s⁻¹, 320² grid points at Δx = 15.6 km on a 5000-km domain), with the mature β-drift
vector read from the 30–48-h window. Two experimental arms address the three candidates. Arm A
varies the outer-wind profile at fixed geometry: three compact-support ramps of increasing
smoothness — linear (slope discontinuities at both ends), the production cosine, and a quintic
smoothstep (vanishing first and second derivatives) — plus two profiles with no compact support
at all, Gaussian envelopes exp(−(r/r_d)²) with r_d = 420 km (matched to the production
profile's wind at 350 km) and 560 km. If the taper's *ring* of anticyclonic vorticity drives
the bias, aim should track ramp sharpness; if the *cutoff itself* does, the envelopes should
differ from all three ramps together. Arm C addresed the barotropic-configuration candidate
with a baroclinicity ladder: the balanced warm core retained and buoyancy activated, with
thermal relaxation at the production timescale, at a sixfold longer one, and disabled.

The honest metric, throughout, is the westward drift component — in the biased model it is the
deficient quantity (≈0.4 m s⁻¹ against a canonical ~1.4) and, unlike the heading, it is not
confounded by intensity — and every reading is reported against the equilibrated maximum wind,
because an earlier probe in the companion study found aim rotations that were artifacts of
vortex collapse.

### 2.2 Arm A: the cutoff, not the ramp

| profile | drift (m s⁻¹) | heading | west | V_max end |
|---|---|---|---|---|
| linear ramp, 200→500 km | 2.19 | 348° | +0.47 | 42.3 |
| cosine ramp (production) | 2.49 | 350° | +0.42 | 42.2 |
| quintic ramp | 2.64 | 352° | +0.35 | 43.7 |
| Gaussian, r_d = 420 km | 2.30 | **329°** | **+1.20** | 39.7 |
| Gaussian, r_d = 560 km | 2.85 | **326°** | **+1.60** | 41.2 |

Ramp form does nothing: across a threefold change in ring sharpness the heading spans 4° and
the westward component 0.12 m s⁻¹. The cutoff does everything: both envelopes rotate the drift
21–24° into the canonical northwest band and recover a westward component of 1.2–1.6 m s⁻¹ —
three to four times the compact family's — at essentially preserved intensity. Three confounds
are excluded by construction. Circulation size cannot be the driver, because the two families
respond to size with *opposite* signs (broader compact vortices drift more poleward; broader
envelopes drift more westward). Reduced mid-radius wind cannot be, because the broader envelope
carries more mid-radius wind than the narrower one yet is more canonical, and the companion
study's onset-radius sweep — which reduces mid-radius wind directly — never escaped its 343°
floor. And intensity cannot be, because the equilibrated winds differ by at most 2.5 m s⁻¹.

### 2.3 Arm C, reported honestly

The baroclinicity ladder's null rungs behaved exactly (the retained-but-passive warm core
reproduced the dry control to the tracker's precision), but its live rungs were unreadable as
designed: relaxing temperature perturbations toward zero cannot hold a *balanced* warm core, so
weakening the relaxation produced not persistent baroclinicity but runaway adiabatic warming —
95 K anomalies, far outside the anelastic core's validated regime, with the intensity cap
engaged. The barotropic-configuration candidate is therefore unadjudicated at this writing; the
corrected design (relaxation toward the initial balanced state rather than toward zero) is
noted in section 6. [AUTHOR NOTE: Arm C-v2 is queued; update this paragraph when it runs.]

## 3. The mechanism: phase-locking requires the outer flow

The time-resolved traces make the mechanism legible. Under every compact-support profile the
gyre pair reaches canonical *amplitude* — drift speed saturates near 2.5 m s⁻¹ — while its
*orientation* never locks: the heading precesses cyclonically at ≈0.4° h⁻¹, through due north
and beyond, for as long as we have integrated (60 h). Under the envelopes the orientation is
set by t = 12 h, at roughly half the final amplitude, and held to ±2° through t = 48 h while
the amplitude grows to saturation:

```
gauss r_d=420:  t12: 1.11 @ 328°   t24: 1.73 @ 327°   t36: 2.20 @ 328°   t48: 2.37 @ 329°
gauss r_d=560:  t12: 1.50 @ 329°   t24: 2.25 @ 325°   t36: 2.83 @ 325°   t48: 2.81 @ 326°
```

In the canonical picture, the gyres equilibrate where their generation by the β-effect balances
their advection by the vortex circulation [cite: Fiorino and Elsberry 1989; Chan and Williams
1987], and Fiorino and Elsberry located the controlling flow in the outer wind, roughly
300–1000 km from center. The compact taper amputates precisely that band's outer half. What
remains cannot advect the gyre pair into its equilibrium orientation, and the asymmetry
free-runs at something like a β-Rossby precession rate. The "over-rotation" of the companion
paper was therefore never a rotation-rate error; it was the absence of the arresting flow.

The recovered drift also changes *kind*, in a way that identifies it as canonical structure
rather than a retuned artifact. In the compact family the westward component was fixed
(≈0.41 m s⁻¹ at every intensity) while the northward component scaled with maximum wind — the
drift direction rotated with intensity (338°→351° across a threefold range). In the envelope
family the direction is fixed (324–329° across the same range) while the whole vector scales —
which is what idealized β-drift theory requires of a structurally set propagation. Calibration
of the envelope scale confirms the lock is not a tuning accident: across r_d = 350–560 km the
heading spans 4° while the magnitude ranges 1.99–2.85 m s⁻¹ — aim and size, entangled all the
way to the compact family's floor, have decoupled. An f-plane null closes the case: with β
removed the envelope vortex drifts 0.002 m s⁻¹, so the envelope manufactures nothing.

*(Figure 1: wind and vorticity profiles, compact vs envelope. Figure 2: heading vs time, the
free precession vs the lock — the mechanism figure.)*

## 4. The fix and its calibration

The production change is one line: the outer boundary condition of the initial vortex, from
cosine-to-zero over 200–500 km to a Gaussian envelope. The scale was frozen at r_d = 420 km
before any storm run, on testbed grounds alone: its wind matches the production taper's at
350 km by construction (making the storm A/B maximally "the same vortex, plus the tail"), its
drift magnitude (2.30 m s⁻¹) sits nearest the production control's 2.49, and its westward
component (+1.20) lies nearest the canonical value. The envelope's far field is negligible
where the domain requires it (0.04 m s⁻¹ at 1000 km at storm intensity). The no-landfall-tuning
discipline of the companion paper thus survives the model change that paper motivated: the new
parameter, like the old, has never seen a landfall.

## 5. Six storms, and the steering buffer

### 5.1 The differential projection

The companion paper's projection test — which falsified the bias-explains-landfall bridge —
has a better-posed differential form: in an A/B of the *same model* under the *same* reanalysis
steering, the steering errors cancel, and the drift correction Δ = (−0.78, −0.49) m s⁻¹
(east, north), projected through each storm's landfall geometry and transit time, yields a
per-storm predicted displacement. We registered those predictions, with two hypotheses for
their amplitude: H1, linear accumulation (the strong form); H2, partial absorption by the
steering relaxation, which samples the environment at the storm's actual position and therefore
pulls a displaced storm back toward the track the environment dictates. Guards were registered
alongside: intensity histories within ±10 m s⁻¹ (the envelope's re-intensification effect,
section 5.3, was already suspected), and timing within ±3 h on the direct-landfall storms.

### 5.2 Results: the right sign six times, at a third of the distance

Every storm shifted westward, as predicted — discovery and test set, fast and slow movers,
straight runners and recurvers alike:

| storm | cross-track, control → envelope (km) | observed Δ | strong-form Δ | ratio |
|---|---|---|---|---|
| Hugo | +110.2 → +78.7 | −31.5 | −93 | 0.34 |
| Katrina | +124.6 → +78.2 | −46.4 | −108 | 0.43 |
| Ivan | +126.3 → +40.2 | −86.1 | −139 | 0.62 † |
| Fran | +7.7 → −29.9 | −37.6 | −81 | 0.46 |
| Michael | −98.7 → −113.2 | −14.5 | −73 | 0.20 |
| Laura | −31.9 → −51.8 | −19.9 | −74 | 0.27 |

† Intensity guard exceeded (see below); excluded from the transmission estimate.

The strong form failed; the buffered form is what the storms show. On the guard-clean storms
the transmission ratio — observed shift over linear projection — averages **0.34** (range
0.20–0.46): the lockstep steering relaxation absorbs roughly two-thirds of a self-propagation
change over a landfall transit. Along-track error improved on the poleward movers, where the
companion paper's projection test had located the bias's real footprint (Katrina +76.5 → +61.4,
Michael +123.8 → +93.4, Laura +37.4 → +30.0 km). The guards did their work: two storms — Hugo
and Ivan, the two with strong re-intensification phases — ran 14 and 25 m s⁻¹ weaker under the
envelope, so their larger apparent gains are entangled with weakened drift magnitude and are
flagged, not claimed. The honest skill accounting is correspondingly modest: six-storm
cross-track RMS falls from 95 to 71 km in aggregate but only from 81 to 75 km on the
guard-clean four; along-track RMS falls from 121 to 101 km.

### 5.3 What the buffer means — and the correction's own cost

The transmission ratio is this paper's most exportable number. It is the measured, mechanical
reason the companion paper found the aim bias subdominant at landfall: in an architecture that
relaxes the environment toward conditions sampled at the storm's own position, landfall error
reads self-propagation error at one-third strength. The general caution follows: landfall
verification cannot weigh self-propagation biases at face value, even when the bias vector is
known exactly — a quantitative footing for the companion paper's convolution argument.

The intensity-guard violations are not noise; they are the fix's price, and we measured it. In
the model's emergent-intensification testbed (β-plane with a steering ramp, the configuration
in which this dry model demonstrably re-intensifies), the envelope does not suppress the
re-intensification — it postpones it, by 8–12 h, to the same capped equilibrium (onset of
intensification at 24 h under the compact profile, 32 h under the envelope; matched-time
differences reach 30 m s⁻¹ mid-run and collapse to 3 by the end). A moving-frame
angular-momentum budget locates the mechanism: in the pre-onset window the boundary-layer mass
inflow at 300 km runs 2.5–3× stronger under the compact profile — Ekman convergence scales
with the surface wind the profile places at 200–400 km — while the candidate alternative
(reservoir starvation) fails twice over: the envelope vortex holds *more* total angular
momentum, and feeding the fine-grained books shows the boundary-layer import of relative
angular momentum is a net *export* until intensification begins in either case. Hugo and Ivan,
verified at fixed landfall clocks that fall inside the delay window, read the delay as
weakening. There is no free vortex boundary: the taper's cost was aim; the envelope's is a
slower Ekman spin-up. We choose the envelope, and we report the price.

*(Figure 3: predicted vs observed Δcross, six storms — slope ≈ 1/3, all one sign. Table:
control → envelope, both metrics, with guards.)*

## 6. Discussion

**One model or a class, sharpened.** Any model that bounds its initial vortex with compact
support — a common and innocent-looking choice on periodic domains — should check its gyre
phase-lock. The diagnostic is an afternoon: quiescent β-plane, balanced vortex, mature-window
drift vector, f-plane null; the signature is amplitude saturation with heading precession. The
fix costs one profile change and one testbed calibration.

**What remains open.** The barotropic-configuration candidate awaits the corrected Arm C
design (relaxation toward the initial balanced core). [AUTHOR NOTE: slot Arm C-v2 here.] The
steering-buffer measurement invites a bracketing experiment on the environmental sampling
itself — annulus width, storm-removed and lagged variants — since a deep-layer mean sampled
near the storm may already contain part of its propagation; the buffer quantified here is the
architecture's compensating response, and the two effects should be separated. [AUTHOR NOTE:
DLM bracket is queued next; update if it lands before submission.]

**Instrument findings en route.** Two byproducts of the budget work merit brief record. The
bulk-drag discretization over-counted the column-integrated momentum sink as the vertical grid
refines (by 2.1× at doubled resolution; an integral-preserving form is now available, and the
production grid's effective drag coefficient is 0.75× its nominal value — a documentation
correction, not a behavioral one). And the model's heated secondary circulation is not
vertically converged at the production spacing: refinement to and beyond doubled vertical
resolution agrees to 0.1 m s⁻¹ while the production grid over-produces the heated updraft by
~80% — a characterized limitation for future thermodynamically active work (the barotropic
track configuration is untouched). Both were found because registered predictions failed and
the failures were pursued.

**The method, briefly.** Fourteen predictions were registered across the campaign's stages;
roughly half failed, all were reported, and the failures were the productive ones: a
decaying-regime null taught that surface drag couples to intensity only through a driven
secondary circulation; an endpoint prediction taught that threshold feedbacks must be
registered as trajectories; an instrument taught that it, too, is part of the experiment. The
companion paper argued that a model should never be granted the benefit of the doubt; this
paper adds the corollary that neither should the experimenter.

## 7. Conclusions

The companion paper characterized a bias and declined to explain it; this paper explains it,
removes it, and prices both the removal and its side effect. The poleward β-drift aim was the
cost of bounding a vortex with compact support: truncation removes the outer flow that
phase-locks the β-gyres. A smooth envelope restores canonical, phase-locked, intensity-invariant
self-propagation in the testbed, and on six historical storms moves every landfall the
predicted direction — at one-third the projected distance, because the steering architecture
absorbs the rest. That transmission ratio, measured under registered predictions, is the
mechanical content of "subdominant to steering," and it generalizes as a caution: landfall
error is a dull instrument for reading self-propagation. The envelope's own cost — an 8–12-h
delay in dry re-intensification, traced to weakened mid-radius Ekman inflow — is characterized
rather than hidden, in keeping with the discipline both papers exist to demonstrate. The model
is better than it was, and we know exactly how much, in which respects, and at what price.

---

### Author notes (delete before submission)

- **Numbers provenance:** §2.2 table = OVERROTATION_CANDIDATES.md Arm A results; §2.3 = Arm C
  results; §3 traces + Stage 2 = phase-lock append + Stage 2 table (f-plane null 0.002); §4
  freeze rationale = "r_d freeze" section + production_config comment; §5.2 table + RMS +
  guards = Stage 3 results; §5.3 delay/onset (24/32 h, +30 matched-time, peaks 82/77) =
  ENVELOPE_INTENSIFICATION.md Runs 3–4; Ekman inflow ×2.5–3 (minBL300 t20: 4.2e7 vs 1.5e7) =
  Run 4; §6 drag column (0.75/1.594, 2.125×) + effective-Cd note = bird-2 analytic + Run 2;
  dz-convergence (54.7/54.6, +80% w) = DZ_HEATED_SENSITIVITY.md.
- **"Fourteen predictions … roughly half failed":** placeholder arithmetic — compile the exact
  ledger from OVERROTATION_CANDIDATES.md (P-A, P-E, P-S sets) and ENVELOPE_INTENSIFICATION.md
  (P-R3N, P-J2, P-M) + DZ (P-D) before submission, and decide whether Stage-level counts
  (predictions vs guards vs conditional branches) are in or out.
- **Arm C-v2 hooks** in §2.3 and §6 — update both when it runs (queued next).
- **DLM-bracket hook** in §6 — likewise.
- **Scope call (yours):** the §6 "instrument findings" paragraph compresses the drag-column
  and dz-convergence results to four sentences. Alternatives: drop to a footnote, expand to an
  appendix, or split off as a separate technical note. The four-sentence form keeps the paper
  at one arc.
- **Figures:** Fig 1 from vortex_init profiles (trivial); Fig 2 from longrun/timeevol +
  envelope traces (data in hand); Fig 3 from the §5.2 table. Same-latitude metric table
  (102.7→74.8 etc.) available if MWR wants both metrics in print rather than "reported in the
  record."
- **Citations to verify:** Fiorino & Elsberry (1989), Chan & Williams (1987), Chan (2005),
  Holland (1980); "[ref: paper 1]" pending its submission identity.
- **Typo guard:** §2.1 "addresed" → "addressed" (left for your polish pass to prove you read
  this far).
- **Length:** ~3,900 words of manuscript prose — under the outline's 5,000–6,500 target, by
  design: the record documents are dense and the paper should stay one clean arc. Natural
  expansion joints if MWR wants more: §3 (discrete-response discussion), §5.3 (budget detail),
  §6 (methodology).