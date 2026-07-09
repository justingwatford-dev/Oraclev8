# Paper 2 outline — the over-rotation mechanism, its removal, and the steering buffer

*Sequel to the BAMS clean-initialization paper (PAPER_outline.md). Status: skeleton, 2026-07-03,
drafted the day the campaign completed. Spine = the question paper 1 deliberately left open
(§5.1's three candidates) answered end-to-end under registered predictions: mechanism →
calibrated fix → six-storm validation — plus the campaign's second, unplanned discovery, the
measured steering buffer. Register: same as paper 1 — drama in the facts, none in the prose.
Primary source for every number: OVERROTATION_CANDIDATES.md (all stages committed).*

Venue: **MWR** recommended (model diagnosis + forecast verification; the phase-lock dynamics
could go JAS, but the paper's arc — artifact, fix, validation, architecture measurement — is
MWR-shaped). Length target ~5000–6500 words. Standard sections, no sidebar needed.

---

## Working titles (pick one — #1 recommended)

1. *The Cost of a Bounded Vortex: A Compact-Support Artifact in Tropical-Cyclone β-Drift and
   Its Removal*
2. *Why the β-Gyres Never Locked: Outer-Wind Truncation and Self-Propagation Aim Error in a
   Tropical-Cyclone Model*
3. *Phase-Locking the β-Gyre: From Characterized Bias to Measured Fix*

## Capsule / abstract beats

- Paper 1 characterized a poleward β-drift aim bias as a bounded model property, showed it was
  subdominant to steering at landfall, and named three untested candidates for its cause. This
  paper tests them.
- **Result 1 (mechanism):** the bias is caused by the compact-support outer-wind taper. The
  β-gyre reaches canonical *amplitude* under any profile, but its *orientation* phase-locks
  only when ambient vortex flow exists at and beyond the gyre radii; truncation leaves the gyre
  pair in free cyclonic precession (~0.4°/h through north). A smooth Gaussian envelope restores
  the lock by t+12 h — canonical aim, held to ±2° for 36 h.
- **Result 2 (structural fix):** with the envelope, β-drift is canonical in magnitude AND aim,
  r_d-robust across a 60% scale range, intensity-invariant across a 3× Vmax range (the
  invariance *changes kind*: fixed-direction/scaling-magnitude, as canonical theory requires),
  and null-clean on the f-plane. One parameter, recalibrated in the testbed only — the
  no-landfall-tuning discipline held through the model change.
- **Result 3 (the steering buffer):** a six-storm A/B under predictions registered before any
  run: the cross-track shift had the predicted (westward) sign in ALL SIX storms, at ~1/3 of
  the linearly-projected magnitude — the lockstep steering relaxation absorbs ~2/3 of a
  self-propagation change over a landfall transit. This measured buffering is the mechanical
  reason self-propagation biases are subdominant at landfall, and it generalizes: landfall
  error is an unreliable instrument for reading self-propagation error even when the bias is
  known exactly.
- **The honest ledger:** skill improves modestly where it should (clean-storm cross RMS −8%,
  along-track improved on the poleward movers) and the guards caught what they were built for —
  two storms ran weaker under the envelope (it trims the mid-radius angular-momentum reservoir
  feeding dry barotropic spin-up), and their larger gains are flagged as intensity-entangled,
  not claimed. The bounded-vortex choice has no free option: the taper's cost was aim, the
  envelope's cost is intensification behavior. We characterize the trade rather than hide it.
- **Method:** every stage gated by registered predictions (five predictions confirmed, two
  failed or unreadable — reported); one live-caught data-layer regression (lost uncommitted
  configs) continuing paper 1's §3 record.

---

## §1 Introduction
*Purpose:* the question paper 1 left open; why "characterize, then hunt" beats "tune."
- Recap paper 1 in one paragraph: two facts, no bridge, three named candidates.
- The claim of this paper: all three candidates adjudicated; mechanism found; fix validated;
  a second quantity discovered en route (the buffer).
- *Feeds from:* paper 1 §4.1/§5.1; OVERROTATION_CANDIDATES.md preamble. *~500 words.*

## §2 Candidate experiments (Arms A and C)
- **§2.1 Design + registered predictions.** The ring-sharpness ladder (linear/cos/smooth5 at
  fixed geometry) + the no-cutoff Gaussian envelopes; the baroclinic persistence ladder; the
  honest metric (westward component, read against Vmax_end — the confound lesson inherited
  from paper 1's formulation probe).
- **§2.2 Arm A result.** Ramp *form* does nothing (4° span); the *cutoff* does everything
  (gauss rows +21–24° into the canonical band, west 0.42 → 1.20–1.60). Confound elimination:
  opposite size-trends in the two families; mid-radius wind ruled out; Vmax_end preserved.
- **§2.3 Arm C, reported honestly.** Nulls passed exactly; live rungs unreadable as designed
  (relax-to-zero cannot hold a balanced warm core; θ′ → 95 K, outside the validated LH82
  regime). Redesign (relax to θ′_ref) stated; candidate 2 neither implicated nor exonerated.
- *Feeds from:* OVERROTATION_CANDIDATES.md Arms A/C + prediction scorecards. *~1000 words.*

## §3 The mechanism: phase-lock requires the outer flow
- The traces: compact = amplitude saturates, heading marches through north at ~0.4°/h, no
  plateau; envelope = heading set by t12 at half amplitude, flat ±2° to t48.
- Physical reading vs Fiorino & Elsberry's outer-wind control (300–1000 km): the truncation
  amputates the ambient swirl that advects the gyre pair into its equilibrium orientation;
  "over-rotation" was never a rotation-rate error — it is the absence of the lock.
- The regime change in intensity dependence (fixed-west/scaling-north → fixed-direction/
  scaling-magnitude) as the signature of recovered canonical structure.
- *Feeds from:* phase-lock trace append; Stage 2 table; paper 1 §4.1 (the characterization
  being explained). *~800 words. Figure: heading-vs-time, compact vs envelope (THE figure of
  the mechanism half).*

## §4 The fix and its calibration
- Gaussian envelope; r_d sweep (aim/size decoupled — 4° across 350–560 km); intensity
  invariance; f-plane null; the freeze (r_d = 420 km, chosen on testbed grounds only, before
  any storm run — the discipline carried through the model change).
- *Feeds from:* Stage 2 + freeze rationale (production_config comment). *~600 words. Table:
  Stage 2 summary.*

## §5 Six-storm validation and the steering buffer
- **§5.1 The differential projection method** — the well-posed form of paper 1's projection
  test (steering errors cancel in the A/B); per-storm registered predictions, H1 vs H2
  framing, guards.
- **§5.2 Results.** Sign 6/6; transmission ratio ~0.34 (clean four); P-S1 core pass;
  strong-form H1 failed; skill accounting both ways (−8% clean / −25% aggregate);
  the Vmax guard catching Hugo/Ivan.
- **§5.3 The buffer.** What a transmission ratio of ~1/3 means: the relaxation architecture
  partially self-corrects self-propagation error; implications for landfall-based attribution
  anywhere (you cannot read drift biases from landfall errors at full strength — quantified);
  connection to paper 1's DLM double-counting question.
- *Feeds from:* Stage 3 registration + scoring tables. *~1200 words. Figure: predicted-vs-
  observed Δcross scatter (slope ≈ 1/3, all same sign) — THE figure of the validation half.
  Table: Stage 3 six-storm c→t.*

## §6 Discussion
- The bounded-vortex trade space: taper → aim error; envelope → DELAYED barotropic
  re-intensification (~8–12 h; Hugo/Ivan scored mid-delay at landfall). Mechanism measured
  (AM-budget Run 4): weaker mid-radius Ekman inflow (surface wind at 200–400 km sets BL mass
  convergence; ×2.5–3 pre-onset), NOT reservoir starvation. No free boundary; characterize
  the trade.
- Generality, sharpened from paper 1 §5.4: any model that bounds its initial vortex with
  compact support should check its gyre phase-lock; the diagnostic is an afternoon.
- What remains: Arm C-v2 (relax-to-θ′_ref); the envelope-intensification study; a second
  basin under the new profile.
- Methodology coda, calm: five registered predictions confirmed, two failed/unreadable, all
  reported; one live-caught uncommitted-config regression (the fran/michael/laura loss) —
  the paper-1 discipline, still earning its keep. *~700 words.*

## §7 Conclusions
- Mechanism, fix, buffer; the bias is gone from the testbed and bounded in the storms; the
  method is the constant. *~300 words.*

---

## Figures & tables inventory

- **Fig. 1** — wind/vorticity profiles: compact taper vs Gaussian envelope (V(r) + implied
  ζ(r); shows the ring and the amputated tail). *Build from vortex_init (trivial).*
- **Fig. 2** — heading vs time: compact (marching through north) vs envelope (locked). *Data
  exists: gate-beta-longrun/timeevol traces + shape/envelope per-window traces.*
- **Fig. 3** — predicted vs observed Δcross scatter, six storms (slope ≈ 1/3). *Data in
  Stage 3 table.*
- **Table 1** — Arm A summary (five profiles). **Table 2** — Stage 2 calibration (seven
  rows). **Table 3** — Stage 3 six-storm control→treatment with guards.
- *Candidate supplemental:* the m=1 gyre snapshot pair (locked vs precessing) if a clean
  envelope gyre figure is wanted — one gate-beta-gyre run at r_d=420 would produce it.

## Open inputs / decisions

- **Venue + title** (MWR + #1 recommended).
- **Paper 1 §5.1 coordination:** add the one-sentence forward reference ("a subsequent
  profile-family test implicated the compact-support cutoff [in prep]") or leave paper 1
  silent. Decide before paper 1 submits — the sequel's existence argues for the sentence.
- **Arm C-v2 before submission?** Candidate 2 is honestly "not adjudicated"; §2.3 can carry
  that as-is, or one relax-to-θ′_ref run closes it. Cheap; recommended but not blocking.
- **Envelope gyre snapshot** (supplemental figure) — one testbed run if wanted.
- Refs to verify: Fiorino & Elsberry (1989) — the outer-wind-control result this mechanism
  leans on; Chan & Williams (1987); paper 1 cross-refs once its section numbers freeze.
