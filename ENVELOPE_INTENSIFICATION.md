# AM-budget study — envelope intensification + the NZ=64 drag artifact

**Date:** 2026-07-03 · **Status:** harness built (`oracle_v8/measure_am_budget.py`); analytic
result derived and component-verified; predictions registered BEFORE any GPU run (do not edit
after; append results). Two birds, one instrument.

## Bird 2 first, because it is already dead on paper

**The NZ=64 vortex spin-down ([LH82 study Phase 3B open flag]) is a drag-discretization
artifact — derived analytically, then confirmed at the component level.**
`SurfaceDragComponent` applies α = Cd·|V′|/dz · max(0, 1−z/H_bl). The 1/dz prefactor is correct
for single-layer application; spread over the boundary-layer profile it over-counts. The
column-integrated sink factor S = Σ_k max(0, 1−z_k/H_bl):

| grid | levels in BL | S (component-measured) |
|---|---|---|
| NZ=32 (dz=625 m) | 312.5, 937.5 m | **0.7500** |
| NZ=64 (dz=312.5 m) | 156.25, 468.75, 781.25 m | **1.5938** |

Halving dz multiplies the integrated surface drag by **2.125** (continuum limit S → H_bl/2dz,
divergent). A doubled momentum sink flipping marginal spin-up into decay is precisely the
Phase 3B observation (64→48 m/s at NZ=64 vs intensification at NZ=32). Fix implemented:
`column_normalized=True` divides by the discrete integral of the weight profile — component-
verified sink = 1.0000 × Cd|V′| on both grids. Default remains False (bit-identical).

**Production note (independent of any run):** at NZ=32 the historical effective column drag is
0.75 × the nominal bulk value — i.e. the calibrated production configuration has effective
C_d ≈ 1.125e-3, not 1.5e-3. Documentation fix at minimum; switching production to normalized
drag is a recalibration decision, not a bug fix, and is NOT proposed here.

## Bird 1 — why does the envelope damp barotropic (re)intensification?

Stage 3's P-S4 guard: Hugo −14.2 and Ivan −25.6 m/s weaker under the envelope; the four
decay/steady storms ≤2.1. Both flagged storms have strong barotropic (re)intensification
phases; the intensify-ladder (J0) showed this spin-up needs no β or steering. Hypotheses:

- **H-A (M supply):** spin-up feeds on angular-momentum import by the drag-driven BL inflow;
  the profiles differ in relative-M at the radii the inflow draws from (compact has more at
  200–350 km; envelope has more at 400 km + the tail — NOT a one-signed difference, which is
  why we measure rather than assert).
- **H-C (inflow structure):** drag ∝ local |V′| ⇒ the profiles drive different Ekman pumping
  and secondary-circulation strength; the import differs via u_r, not via M.
- **H-D (coupling):** the phenomenon needs β/steering and won't reproduce on the f-plane
  (would contradict J0's precedent; if so, re-run rows with β on).

## The experiment (`measure_am_budget.py`, six rows, f-plane, init Vmax 64, 48 h)

A compact/NZ32 (control) · B gauss-420 · C gauss-560 · D compact/NZ64 · E compact/NZ64
normalized-drag · F compact/NZ32 normalized-drag. Instrument every 2 h: low-level Vmax,
relative-AM reservoirs (r<300/500/800 km), ring imports at 300/500 km (column + BL branch),
BL mass inflow, drag and cap torques inside 500 km, max|w|, center guard.

## Registered predictions (Claude, 2026-07-03)

- **P-EI1 (reproduction):** row A spins up on the f-plane (early dip then growth; Vmax_end ≥ 60,
  J0 precedent) and row B ends ≥ 10 m/s below row A — the storm phenomenon reproduced in
  isolation. Confidence ~65%. If A itself fails to spin up, H-D routing: repeat A/B with β.
- **P-EI2 (mechanism):** the A−B gap is accounted (≥70%) by the ring-import difference at
  300–500 km (with the BL branch carrying most of it), drag-torque difference secondary.
  Confidence ~50% — H-C is a live alternative; the instrument separates them (import split
  into M-carried vs mass-inflow-carried terms).
- **P-EI3 (tail dial):** row C lands between A and B or above B in Vmax_end — intensification
  increases with envelope scale. Confidence ~55%.
- **P-N1 (bird 2 numerical):** row D decays from 64 (Vmax_end ≤ 50, reproducing Phase 3B);
  row E recovers to within ~5 m/s of row F; the D−E gap ≥ 10 m/s. Confidence ~75%.
- **P-N2 (the 0.75→1.00 effect):** row F ends 3–10 m/s below row A (33% more column drag at
  the production grid ⇒ lower equilibrium). Whatever F−A is, it quantifies the hidden
  effective-Cd calibration noted above.

**Decision rules:** P-EI1 pass ⇒ read P-EI2's budget attribution and classify H-A vs H-C;
P-EI1 fail ⇒ β-coupled, redesign. P-N1 pass ⇒ NZ=64 chip CLOSED (artifact found, fixed,
LH82-study caveat 3 updated); P-N1 fail ⇒ the drag factor is real but not sufficient — look
at vertical advection/diffusion of momentum at dz/2 next.

## Cost & run

```powershell
python -m oracle_v8.measure_am_budget          # all six rows, ~2.5 h GPU
# $env:AMB_ROWS="AB" etc. for subsets; $env:AMB_HOURS to shorten
```

## Results

### Run 1 (2026-07-03): NOT SCORED — instrument fault, caught by the instrument's own books

All six rows integrated stably for 48 h, but the harness's headline metric read
max speed at **k=0** — the drag-drained surface level — and reported Vmax_end 2.9–5.4
everywhere. The reservoir column contradicted it: res500 within ~10–20% of initialization in
every row — the circulations were alive; the metric was measuring the drag layer's local
equilibrium (τ_drag ≈ 23 h at 5 m/s), which decouples in a quiescent testbed because this model
has no vertical momentum diffusion. Every prior intensity number (J0's 45→76, gate-beta's 42,
the LH82 study's 64→48) used different instruments (`low_level_vmax`, max z<3 km; or max|u|).
**The registered predictions are NOT scored against this run** — they remain frozen; the
instrument was invalid for them. Qualitative note, explicitly not banked: even the surface
metric showed the bird-2 ordering (D 2.9 < E 4.0 ≈ F 4.2 < A 5.2, and E ≈ F is the fix's
signature) — suggestive, unusable at that compression.

**Harness fix (same day):** headline metric replaced with the production instrument
(`low_level_vmax`), with three instruments now reported side by side (production z<3 km max;
k=0 surface; max|u| for LH82-study comparability), a BL/above split of res500, and a final
Vmax(z) profile per row (the decoupling picture). Re-run required; predictions unchanged.

### Run 2 (2026-07-03, GPU, fixed instrument): SCORED — most predictions failed, informatively

| row | Vmax_end (z<3km) | v_sfc | max\|u\| | res500 |
|---|---|---|---|---|
| A compact, nz=32 | 39.5 | 5.2 | 39.3 | 2.204e22 |
| B gauss-420 | 37.4 | 5.4 | 37.3 | 2.454e22 |
| C gauss-560 | 38.2 | 5.3 | 38.1 | 3.250e22 |
| D compact, nz=64 | 41.4 | **2.9** | 41.0 | 2.174e22 |
| E compact, nz=64, normalized | 41.2 | **4.0** | 40.9 | 2.193e22 |
| F compact, nz=32, normalized | 39.4 | **4.2** | 39.3 | 2.191e22 |

**Scorecard:** P-EI1 **FAILED** (A did not spin up — 39.5, not ≥60; the J0 precedent used
Ivan-specific Rmax/B overrides; the standard testbed vortex decays to ~38–41 under every
profile, matching gate-beta's 42). P-EI2 **unscoreable** (no gap to attribute: A−B = 2.1 m/s,
the same small equilibrium offset Stage 2 found, NOT the storm effect). P-EI3 trivially
satisfied on a 2-m/s spread — uninformative. P-N1/P-N2 **FAILED in the vortex metric** (D ≈ E ≈
A at low levels) — but **CONFIRMED at the surface**: v_sfc orders exactly as the drag analysis
demands (hist-nz64 2.9 « norm-nz64 4.0 ≈ norm-nz32 4.2 < hist-nz32 5.2; E≈F is the fix's
signature, now demonstrated in the full numerics, not just the component test).

**The finding Run 2 actually delivered: surface drag couples to vortex intensity only through
a driven secondary circulation.** In a decaying quiescent vortex the boundary layer drains
locally and the flow above does not care — the doubled column drag at nz=64 is real and
measurable at the surface yet invisible in low-level Vmax. Both target phenomena live in
INTENSIFYING regimes: the NZ=64 spin-down was observed in the heated buoyancy-on LH82
configuration; Hugo/Ivan's envelope damping occurred during β+steering re-intensification.
Also banked: the envelopes carry MORE total angular momentum than the compact profile
(res500 ordering C > B > A) and intensity does not follow the reservoir — the naive
trimmed-reservoir hypothesis (H-B flavor of H-A) is dead independently.

**Registered-confidence audit:** P-EI1 was given 65%, P-N1 75% — both wrong. The common error:
anchoring on intensity outcomes (J0's spin-up, Phase 3B's decay) without matching the REGIME
that produced them.

### Redesigned experiments (Run 3, to be registered before running)

1. **Bird 2, properly aimed:** the LH82 Phase-3B configuration itself (Q=1e-2, buoyancy ON,
   upwind5h, NZ=64, dt=15) with historical vs column_normalized drag, plus the NZ=32 reference.
   If normalized drag converts the 64→48 spin-down back toward the NZ=32 behavior (74–84), the
   LH82 caveat-3 flag closes. ~3 runs at 128², minutes each on GPU.
2. **Bird 1, properly aimed:** the J2-style intensifying rung (β + steering — the configuration
   that demonstrably intensifies in this model) with compact vs gauss-420, AM budget riding
   along. If the envelope damps THAT intensification, the budget says which term carries it.

## Run 3 — registered predictions (2026-07-04, harnesses built; frozen BEFORE any GPU run)

**Bird 2 — `measure_nz64_drag.py`** (heated Phase-3B regime: Q=1e-2, buoyancy ON, upwind5h,
ε=0, dt=15; 2×2 = {nz 32/64} × {drag historical/normalized}; metric = max|u|, the Phase-3B
instrument):
- **P-R3N1:** nz=64 historical reproduces the spin-down (max|u| ≤ 55; recorded 48.1). ~90%
  (deterministic re-run through a default-bit-identical drag refactor).
- **P-R3N2 (the kill shot):** nz=64 normalized recovers intensification — max|u| ≥ 65 AND
  ≥ 15 m/s above nz=64 historical. Confidence ~60%: the regime is now right (heated BL inflow
  is the pathway a doubled BL sink strangles), but dz/2 also changes vertical advection and
  heating-layer resolution, which may carry part of the recorded gap.
- **P-R3N3:** nz=32 normalized lands 3–12 m/s below nz=32 historical (the 0.75→1.00 column-drag
  effect, which the heated regime SHOULD express, unlike Run 2's decay regime). ~65%.
- **Decision:** P-R3N2 pass ⇒ the LH82 caveat-3 NZ=64 flag CLOSES (mechanism + validated fix);
  fail ⇒ drag factor real but insufficient — remaining suspects (vertical advection,
  heating-layer resolution) inherit the flag, better localized.

**Bird 1 — `run_translation_test.py gate-j2-profile`** (J2 rung: Ivan structure Vmax 66.9
B 1.5 Rmax 75 km, β + steering ramp (−1.1,3.9)→(+0.6,6.6), cap 70, 320²/5000 km, 52 h;
compact taper vs frozen production envelope r_d=420):
- **P-J21:** compact J2 re-intensifies, peak ≥ 75 (ladder reference ~80). ~80%.
- **P-J22:** envelope J2 peaks AND ends ≥ 10 m/s below compact — the Hugo/Ivan damping
  reproduced in isolation. ~55% (the storms say yes; J2's ramp is not the storms'
  position-dependent ERA5 sampling — a null here would localize the damping to the
  steering-feedback pathway, informative either way).
- **P-J23 (conditional):** if P-J22 passes, the mechanism is NOT total-reservoir starvation
  (Run 2: envelope carries MORE total AM) — the follow-up is a moving-frame AM budget to split
  BL-import vs inflow-strength. No numeric prediction registered.
- **Guards:** final y vs the β-taper zone (~4000 km); early-time (t ≤ 24 h) Vmax comparison
  before the profiles' tracks diverge.

### Run 3 results (appended 2026-07-04, GPU runs by Justin; scored vs the frozen registrations)

**Bird 2 (`measure_nz64_drag`, heated Phase-3B regime):**

| row | max\|u\| |
|---|---|
| nz=32 historical | 84.3 (recorded ref 84.3 — exact) |
| nz=32 NORMALIZED | 69.6 |
| nz=64 historical | 48.1 (recorded ref 48.1 — exact) |
| nz=64 NORMALIZED | 54.7 |

**P-R3N1 PASS exactly** (both historical rows bit-reproduce the Phase-3B records — also the
strongest available regression proof that the drag refactor's default is bit-identical).
**P-R3N2 FAIL** (54.7, not ≥65; +6.6 of the 36.2 gap). **P-R3N3 FAIL by overshoot** (−14.7,
band was 3–12; the heated regime expresses the column-drag factor at double the predicted
strength — vs zero in Run 2's decay regime).

**Verdict — the flag splits instead of closing.** Drag-matched (normalized) cross-grid gap:
69.6 − 54.7 = **14.9 m/s of genuine, drag-independent dz-sensitivity** in heated
intensification. So: (a) the drag discretization artifact is real, in-regime large, quantified,
and fixed (`column_normalized`); (b) it is NOT the dominant cause of the NZ=64 spin-down;
(c) the residual belongs to dz/2 itself — vertical advection and heating-layer resolution
inherit a better-localized flag. Production note stands: any future buoyancy-on production
config must decide its drag normalization deliberately (the heated regime is strongly
sensitive; the barotropic track config is not).

**Bird 1 (`gate-j2-profile`):**

| profile | init | peak | end | y_end |
|---|---|---|---|---|
| compact (J2 reference) | 63 | 82 | 76 | 3809 |
| gauss r_d=420 | 61 | 77 | 73 | 3745 |

**P-J21 PASS** (peak 82; J2 reproduced with the decay-then-reintensify shape). **P-J22 FAIL as
registered** (peak gap +5.6, end gap +2.8) — and the registered metric asked the wrong
question. Time-resolved: compact re-intensifies at t≈24–28; the envelope holds 45–49 until
t≈36–40, then intensifies and converges. Matched-time gaps: **+30 m/s (t32), +28 (t36),
+2.8 (t52)**. **The envelope does not suppress the barotropic re-intensification — it
POSTPONES it ~10–12 h and then delivers nearly the same storm.** This resolves the Stage-3
Hugo/Ivan "damping": their runs end at landfall, inside the delay window — a delayed
intensifier scored at a fixed clock reads as a weakened one. Untested hypothesis for the
delay, queued for a moving-frame AM budget: reservoir LOCATION (compact holds its relative AM
at 200–350 km where BL inflow reaches quickly; the envelope's sits at 400–900 km — a longer
supply chain, same cap-limited equilibrium). Guards: both tracks stayed inside the β-taper
interior (3809/3745 < ~4000, margin thin — flag for any longer rerun); early-time gaps ≤4 m/s.

**Confidence audit (Run 3):** P-R3N1 90% ✓, P-J21 80% ✓, P-R3N2 60% ✗, P-R3N3 65% ✗
(direction right, band half the true size), P-J22 55% ✗ (wrong metric — peak/end instead of
onset time). Recurring lesson now twice-paid: register predictions about TRAJECTORIES, not
endpoints, when the phenomenon is a feedback with a threshold.

### Study status after Run 3

- **Envelope-intensification (bird 1): mechanism-class RESOLVED at the phenomenon level** —
  delay, not suppression; paper-2 §6's trade-space paragraph should say "delays dry barotropic
  re-intensification (~10 h in the J2 testbed)" rather than "damps." Optional next: the
  moving-frame budget to pin the supply-chain hypothesis.
- **NZ=64 (bird 2): drag artifact found/fixed/quantified; residual 14.9 m/s dz-sensitivity
  remains open** (vertical advection / heating-layer resolution) — logged in
  LH82_SMALL_PERTURBATION_FINDINGS.md caveat 3.

## Run 4 — the moving-frame budget (registered 2026-07-07, BEFORE any GPU run)

**Question:** what sets the envelope's ~10–12 h re-intensification delay? Hypothesis on the
table (from Run 3): reservoir LOCATION — the compact profile's relative AM sits at 200–350 km
where the BL inflow reaches it quickly; the envelope's sits at 400–900 km, a longer supply
chain to the same cap-limited equilibrium. Alternative: inflow STRENGTH (the profiles drive
different Ekman pumping; the difference is in how much air flows in, not what it carries).

**Design:** `measure_am_budget.py` rows G/H (`$env:AMB_ROWS="GH"`) — the J2 rung exactly as
Run 3 ran it (β + ramp, Ivan structure, 52 h, 320²), compact vs frozen envelope-420, with the
AM budget in the moving frame (vorticity-center tracked, geometry recentered, perturbation-wind
fields so the steering cannot leak into the rings via gyre asymmetries). New summary metric:
onset(>60 m/s, post-dip) — a trajectory quantity, per the twice-paid Run-2/3 lesson.

**Registered predictions (Claude, 2026-07-07):**
- **P-M1 (integrity):** rows G/H reproduce Run 3's trajectories (onsets ≈ 26–30 h and
  ≈ 38–42 h; peaks within ~3 m/s of 82/77) — the budget is read-only and must not perturb the
  runs. ~85%.
- **P-M2 (the hypothesis):** in each row, the rise of impBL300 (BL relative-M import at 300 km)
  LEADS that row's own Vmax onset by ≤ 6 h; and in the pre-onset window (t = 12–24 h) the
  compact row's impBL300 exceeds the envelope's by ≥ 2×. ~50% — honest coin-flip; that is why
  we measure.
- **P-M3 (the discriminator):** in the same window the BL mass inflow (minBL300) is comparable
  between rows (within ~30%) while impBL300 differs — the difference is in what the inflow
  CARRIES (reservoir location), not how much flows (inflow strength). If instead minBL300
  differs ≥ 2×, the mechanism is inflow strength. ~45% on the location branch.
- **Decision:** P-M2+P-M3 pass ⇒ supply-chain mechanism confirmed; paper-2 §6 gets one
  mechanism sentence with a measured lead time. P-M2 passes but P-M3 fails toward
  inflow-strength ⇒ equally clean, different sentence. Both fail ⇒ the delay is not
  BL-import-controlled — look at the gyre/asymmetry pathway next (the β-gyres themselves
  redistribute AM), and say so honestly.

**Cost:** 2 × 52 h × 320² + budget ≈ 25–35 min GPU
(`$env:AMB_ROWS="GH" ; python -m oracle_v8.measure_am_budget`).

### Run 4 results (appended 2026-07-07, GPU run by Justin; log `AM-budget_gh.txt`)

| row | peak | end | onset(>60) | y_end |
|---|---|---|---|---|
| G compact | 80.7 | 80.4 | 24 h | 3758 |
| H gauss-420 | 78.6 | 73.2 | 32 h | 3867 |

**The delay reproduced (8 h) with near-identical peaks — measured twice now.**

**Scorecard:**
- **P-M1 PARTIAL:** peaks/ends reproduce Run 3 within ~3 m/s; absolute onsets ran 2–6 h
  earlier than the registered bands (envelope 32 vs 38–42). Two config seeds are mine — the
  harness hardcoded F = 5.7e-5 and V0 = 66.9 where the ladder uses IVAN's exact 5.6985e-5 and
  66.878 (~0.03% each); near a threshold, deterministic chaos amplifies seeds into hours. The
  A/B *within* Run 4 shares one config exactly, so the 8-h delay is clean.
- **P-M2 FAIL, instructively:** pre-onset, the BL relative-M flux at 300 km is a net EXPORT in
  both rows (gyre-asymmetric drainage), flipping to strong import only AT each row's onset
  (coincident within the 2-h cadence). The M-import does not lead intensification — it is
  intensification, viewed in the books.
- **P-M3 RESOLVED — INFLOW-STRENGTH branch:** pre-onset (t = 12–24 h) BL mass inflow
  (minBL300) is **2.5–3× stronger in the compact row** (t20: 4.2e7 vs 1.5e7; t24: 4.8e7 vs
  1.8e7), rising steadily from t≈8 (compact) vs t≈12 (envelope). NOT comparable-within-30%;
  the reservoir-location hypothesis is rejected (consistent with Run 2's total-AM finding).

**MECHANISM (bird 1, closed at phenomenon+mechanism level):** the envelope's re-intensification
delay traces to weaker drag-driven Ekman pumping at mid-radii — the surface wind the profile
places at 200–400 km sets the BL mass convergence (compact ~45 m/s at 200 km vs envelope ~36 →
×2.5–3 inflow), so the core spins up later; once the intensification feedback closes, the same
cap-limited equilibrium is reached. **Paper-2 §6 sentence: "the envelope delays dry barotropic
re-intensification (~8–12 h) by weakening mid-radius Ekman inflow, not by starving the
angular-momentum reservoir."**

**Instrument caveat:** row G's vorticity-center jumps ±200 km post-onset (t26/38/52) with
budget sign-flips riding along — the known ζ²-centroid fragility at high intensity
(gate-beta-longrun note). Pre-onset windows (where all scoring above lives) are clean in both
rows; post-onset G budget lines are not load-bearing.

**Confidence audit:** P-M1 85% → partial; P-M2 50% → fail (export, not import); P-M3 ~45% on
location → strength branch confirmed. The discriminator design (P-M3) is the win: it was built
to decide, and it decided.

## Study status — CLOSED (2026-07-07)

Bird 1: envelope delay = weaker mid-radius Ekman inflow (measured, ×2.5–3 pre-onset). Bird 2:
drag discretization found/fixed/quantified; residual ~15 m/s dz-sensitivity of heated
intensification remains the one open flag (vertical advection / heating-layer suspects).
Four runs, twelve registered predictions, five confirmed, six failed-and-reported, one
unreadable — every failure narrowed the search. Follow-ups live in the resumption menu
(am-budget-study memory).
