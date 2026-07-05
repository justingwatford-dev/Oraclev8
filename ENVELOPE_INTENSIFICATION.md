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

### Run 3 results

*(append after the GPU runs; predictions frozen)*
