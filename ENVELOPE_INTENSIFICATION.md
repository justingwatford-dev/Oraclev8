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

*(append after the GPU runs; predictions above frozen)*
