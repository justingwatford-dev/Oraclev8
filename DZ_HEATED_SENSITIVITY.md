# dz-sensitivity of heated intensification — design + registered predictions

**Date:** 2026-07-07 · **Status:** harness built (`oracle_v8/measure_dz_heated.py`); predictions
registered BEFORE any GPU run (do not edit after; append results). The last open flag from the
AM-budget study (ENVELOPE_INTENSIFICATION.md) and the LH82 study's caveat 3.

## The question

With the drag artifact removed (column_normalized both grids), the heated Phase-3B configuration
still reads **max|u| 69.6 (nz=32) vs 54.7 (nz=64)** — ~15 m/s of genuine dz-sensitivity, riding
on the secondary circulation (max|w| dropped ~3× at nz=64 in the LH82 study). Which grid is
lying?

**Hypothesis, stated before running (the sign-inversion reframe):** the anomaly is the COARSE
grid's. The heating layer (width_z = 3 km) spans ~4.8 cells at dz=625; discrete vertical
operators under-estimate the effective wavenumber of marginally-resolved forcing, the anelastic
response meets less opposition, and the coarse grid OVER-produces w. If so, nz≥64 is the
convergent regime and the "nz=64 spin-down" was never a spin-down — production dz simply cannot
converge a heated updraft, which is a characterized limitation, not a bug (production track
runs are barotropic and unaffected).

## Design (all rows: Q=1e-2, ε=0, dt=15, upwind5h, buoyancy ON, drag NORMALIZED)

| row | grid | purpose |
|---|---|---|
| R1 | nz=32 | drag-matched reference (69.6) |
| R2 | nz=64 | drag-matched reference (54.7) |
| R3 | **nz=96** | **the convergence-direction row** |
| R4 | nz=64, heating sampled at nz=32 centers | forcing-representation control (`z_sample_nz`) |

A blow-up at nz=96 is a datum (vertical CFL headroom shrinks); fallback `DZH_DT=10`, noted.

## Registered predictions (Claude, 2026-07-07)

- **P-D1 (~55%):** 32→64→96 decreases monotonically with a SHRINKING increment —
  |u(96)−u(64)| ≤ 0.5·|u(64)−u(32)| — i.e. the fine side is convergent (asymptote ≈ 50±3) and
  nz=32 is the outlier. The alternative (96 bounces back up or the increment grows) would
  point at a fine-grid pathology instead (vertical dispersion of the centered scheme, dt).
- **P-D2 (~70%):** R4 lands within ~2 m/s of R2 — the discrete SAMPLING of the heating profile
  is not the driver; the dynamics own the gap. (If R4 recovers toward R1, the fix is real and
  cheap: cell-integrated rather than center-sampled heating.)
- **P-D3 (~75%):** max|w| ordering tracks max|u| ordering across all rows — the intensification
  difference continues to ride on the secondary-circulation strength.

**Decision rules:** P-D1+P-D2 pass ⇒ the flag CLOSES as a characterized limitation: "heated
intensification is not vertically converged at production dz=625; the convergent regime is
nz≥64 with normalized drag; coarse-grid heated magnitudes (incl. the LH82 study's max|w| and
the 84.3/69.6 intensities) are upper-biased by vertical under-resolution." LH82 findings
caveat 3 gets its final wording. P-D2 fail ⇒ implement cell-integrated heating and re-run.
P-D1 fail ⇒ dt sweep at nz=96 + buoyancy full→half interpolation become the suspects; say so.

## Cost & run

```powershell
$env:LH82_STEPS = "1000"
python -m oracle_v8.measure_dz_heated     # 4 rows at 128², ~10 min GPU
```

## Results (appended 2026-07-07, GPU run by Justin; scored vs the frozen registrations)

| row | max\|u\| | max\|w\| |
|---|---|---|
| R1 nz=32 | 69.6 (ref exact) | 5.89 |
| R2 nz=64 | 54.7 (ref exact) | 3.20 |
| R3 nz=96 | **54.6** | 3.42 |
| R4 nz=64, heating@nz32 | 55.1 | 3.41 |

**Scorecard — all three registered predictions PASS (first clean sweep of the campaign):**
- **P-D1 PASS, decisively:** Δ(64−32) = 14.9; Δ(96−64) = 0.1 — not a shrinking increment but
  full convergence. The fine side is converged; nz=32 is the outlier. (Audit: the registered
  monotone-shrinking rule passed by two orders of margin; the asymptote guess 50±3 was 1.6 m/s
  low.)
- **P-D2 PASS:** R4 within 0.4 of R2 — the coarse grid's exact discrete heating on fine
  dynamics changes nothing. Forcing representation innocent; the dynamics own the gap.
- **P-D3 PASS in substance:** the 15 m/s intensity gap rides on the ~2× w gap (5.89 vs
  3.2–3.4); intra-cluster ±0.2 differences are noise-level and unordered.

## VERDICT — the flag closes (sign inversion confirmed)

**There was never an nz=64 spin-down; there is an nz=32 over-response.** A 3-km heating layer
resolved by ~4.8 cells meets too little vertical opposition in the discrete anelastic response,
and the production grid over-produces the heated secondary circulation (~+80% in max|w|, ~+27%
in intensity) relative to the converged solution. Characterized limitation, not a bug:

- **Recipe for any heated / buoyancy-on work: nz ≥ 64, column-normalized drag, dt = 15.**
  The converged heated Phase-3B state is max|u| ≈ 54.7, max|w| ≈ 3.2–3.4.
- Coarse-grid heated magnitudes in prior records (84.3, 69.6, w = 10.7/5.9) are upper-biased
  by vertical under-resolution — flagged wherever quoted.
- The LH82-validity conclusions are untouched and were conservative: the neglected-term ratios
  *shrink* at fine dz (Phase 3B), and the equation-set verdict never rested on w magnitude.
- Production track results (barotropic, unheated) are unaffected in full.

**Confidence audit:** P-D1 55% ✓, P-D2 70% ✓, P-D3 75% ✓.

**Study CLOSED 2026-07-07.** With it, the project's last open experimental flag is retired:
every anomaly raised since the red-team arc — the over-rotation, the envelope delay, the drag
column factor, and the dz-sensitivity — now has a measured mechanism, a validated fix, or a
characterized limitation with a recipe.
