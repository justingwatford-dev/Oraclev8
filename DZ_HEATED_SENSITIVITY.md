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

## Results

*(append after the GPU runs; predictions frozen)*
