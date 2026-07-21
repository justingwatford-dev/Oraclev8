# A/B expansion — three new storms, registered predictions

*Registered 2026-07-18, BEFORE any ERA5 download or model run. Extends the six-storm A/B of
PAPER2 with three storms chosen for geometric stress (paper 1 §5.5's own list: fast movers,
slow movers, new landfall headings), under the campaign's rules: predictions and confidences
frozen now; scoring uses the shared verification module's landfall-fix along/cross
decomposition; no threshold moves after data. The decisive addition over the original six:
the per-storm transmission decomposition (measure_transmission_decomp.py) now exists, so each
new storm is an **out-of-sample test of the §5.3 mechanism**, not merely another sample.*

## The storms (configs added to `era5_steering.py` STORM_CONFIGS this commit)

| storm | id | init | transit T | landfall | est. heading | why chosen |
|---|---|---|---|---|---|---|
| Charley 2004 | AL032004 | 2004-08-12 18Z | 25.75 h | Cayo Costa FL, 13/1945Z | ~015° (NNE) | fast; lands on Michael's axis; compact storm (tests Rmax=75 km floor); **Cuba crossing at t≈10.5 h is an unrepresented land feature, registered here, not discovered later** |
| Florence 2018 | AL062018 | 2018-09-13 00Z | 35.25 h | Wrightsville Beach NC, 14/1115Z | ~290° (W, zonal) | slow, decelerating; the record's first zonal-mover geometry — the westward drift correction should read as **along-track**, inverting the usual axis |
| Ida 2021 | AL092021 | 2021-08-28 12Z | 28.92 h | Port Fourchon LA, 29/1655Z | ~330° (NW→N) | near-twin of the Katrina/Laura falsifier pair; near-pure cross-track geometry |

## Strong-form projections (mature Δ = (−0.78, −0.49) m s⁻¹ through estimated headings)

| storm | strong-form Δcross | strong-form Δalong | geometry note |
|---|---:|---:|---|
| Charley | −58 km | −63 km | mixed-axis |
| Florence | −92 km | **+72 km** | along-dominated (ahead); cross is southward |
| Ida | −96 km | −4 km | cross-dominated, near-clean |

Headings/transits above are registration-time estimates from the best track; scoring uses the
verification module's decomposition and each run's actual geometry.

## Registered predictions

- **P-N1 (70%) — blind-skill extension.** All three control (compact-taper) runs are
  scoreable, and at least two of three land with |cross-track| ≤ 150 km at the landfall fix.
  Stated lean: Charley is the likeliest casualty (Cuba + compact core), Ida the cleanest.
- **P-N2 (80%) — sign, six-for-six becomes nine-for-nine (dominant axis).** Under the
  envelope, every storm shifts in the direction of its projected correction on its dominant
  axis: Charley cross-track westward, Florence **along-track backward** (the +72 km ahead
  projection means envelope-minus-control moves toward "ahead"; observed shift should carry
  that sign), Ida cross-track westward.
- **P-N3 (55%) — transmission band.** For guard-clean storms, the dominant-axis transmission
  ratio (observed shift / strong-form) falls in [0.15, 0.55], consistent with the six-storm
  0.20–0.46.
- **P-N4 (50%) — the mechanism, out of sample.** The per-storm decomposition, fed each run's
  actual intensity histories after the fact, lands within ±0.15 of the observed dominant-axis
  ratio for every guard-clean storm. This is the §5.3 mechanism's first true out-of-sample
  test; a Michael-style axis-specific miss on Charley (the other NNE-heading storm) would be
  evidence that the residual is a class, not a one-off — informative either way.
- **P-N5 (60%) — Florence's geometry.** |observed along-track shift| > |observed cross-track
  shift| for Florence: for a zonal mover the westward correction reads as timing, which is the
  transit-attenuation story's geometric corollary.
- **Guards (registered as before):** intensity histories within ±10 m s⁻¹ between A and B;
  timing within ±3 h on direct landfalls. Guard-violating storms are flagged and excluded
  from transmission estimates, as in the original six.

## Execution order (after this commit)

1. Justin: `python -m oracle_v8.era5_steering --download --storm charley` (and florence, ida)
   — needs CDS credentials; ~10 MB each.
2. Control + envelope runs per storm (`run_storm` by name; envelope via the frozen production
   profile, control via `ORACLE_OUTER_ENVELOPE_M=taper`), logs checked in.
3. Score P-N1..P-N5, append outcomes here, fold results into PAPER2 §5 and the ledger.
