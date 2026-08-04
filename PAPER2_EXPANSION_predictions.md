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

---

# OUTCOMES (scored 2026-07-21; runs on the GPU backend — woe_env cupy/RTX 5070, both arms
# per storm on the same backend so the A/B is backend-consistent; scoring script
# `oracle_v8/measure_expansion_score.py`; logs `oracle_v8/Logs/{Charley,Florence,Ida}/`)

**Tally: 3 confirmed / 2 failed — and the P-N4 failure is the campaign's most valuable
result.** All guards clean on all three storms (max A/B intensity divergence 2.6 m s⁻¹,
A/B timing ≤ 0.8 h): every ratio below counts.

| storm | cross c→e (km) | along c→e (km) | axis | obs shift | strong | obs ratio | decomp pred |
|---|---:|---:|---|---:|---:|---:|---:|
| Charley | −31.8 → −45.0 | −160.7 → −172.8 | cross | −13.2 | −59 | **0.22** | 0.42 |
| Florence | +62.7 → +65.7 | −124.3 → −78.3 | along | +46.0 | +73 | **0.63** | 0.56 |
| Ida | +36.1 → +13.1 | −14.0 → +0.6 | cross | −23.0 | −97 | **0.24** | 0.31 |

- **P-N1 CONFIRMED** (3/3; every |cross| ≤ 150 both arms). The lean was wrong in the best
  way: Charley's placement held (−32/−45 km cross) — its damage is along-track (−161/−173 km,
  +4.4/+5.2 h late), the dry capped model unable to follow the observed rapid acceleration
  and RI into landfall, compounded by the registered Cuba crossing. Ida's envelope landfall:
  **13.1 km total error** (+2.6 km same-latitude) — the best landfall in the nine-storm record.
- **P-N2 CONFIRMED — nine for nine on sign.** Charley west, Florence along-forward, Ida west.
- **P-N3 FAILED, informatively** — Florence's along-axis ratio 0.63 exceeds the [0.15, 0.55]
  band. The band was generalized from cross-axis experience; a long, slow transit spends more
  of its time at mature drift, and the along axis evidently transmits harder. The band was
  the wrong generalization, not the storm.
- **P-N4 FAILED — 2/3 closed, and the miss is the discovery.** Florence within 0.07, Ida
  within 0.07 — the §5.3 mechanism survives its first true out-of-sample test on both
  quasi-straight/zonal movers. Charley misses by 0.20 (0.42 predicted, 0.22 observed) —
  **the Michael shortfall (0.45 predicted, 0.20 observed) replicated, under frozen
  predictions, on the record's only other sharply-recurving NNE landfall.** The residual is
  a reproducible class (n = 2), not a one-off. Post-hoc note (labeled): Charley's along-axis
  ratio is also low (≈0.19), so its attenuation is uniform rather than cross-specific —
  the leading untested candidate is now the *fixed-heading projection itself*: for a track
  that curves sharply into landfall, projecting the accumulated correction through the final
  heading overstates what a transit-integrated correction can deliver. Michael and Charley
  are the two curved-landfall geometries; every quasi-straight mover closes.
- **P-N5 CONFIRMED, decisively** — Florence |along| 46.0 km vs |cross| 3.0 km: the zonal
  mover reads the westward correction as timing, the transit-attenuation story's geometric
  corollary, predicted and observed.

**Running program tally: 48 registered predictions (32 campaign + 11 hardening + 5
expansion): 24 confirmed, 17 failed, 7 partial/unscoreable/against-lean.**
