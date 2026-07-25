# Reviewer artifact requests — pointers and honest answers

*Response package, 2026-07-24. Four artifacts were requested; each is a repo path, with an
honest note where the review's suspicion was correct.*

## 1. Taper/envelope construction and the Γ(r) computation

`oracle_v8/vortex_init.py` — `HollandVortexInit.tangential_wind` contains all outer boundary
conditions in one method: the compact-support ramps (`wind_taper`, shapes cos/linear/smooth5),
the Gaussian envelope (`outer_envelope_m`), and, added this commit, the ring+tail
discriminator (`ring_plus_tail`). The implied-vorticity panel of Paper 2 Fig. 1 is computed in
`figures/make_paper2_figures.py` (ζ = (1/r) d(rV)/dr on the profile grid).

**On the ring-vs-cutoff confound: the reviewer is right.** Compact support fixes the removed
circulation, so the three Arm A ramps share integrated ring strength and differ only in edge
sharpness; both hypotheses predicted Arm A's outcome. We checked for an existing
discriminating run: there is none — `gate-beta-renv` (the closest) moves the ring outward but
confounds it with circulation size. The clean single-variable test the review specified is
registered (`PAPER2_RING_predictions.md`, predictions frozen pre-run) and implemented as
`gate-beta-ringtail`; Paper 2 §2.1 now states the design limitation explicitly.

## 2. Testbed drift-vector reader (resolution sweep)

`oracle_v8/run_translation_test.py` — `_mature_drift` (the 30–48-h window vector; the signed
circular-delta heading readout is the post-wraparound-bug version — the bug and its catch are
documented in `GATE_BETA_RES_SWEEP.md` and recounted in Paper 1 §3.3). The resolution sweep
harness and thresholds: `GATE_BETA_RES_SWEEP.md`; the run log: `GATE_BETA_RES.txt` series
(320/480/640 → 350°/351°/352°, 2.49/2.39/2.29 m s⁻¹).

**On convergence: the reviewer is right.** Successive magnitude differences are 0.10 and
0.10 m s⁻¹ (ratio 1.00) — no asymptotic regime. Both manuscripts now say so: Paper 1 §4.1
claims convergence only for the *heading* and flags the magnitude's undemonstrated resolution
dependence; Paper 2 §5.1 states that Δ and every ratio built on it inherit that caveat.

## 3. Transmission ratios — "are the aggregate and the mean computed in the same place?"

**They were not.** The guard-clean mean lived in the Stage-3 A/B append
(`OVERROTATION_CANDIDATES.md`), the "0.42 product" in the AM-budget accounting, and Fig. 3's
line in the figure script — and the manuscript conflated populations across them (the §5.3
validation compared the guard-clean prediction mean to the all-six aggregate). All variants
are now computed in one place, `oracle_v8/measure_transmission_summary.py`:

```
guard-clean mean of ratios (exported)      0.340
guard-clean aggregate                      0.352
all-six aggregate (incl. flagged storms)   0.415
through-origin fit (Fig. 3)                0.365
decomposition prediction mean              0.445  -> honest validation: 0.45 vs 0.34
nine-storm guard-clean mean                0.350
```

The manuscripts now export the guard-clean mean (0.34), state the population rule (a
guard-excluded storm never enters a validation target), and report the honest 0.44-vs-0.34
overshoot with its driver (Michael). The registered P-D4 verdict stands as scored; an
amendment in `PAPER2_TIGHTENING_predictions.md` records the mis-specification.

## 4. Per-storm decomposition and expansion scoring

`oracle_v8/measure_transmission_decomp.py` (six-storm, registered P-D1..P-D4),
`oracle_v8/measure_expansion_score.py` (Charley/Florence/Ida, registered P-N1..P-N5), with
run logs under `oracle_v8/Logs/{Storm}/` and every registration/outcome in
`PAPER2_TIGHTENING_predictions.md` and `PAPER2_EXPANSION_predictions.md`.
