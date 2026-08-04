# Red-team audit — structural init and six-storm pivot

Status: 2026-06-23. This file maps the external red-team critique against the
actual repository state. It is intentionally conservative: if a claim is plausible
but untested, it becomes a hold item rather than a dismissed concern.

## Bottom line

The red team found real blockers in the interpretation layer. Do not run or
interpret the structural-init treatment yet. The immediate state is:

- The six-storm result supersedes the three-storm eastward headline.
- The structural recipe's storm-specific radii and derived `R_env / taper / B`
  values are useful as a data sheet.
- Directional landfall predictions are suspended until a component-decomposed
  size sweep and an ERA5 DLM sensitivity test are complete.
- Any future structural-init implementation needs a pipeline-regression mode
  that reproduces the current control constants through the new pathway.

## Findings ledger

| # | Verdict | Repository check | Action |
|---|---------|------------------|--------|
| 1 | **Valid blocker.** | `CHEATSHEET.md` mixes R_env/taper-start magnitude-and-heading results with f-plane/Vmax decompositions where the west component stays small. The current recipe predicted signs without a storm-relevant N/W component sweep. | `STRUCTURAL_INIT_RECIPE.md` now suspends predictions pending a fixed-Vmax, component-decomposed size sweep. |
| 2 | **Partly valid blocker.** | Red team says DLM is sampled at/near the storm center. Actual code: `ERA5Steering.get_dlm()` uses a 3–7° annulus, not a point sample. Still, the annulus is not a storm-removed environmental field. | Added ERA5 DLM sensitivity as a hold item. At minimum, bracket Ivan with wider-annulus or storm-removed/lagged steering. |
| 3 | **Valid.** | Five of six recipe fits hit `B=2.50`; Fran is the only non-clamped B. | Recipe now frames this as a data-informed outer-size experiment, not a clean shape experiment. |
| 4 | **Valid, already evidenced.** | `SIX_STORM_VALIDATION_NOTE.md` explicitly says the systematic eastward landfall claim does not generalize. `PAPER_track_error_characterization.md` still had an internal overclaim. | Paper paragraph patched; cheatsheet headline patched. |
| 5 | **Valid.** | Prediction bands were broad and not falsifiable. | Predictions suspended; recipe now requires one-sided pass/fail bounds after hold tests. |
| 6 | **Valid.** | Ivan and Fran hit the 800 km `R_env` clamp. | Recipe now treats those as lower-bound structural responses and keeps predictions suspended. |
| 7 | **Partly valid.** | The original "clean comparison" wording overstated the treatment; 75 km Rmax and B clamps mean both vortices are approximations. | Recipe now calls this "data-informed outer-size vortex vs generic constant-size vortex" and moves any self-consistent size-only arm to optional later work. |
| 8 | **Valid.** | New EBTRK parser, structural fit, config routing, and domain handling would create new code paths. | Pipeline regression is now a required hold item before treatment interpretation. |
| 9 | **Partly valid.** | Vmax effects are complex: older f-plane/Vmax decomp shows component changes, later gyre instrumentation reports mature steering bearing nearly flat across Vmax. | Hold test must run at storm-relevant Vmax and report N/W components before predictions. |
| 10 | **Valid.** | No EBTRK parser exists yet. The recipe only had a prose warning about fixed-width radii. | Recipe now specifies fixed-width quadrant parsing, positive-radii inclusion, missing sentinels, and source-code logging. Parser still pending. |
| 11 | **Partly valid.** | Hugo's `4444` source code applies to RMW/eye/POCI/ROCI, not directly to wind radii, but the record is still legacy and POCI/ROCI are unphysical. | Recipe keeps Hugo wind radii usable but confidence-flagged; POCI/ROCI must not feed pressure/domain calculations. |
| 12 | **Valid future-code risk.** | Current `choose_domain()` depends on global `R_ENV_M`; structural init is not implemented yet, but a naive implementation could change domains. | Recipe now requires fixed domain/grid per storm across control/treatment. |
| 13 | **Documentation risk.** | `CHEATSHEET.md` notes the m=1 phase sign bug, but not every downstream narrative clearly distinguishes post-fix from pre-fix. | Backlog: add a short gyre-figure provenance note when the mechanism section is rewritten. |
| 14 | **Valid.** | "Nonzero quadrant" was underspecified. | Recipe now says available positive finite quadrants only; log raw quadrants and missing handling. |
| 15 | **Valid.** | Log-radius fitting objective was underspecified. | Recipe now gives the nonlinear objective function. |
| 16 | **Valid.** | "Both outcomes are publishable" weakened falsifiability. | Recipe now commits: if treatment changes are all within ±20 km or signs fail, do not claim a β-drift landfall-error mechanism. |
| 17 | **Valid.** | `CHEATSHEET.md` had a six-storm banner but a stale "FULL ARC CLOSED" version line and old body text. | Header/version patched; deeper rewrite still needed before paper packaging. |
| 18 | **Valid nuance.** | Rmax=75 km is a resolution floor, not observed structure. | Principle and honesty ledger now say outer radii are observed; Rmax is a documented resolution floor unless observed Rmax exceeds it. |
| 19 | **Valid lower-priority doc issue.** | Acronyms such as DLM/ROCI are scattered and not always defined near first use. | Backlog: add glossary before publication. |
| 20 | **Valid.** | Fran prediction was over-specific despite fast translation caveat. | Predictions suspended. |
| 21 | **Mostly already handled.** | No top-level `run_hugo.py`/`run_katrina.py` live in `oracle_v8`; legacy scripts are under `oracle_v8/legacy`, have deprecation banners, and require `ORACLE_ALLOW_LEGACY=1`. | Fixed one stale legacy Ivan comment (`0.50` → `0.40`). |
| 22 | **Partly stale, still a doc cleanup.** | f-plane decomposition reports near-zero f-plane drift, so a self-translation leak is not currently evidenced. Older "energy audit next" language remains in the cheatsheet body. | Backlog: remove stale "leak suspect" phrasing during cheatsheet rewrite. |

## Required next tests before structural treatment

1. Gate-beta size/taper sweep at fixed storm-relevant Vmax, reporting north and west components.
2. ERA5 DLM sensitivity, at least Ivan, using wider annulus or storm-removed/lagged environmental steering.
3. Structural-init control regression: force `R_env=500 km`, `taper_start=200 km`, `B=1.5`, `Rmax=75 km`
   through the new path and verify it reproduces the existing control route.
4. Fixed-domain check: record `nx`, `Ly`, and `dx` for control and treatment; primary comparison must keep
   them identical per storm.
5. Prediction refreeze: one-sided signs/bounds written before any treatment run.
