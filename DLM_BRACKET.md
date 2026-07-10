# DLM double-counting bracket — design + registered predictions

**Date:** 2026-07-07 · **Status:** knobs built into run_storm (env-var, bit-identical when
unset); predictions registered BEFORE any GPU run (do not edit after; append results). The last
construction item from the post-paper menu; feeds paper-1 §5.3's deferred caveat and paper-2
§6's bracketing hook. Discharges red-team hold item 2.

## The two entangled questions

1. **Is the DLM clean?** The 3–7° annulus, sampled near the storm, may contain part of the
   storm's own outer circulation — "environmental steering" that is partly the storm steering
   itself (the double-counting suspicion, raised for Ivan in the V8.6 era).
2. **Is the steering buffer position-feedback?** Paper 2's transmission ratio (~0.34) is
   attributed to the relaxation sampling the environment at the storm's actual position —
   a displaced storm samples a different environment and is pulled back. That attribution has
   not been tested directly.

## Design (8 runs; all knobs banner-printed — CHECK THE BANNER, per house tradition)

**Arm 1 — annulus sweep (question 1):** Ivan, Katrina, Laura × annulus {5–9°, 7–11°} = 6 runs,
under the production envelope profile. Controls = the existing `*_gauss_envelope` logs (3–7°).
```powershell
$env:ORACLE_DLM_INNER="5"; $env:ORACLE_DLM_OUTER="9";  python -m oracle_v8.run_storm <storm>
$env:ORACLE_DLM_INNER="7"; $env:ORACLE_DLM_OUTER="11"; python -m oracle_v8.run_storm <storm>
```

**Arm 2 — obs-anchored steering (question 2):** Katrina × {envelope, compact} with the DLM
sampled along the OBSERVED track (feedback loop severed) = 2 runs. **Clear the annulus vars
first** (`Remove-Item Env:ORACLE_DLM_INNER, Env:ORACLE_DLM_OUTER`).
```powershell
$env:ORACLE_STEER_ANCHOR="obs"; python -m oracle_v8.run_storm Katrina
$env:ORACLE_OUTER_ENVELOPE_M="none"; python -m oracle_v8.run_storm Katrina   # compact arm
```
The measurand in Arm 2 is the **A/B delta only** — absolute errors under obs-anchoring are a
different architecture and NOT comparable to production skill numbers.

## Registered predictions (Claude, 2026-07-07)

- **P-DLM1 (~65%):** the direct storms (Katrina, Laura) are annulus-robust — landfall-fix
  cross-track within 15 km and timing within 1 h of their 3–7° references, at both wider
  annuli. The production sampling is then clean-enough, and paper-1 §5.3's caveat closes with
  a number.
- **P-DLM2 (Ivan, split lean):** Ivan — slowest mover, largest circulation, the original
  suspect — shifts more than the direct storms. Branches: <15 km (clean, ~35%), 15–40 km
  (mild contamination, quantified, ~50%), >40 km or >1.5 h timing (the double-counting
  suspicion was right and matters, ~15%). Direction lean if it shifts: wider annulus → less
  storm-signal → *less* poleward over-run (timing less early).
- **P-DLM3 (~60%):** with the feedback severed, the Katrina envelope-vs-compact cross-track
  delta rises from the buffered −46 km toward the strong form −108: |Δ| ≥ 90 km ⇒ **the buffer
  is position-feedback, confirmed** — paper 2 §5.3's mechanism sentence gets its direct test.
  |Δ| < 70 ⇒ the buffer lives elsewhere (e.g., in the relaxation timescale itself) — rewrite
  the sentence honestly. (70–90: gray zone, decompose per-segment.)
- **P-DLM4 (guards):** Vmax histories within ±10 m s⁻¹ of the matching-profile references;
  every run's banner shows the intended knobs (the stale-env hazard is real and has bitten
  once); obs-anchored runs additionally sanity-checked against the F-V5 lesson — the observed
  track here is HURDAT2, used deliberately as a forcing location, never as a verification
  reference.

**Decision rules:** P-DLM1+2 pass clean ⇒ annulus caveat closes; contamination branches ⇒
paper-1 §5.3 and paper-2 §6 get the measured number instead of the caveat. P-DLM3's first
branch ⇒ transmission ratio decomposed (feedback confirmed as the mechanism); second branch ⇒
the ratio stands but its attribution changes. All branches publishable; none blocks either
paper.

## Results

*(append after the GPU runs; predictions frozen)*
