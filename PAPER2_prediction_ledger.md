# Paper 2 — registered-prediction ledger (supplementary material candidate)

*Compiled 2026-07-07 from the campaign records (OVERROTATION_CANDIDATES.md,
ENVELOPE_INTENSIFICATION.md, DZ_HEATED_SENSITIVITY.md). Every prediction below was written and
committed BEFORE the run it concerns; outcomes quote the scored appends. This file replaces the
"fourteen predictions … roughly half failed" placeholder in PAPER2_draft.md and is a candidate
supplementary exhibit — it is the methodology section, in table form.*

## The ledger

| # | ID | registered claim (compressed) | conf. | outcome |
|---|---|---|---|---|
| 1 | P-A1 | ramp forms agree within ±3° / ±0.08 m s⁻¹ west | 70% | **CONFIRMED** (4° / 0.12 span) |
| 2 | P-A2 | envelopes track size, ≲5° from control | 60% | **FAILED** (21–24° NW rotation) |
| 3 | P-A3 | cutoff mechanism: ≥10° NW, west +0.2 | 20% | **CONFIRMED** — the surprise branch |
| 4 | P-C1 | passive null ≡ dry control | 95% | **CONFIRMED** (exact, 0.0°) |
| 5 | P-C2 | τ=30 min row ≈ control | 85% | **CONFIRMED** |
| 6 | P-C3 | baroclinic verdict (40/50/10 split) | — | **UNREADABLE** (the 10% branch; θ′→95 K) |
| 7 | P-E1 | heading NW-band across r_d 350–560 | 75% | **CONFIRMED** (4° span) |
| 8 | P-E2 | heading ±8° across Vmax 21–64 | 70% | **CONFIRMED** (4–5°) |
| 9 | P-E3 | f-plane null ≲ 0.05 m s⁻¹ | 90% | **CONFIRMED** (0.002) |
| 10 | P-S1 | along-track improves Katrina/Michael/Laura; Michael timing ≥1.5 h | 70% | **PARTIAL** (core 3/3 ✓; timing ✗) |
| 11 | P-S2 | H1 strong form (60–150 km west shifts) | 45% | **FAILED** |
| 12 | P-S3 | H2 buffer: < half strong form, one sign | 35% | **CONFIRMED** (ratio 0.34, sign 6/6) |
| G | P-S4 | guards: ΔVmax ≤ 10; timing ±3 h | guards | timing held; Vmax guard violated by Hugo/Ivan — **caught real confounds (functioned)** |
| 13 | P-EI1 | f-plane row A spins up ≥60; B ≥10 below | 65% | **FAILED** (regime error) |
| 14 | P-EI2 | ring import accounts ≥70% of A−B gap | 50% | **UNSCOREABLE** (no gap to attribute) |
| 15 | P-EI3 | C between A and B or above B | 55% | **CONFIRMED (trivial** — 2 m s⁻¹ spread) |
| 16 | P-N1 | D ≤ 50; E ≈ F; D−E ≥ 10 | 75% | **FAILED** in vortex metric (surface ordering confirmed) |
| 17 | P-N2 | F 3–10 m s⁻¹ below A | *(none stated)* | **FAILED** (0.1) — also a compliance slip: registered without a confidence |
| 18 | P-R3N1 | nz=64 historical reproduces 48.1 | 90% | **CONFIRMED** (bit-exact) |
| 19 | P-R3N2 | kill shot: normalized ≥65, ≥15 above | 60% | **FAILED** (+6.6 of 36.2) |
| 20 | P-R3N3 | nz=32 normalized 3–12 below historical | 65% | **FAILED by overshoot** (−14.7) |
| 21 | P-J21 | compact J2 peak ≥ 75 | 80% | **CONFIRMED** (82) |
| 22 | P-J22 | envelope ≥10 below on peak AND end | 55% | **FAILED as written** (delay, not suppression — wrong metric) |
| C | P-J23 | conditional, no numeric | — | not counted |
| 23 | P-M1 | Run-3 trajectories reproduce (onsets in bands) | 85% | **PARTIAL** (peaks ✓; onsets shifted by config seeds) |
| 24 | P-M2 | BL M-import rise leads onset; compact ≥2× | 50% | **FAILED** (export pre-onset; flip coincident) |
| 25 | P-M3 | discriminator; stated lean: reservoir location | 45% | **RESOLVED against the lean** — inflow strength (×2.5–3); the instrument decided |
| 26 | P-D1 | 32→64→96 monotone, shrinking increment | 55% | **CONFIRMED** (increment 14.9 → 0.1) |
| 27 | P-D2 | forcing-sampling control ≈ R2 within 2 | 70% | **CONFIRMED** (0.4) |
| 28 | P-D3 | max\|w\| ordering tracks max\|u\| | 75% | **CONFIRMED** (in substance) |
| 29 | P-C2v1 | held-core null ≡ control | 90% | **CONFIRMED** (exact) |
| 30 | P-C2v2 | bounded θ′/w, no cap-pinning | 75% | **PARTIAL** (bounded ✓; cap-pinned — maintained core is an energy source) |
| 31 | P-C2v3 | flat < 5° / < 0.15 ⇒ exonerate | 70% | **FAILED as written**; verdict obtained by the unregistered SIGN branch (poleward ⇒ exonerated as cause) |
| P | P-C2v4 | matched-intensity cleanup (dry cap-pinned row ≥354°, west ≤0.25) | 60% | **PENDING** (gate-beta rerun in progress) |

## Tally (31 scored items; guard set, conditional, and pending item counted separately)

**15 confirmed** (one trivially) · **10 failed** · **3 partial** · **2 unscoreable as designed**
· **1 discriminator resolved against the stated lean.**

## Calibration analysis

- **High-confidence bucket (≥80%, 7 items):** six confirmed, one partial — verified at
  ~93% against stated ~87%. Well calibrated, slightly under-confident.
- **Middle band (65–75%, 11 items):** five confirmed, two partial, four failed — realized
  ~59% against stated ~70%. **Overconfident**, and the misses are not random: they concentrate
  in exactly the two lessons the campaign recorded en route — regime matching (P-EI1, P-N1)
  and endpoint-versus-trajectory registration (P-J22, P-C2v3's missing sign branch).
- **Long shots (≤60%, 11 items):** behaved like long shots — and the two deliberate surprise
  branches both HIT: **the 20% cutoff branch (P-A3) and the 35% steering-buffer branch (P-S3)
  are the paper's two principal results.** The headline findings entered the campaign as its
  least-favored registered hypotheses, which is the strongest available evidence that the
  hypothesis space was honest rather than stacked.
- **Compliance slips, both minor, both recorded:** one prediction registered without a stated
  confidence (P-N2); one directional quantity registered without a sign-resolved branch
  (P-C2v3).

## Recommended manuscript wording (applied to PAPER2_draft.md this date)

§1: "Of thirty-one quantitative predictions registered across the campaign, fifteen were
confirmed, ten failed, and six were partial, unscoreable, or resolved against our stated lean…"

§6: exact counts + the calibration sentence + the long-shot observation (see draft).
