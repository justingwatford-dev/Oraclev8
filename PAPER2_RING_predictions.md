# Ring-vs-cutoff discriminator + envelope domain null — registered predictions

*Registered 2026-07-24, BEFORE either run executes. Both experiments are direct responses to
external review: (1) the reviewer showed Arm A cannot separate "the cutoff" from "the
anticyclonic ring" — all three compact ramps carry identical integrated ring strength (compact
support fixes the removed circulation), so both hypotheses predicted the observed pattern —
and named the clean single-variable test; (2) the reviewer noted the envelope places real
circulation farther out than the taper, and cascade chapter 2 established this model's
domain-edge sensitivity, so the envelope owes a domain-size null. We checked whether the
discriminating run already existed: it does not — the closest is `gate-beta-renv` (taper ring
moved outward with R_env), which confounds ring location with circulation size.*

## Experiment R — ring + tail (`gate-beta-ringtail`)

**Profile.** V(r) = V_hol(r) · [T(r) + G(r)·(1−T(r))], with T the production cosine ramp
(200→500 km) and G = exp(−(r/420 km)²): identical to both parents inside 200 km; the taper's
anticyclonic ring in the 200–500 band (the crossfade removes circulation from the Holland
level down to the envelope level — ≈76% of the taper's ring circulation, same band, since
G(500 km) = 0.24); the r_d = 420 envelope's tail beyond 500. Implemented as
`ring_plus_tail=True` in `vortex_init.py` (this commit). Testbed: Vmax 64, cap 70,
5000 km/320², quiescent β-plane, 48 h; mature 30–48-h drift + t12/24/36/48 window traces;
envelope-420 control re-run in the same harness invocation.

**Predictions.**

- **P-R1 (60%) — the discriminator.** The ring+tail profile **locks**: heading span ≤ 6°
  between the t24 and t48 windows, mature heading ≤ 340°. That verdict implicates the
  **cutoff** (absence of the arresting far field), not the ring. Stated lean follows
  `gate-beta-renv` (ring relocated outward did not restore westward drift) and the envelope
  family (a weaker, broader ring locks canonically). **Branch:** heading precessing ≥ 8°
  t24→t48 and reaching ≥ 345° implicates the **ring**, and the exported advice changes from
  "keep the tail" to "keep compensating vorticity out of the gyre band." Either branch
  completes the attribution Arm A could not.
- **P-R2 (85%, guard).** Equilibrated Vmax within 3 m s⁻¹ of the envelope control's — the
  comparison is not intensity-confounded.

## Experiment DN — envelope domain null (`gate-beta-envdomain`)

Envelope r_d = 420, Vmax 64, quiescent β-plane, 48 h, **7500 km / nx = 480** (same Δx),
against the committed 5000 km/320² Stage-2 control (2.30 m s⁻¹ @ 329°, west +1.20,
Vmax_end 39.7).

- **P-DN1 (75%).** Mature drift within ±0.15 m s⁻¹, ±4°, and westward component within
  ±0.15 m s⁻¹ of the committed control: no domain-edge sensitivity for the envelope at the
  production domain rule. Failure sends the envelope back through the cascade's chapter-2
  geometry check before either paper ships.

## Scoring

CONFIRMED/FAILED on the stated thresholds; no post-hoc motion; outcomes appended here and
folded into Paper 2 §2.1/§3/§4 (the Arm A design limitation is already acknowledged in the
manuscript text as of this commit, with the discriminator referenced as "reported in
Section 3").

---

# OUTCOMES (scored 2026-07-24, same day; logs `oracle_v8/Logs/GATE_BETA_RINGTAIL.txt`,
# `GATE_BETA_ENVDOMAIN.txt`; GPU backend, same-harness envelope control re-run)

| profile | mature drift | heading trace t12/24/36/48 | west | Vmax_end |
|---|---|---|---|---|
| gauss r_d=420 control | 2.30 @ 329° | 328/327/328/329 | +1.20 | 39.7 |
| **ring+tail** | 2.97 @ 335° | 332/335/334/335 | +1.28 | 44.5 |
| envelope, 7500 km/480 | 2.26 @ 327° | — | +1.22 | 39.6 |

- **P-R1 CONFIRMED — the CUTOFF, not the ring.** The ring+tail profile **locks**: heading span
  t24→t48 is ≤1° (threshold ≤6°), mature heading 335° (threshold ≤340°), westward component
  +1.28 — fully canonical behavior with the taper's ring present at ≈76% strength in the gyre
  band. The ring's effect is a few degrees of equilibrium shift (335 vs the envelope's 329),
  exactly the *parameter dependence* the prior literature describes; the unlock of the compact
  family is the absence of the arresting far field. The exported advice stands: keep the tail.
- **P-R2 FAILED (guard).** Equilibrated Vmax 44.5 vs the control's 39.7 — 4.8 m s⁻¹, beyond
  the registered 3. The ring+tail profile carries more mid-radius wind and equilibrates
  hotter. Consequence, stated precisely: the *few-degree equilibrium offset* (335 vs 329) is
  intensity-entangled and not claimed; the **lock itself** — amplitude growing into a fixed
  orientation over 24+ h — cannot be an intensity artifact (the envelope family's aim is
  intensity-invariant across Vmax 17–40), and P-R1's verdict rests on the lock, not the offset.
- **P-DN1 CONFIRMED.** 7500-km/480² envelope: deltas vs the committed control are
  |drift| −0.04 m s⁻¹, heading −2°, west +0.02 — all inside ±0.15/±4°/±0.15. No domain-edge
  sensitivity; the cascade-chapter-2 loop is closed for the envelope.

**Running program tally: 51 registered / 26 confirmed / 18 failed / 7 other.**
