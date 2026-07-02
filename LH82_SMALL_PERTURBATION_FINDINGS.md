# LH82 small-perturbation study — findings

**Date:** 2026-07-01 · **Status:** Phase 1–3 complete; conclusion banked.

## Question

The LH82 anelastic system rests on the small-perturbation assumption
**|θ′| ≪ θ̄** (θ̄ ≈ 300 K): it uses the base-state density ρ̄(z) in continuity
and a linearized buoyancy `b = g·θ′/θ̄`, dropping higher-order terms. LH82 flags
the eyewall — where θ′ can reach 10–20 K (≈ 3–7 %) — as where this may break.
The empirical question: **does our eyewall reach that regime, and how large are
the neglected terms when it does?**

## Setup

- Grid 128×128×32, dx = 15.625 km, dz = 625 m, dt = 30 s (15 s in one sensitivity run).
- Config: **LH82 anelastic + `upwind5h` advection + BuoyancyComponent + prescribed
  annular eyewall heating (`DiabaticHeatingComponent`, ring at r ≈ Rmax = 75 km,
  z_peak = 5 km) + surface drag + hyperdiffusion + Newtonian cooling (τ = 1800 s)**.
  Helmholtz divergence damper ε swept (0, 0.1, 0.5).
- Balanced Holland vortex (Vmax 64 m/s, Rmax 75 km); the 42 K balanced init core
  relaxes within ~5τ, so steady-state θ′ is heating-controlled (θ′_eq ≈ Q·τ_cool).
- Small parameter reported as **max(θ′/θ̄)** over the domain (worst case).
- Harnesses: `oracle_v8/measure_lh82.py` (Phase 1), `oracle_v8/measure_lh82_phase2.py`
  (Phase 2), `oracle_v8/measure_lh82_phase3.py` (Phase 3, grid-parametrized).
  Run on GPU with `$env:LH82_STEPS = 1000`.

## Phase 1 — regime map (ε = 0, dt = 30, 1000 steps)

| Q (K/s) | max\|u\| | max\|w\| | max\|θ′\| | **max(θ′/θ̄)** | outcome |
|---|---|---|---|---|---|
| 2.5e-3 | 39.9 | 0.35 | 3.1 K | 0.99 % | survived |
| 5.0e-3 | 48.0 | 1.08 | 6.1 K | 1.91 % | survived |
| 1.0e-2 | 74.1 | 5.01 | 12.6 K | **4.00 %** | survived |
| 1.5e-2 | — | — | — | — | **BLEW UP @705 (5.9 h)** |
| 2.0e-2 | — | — | — | — | **BLEW UP @466** |

- θ′/θ̄ climbs ~linearly with forcing into LH82's flagged band (4 % at Q=1e-2).
- **At Q=1e-2 the eyewall is physically credible**: max|w| = 5 m/s and the vortex
  *intensifies* 64 → 74 m/s (secondary-circulation angular-momentum convergence,
  WISHE-like) — not a contrivance.
- The model appeared to destabilize just past ~5 % → the key question for Phase 2:
  physical (LH82 breakdown) or numerical?

## Phase 2A — the blow-up is NUMERICAL, not LH82 breaking down

Re-ran the blow-up cases against numerical knobs:

| config | outcome | max\|u\| | max\|w\| | max(θ′/θ̄) |
|---|---|---|---|---|
| Q=1.5e-2  ε=0.0  dt=30 | BLEW UP @705 | — | — | — |
| Q=1.5e-2  ε=0.1  dt=30 | survived | 67.1 | 1.38 | 6.41 % |
| Q=1.5e-2  ε=0.5  dt=30 | survived | 65.8 | 0.32 | 7.46 % |
| Q=1.5e-2  ε=0.0  dt=15 | survived | 103.6 | 10.67 | 7.17 % |
| Q=2.0e-2  ε=0.5  dt=30 | survived | 71.4 | 0.53 | **9.74 %** |

**Both damping (ε > 0) and halving dt independently cure the blow-up.** A genuine
physical breakdown would survive neither. So the "~5 % boundary" was the ε = 0
setup outrunning the upwind dissipation at strong forcing — a numerical artifact,
correctly identified rather than mistaken for physics. **LH82 runs stably to
θ′/θ̄ ≈ 10 %**, well past its own flagged 3–7 % band. We did not find its edge.

## Phase 2B — neglected-term magnitudes (stable Q=1e-2, θ′/θ̄ = 4 %)

```
buoyancy relative error   θ′/θ̄                    = 4.00 %
continuity error          max|∇·(ρ′u)| / max|∂(ρ̄w)/∂z|  = 3.25 %
```
(ρ′ ≈ −ρ̄·θ′/θ̄, dry constant-pressure estimate.) Both neglected-term ratios are
**~3–4 %**, consistent with the leading-order expectation that they scale as θ′/θ̄.

## Phase 3 — tightening the two soft spots before banking

Phase 2 left the ~10 % point resting on heavy damping (ε = 0.5 → max|w| = 0.53)
and the 4 % ratios unchecked against resolution (the heating ring, width_r = 30 km,
is only ~2 cells at dx = 15.6 km). Phase 3 addressed both.

### 3A — the ~10 % state with a live(r) circulation (Q = 2e-2, dt = 15)

| config | outcome | max\|u\| | max\|w\| | max(θ′/θ̄) | continuity err |
|---|---|---|---|---|---|
| ε = 0.0 | **BLEW UP @1410 (5.9 h)** from a ~6 % state | — | — | — | — |
| ε = 0.1 | survived | 74.2 | 1.55 | **9.26 %** | **7.90 %** |

- At Q = 2e-2, halving dt alone no longer rescues ε = 0 (it did at Q = 1.5e-2):
  the undamped-upwind ceiling sits between Q = 1.5e-2 and 2e-2 even at dt = 15.
  But the run went unstable while θ′/θ̄ was only ~6 % (heartbeats: 5.66 % @500,
  6.08 % @1000 with max\|u\| run up to 108.8) — *below* the 9.26 % at which the
  ε = 0.1 run sits happily in steady state. The instability is not
  θ′-amplitude-driven: numerical, again.
- Mild damping now carries the ~10 % claim with a moderately live circulation:
  max\|w\| = 1.55 (3× the ε = 0.5 value), and the vortex intensifies 64 → 74.
  Heartbeats (8.45 → 8.63 → 8.77 → 9.26 %) show a clean approach to equilibrium.
- **Second point on the scaling line**: continuity error 7.90 % at θ′/θ̄ = 9.26 %
  → ratio 0.85, vs 0.81 at the 4 % point. The neglected terms grow **linearly**
  (≈ 0.8–0.9 × θ′/θ̄) over a 2.3× amplitude range — no superlinear growth, no
  sign of an approaching breakdown.

### 3B — grid convergence of the 4 % ratios (Q = 1e-2, ε = 0)

| grid | dt | max\|u\| | max\|w\| | max(θ′/θ̄) | continuity err |
|---|---|---|---|---|---|
| 128²×32 (base = Phase 2B) | 30 | 74.1 | 5.01 | 4.00 % | 3.25 % |
| 128²×32 (dt control) | 15 | 84.3 | 10.68 | 3.87 % | 3.16 % |
| 256²×32 (dx/2) | 15 | 68.0 | 3.49 | 3.40 % | 3.12 % |
| 128²×64 (dz/2) | 15 | 48.1 | 3.80 | 2.94 % | 1.49 % |

- The baseline rerun reproduces Phase 2B exactly (4.00 / 3.25) — deterministic.
- The dt control leaves both ratios essentially unchanged while doubling
  max\|w\| (5.0 → 10.7): the **ratios are robust to the numerics that the
  circulation is sensitive to**.
- dx/2 (ring 2 → 4 cells): 3.40 / 3.12 % — converged (≤ 15 % relative).
- dz/2 moves both ratios **down** (continuity halves to 1.49 %) — but alongside
  a genuinely different vortex equilibrium: at NZ = 64 the vortex *decays*
  64 → 48 m/s under the same forcing that intensifies it at NZ = 32 (suspect the
  surface-drag layer H_BL = 1000 m spanning 1.6 vs 3.2 cells; not chased here).
  Direction is conservative either way: **refinement reveals smaller neglected
  terms, not larger** — the coarse-grid headline numbers are upper bounds.

## Conclusion

**In this dry idealized eyewall, LH82's small-perturbation assumption holds up
well.** The neglected terms grow **linearly** (≈ 0.8–0.9 × θ′/θ̄) from the
realistic 4 % eyewall to the strongest steady state we can hold (9.26 % at
ε = 0.1), with no superlinear growth or qualitative misbehavior, and the
headline ratios are grid-robust (refinement makes them smaller, not larger).
There is no visible breakdown in the reachable regime — a *negative result for
"LH82 fails," which validates LH82 here*. Accordingly, **the motivation to
derive/implement PseudoIncompressible is weak on current evidence** — this
diagnostic was designed to decide exactly that, and it says PI is not worth
the derivation yet.

### Why we stopped at ~10 % (scope decision, 2026-07-01)

- LH82 is an asymptotic expansion: errors scale linearly with θ′/θ̄ and there is
  no sharp edge to find. Pushing θ′ higher confirms an extrapolation (now
  verified to 9.3 %), and θ′/θ̄ > 10 % is outside observed eyewall amplitudes —
  validity there is decision-irrelevant for how this model is used.
- A moist extension at dx = 15.6 km would swap the explicit Q knob for closure
  knobs (θ′ no more "genuine") and add virtual-temperature / water-loading terms
  of the *same order* as the LH82 correction, destroying the clean isolation.
  Moisture is a model-roadmap decision, not an LH82 follow-up.
- The one failure mode this diagnostic cannot see — accumulated ∇·(ρ′u) mass
  error feeding back on the solution — is only measurable by implementing PI and
  diffing solutions, which is precisely what current evidence says isn't
  warranted.

## Caveats

1. **Dry & idealized** — prescribed heating, no moisture/latent heating. Real moist
   eyewalls could push θ′ higher or add effects this setup can't show; the ~10 %
   ceiling is a dry-model ceiling.
2. **Leading-order diagnostic** — captures the dominant neglected term, not the
   full LH82-vs-PI difference (exact pseudo-density form). A rigorous split still
   needs PI; the point is that it isn't urgent.
3. **Numerics strongly shape the circulation** — max|w| ranged 0.32 (ε=0.5) to
   10.7 (ε=0, dt=15) at the *same* forcing (20×), and at NZ=64 the vortex
   *decays* (64 → 48) where NZ=32 intensifies (→ 74–84) — even the sign of
   intensification is resolution-sensitive (suspect the H_BL=1000 m drag-layer
   discretization; unexplained). Phase 3B shows the *diagnostic ratios* are
   robust to all of this, but updraft strength / intensification must be flagged
   if either ever becomes a paper quantity.
4. The ~10 % ceiling is numerical, not a proven physical limit: at dt=15 the
   undamped (ε=0) ceiling sits between Q=1.5e-2 and 2e-2, while ε=0.1 holds
   9.3 % in steady state. Finer grid / flux-form conservative advection would
   push it, but there is no scientific need (see scope decision).

## Methodological note

This nearly produced a false "LH82 breaks at 5 %" result; the physical-vs-numerical
disambiguation (vary ε *and* dt) caught it as numerics. Same discipline as the
red-team pass: measure, and rule out the numerics before claiming physics.
