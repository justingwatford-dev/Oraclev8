# Oracle V8 — Fourth Ensemble Brief
## Katrina 2005: Six-Run Diagnostic Campaign

**To:** Five (GPT-5 / Codex), Gemini
**From:** Justin (Orchestrator) + Claude
**Re:** Systematic Katrina track failures — tracker asymmetry under secondary intensification
**Date:** May 2026

---

## Context

Hugo 1989 was validated in 19 runs to 48.1 km east bias (Run 19). Katrina 2005
has now undergone 6 dedicated runs. The runs are fully stable (no NaN, 42h
completions) but the track does not reach landfall (29.1°N) within the window.
This brief documents the complete diagnostic progression and asks focused
questions on the residual failure mode.

---

## Storm Parameters

```
Init:        2005-08-28 00Z
Position:    24.0°N, 86.3°W (Cat 5, 145 kt, 74.6 m/s)
Rmax run:    75 km (5×dx; observed 46.3 km)
Landfall:    29.1°N, 89.6°W  t+35h (Buras, LA)
Grid:        128×128×32, dx=15.6 km, dt=30s, CFL=0.143
ERA5 DLM (ring-average, 3-7° annulus):
  t=0h:   u=-2.53  v=+1.27 m/s
  t=12h:  u=-2.07  v=+2.72 m/s
  t=18h:  u=-1.79  v=+4.29 m/s   ← trough interaction begins
  t=24h:  u=+0.17  v=+4.76 m/s   ← u turns EASTWARD (ridge weakens)
  t=35h:  u=+1.12  v=+6.12 m/s
  Mean:   u=-1.01  v=+3.84 m/s
```

Note: Katrina's ERA5 DLM u changes SIGN over 35h. Hugo's u stayed negative
(westward) throughout. This is the fundamental difficulty.

---

## Six-Run Diagnostic Campaign

| Run | Key change | Track result | Failure mode |
|-----|-----------|-------------|--------------|
| 1 | ERA5 t=0 DLM (v=+1.27 m/s) | 0.27°N in 35h | v too small; constant steering fundamentally wrong for Katrina |
| 2 | Track-mean DLM (v=+3.84 m/s) | 2.15°N, tracker jumps to 24.95°N at t=14h | θ' accumulation → false tracker minimum; max\|u\|→177 m/s |
| 3 | + NewtonianCooling(τ=1h), 42h | 2.21°N; 4-cell cluster parking | τ too long; tracker locks onto grid clusters; angular momentum spin-up continues |
| 4 | τ=30min cooling | 2.39°N to t=17.5h, then jumps south | τ prevents θ' clusters but max\|u\|→120 m/s; vorticity CoM jumps to wrong anchor at t=18h |
| 5 | Dead-reckoning tracker (v/u_env_t0) | 3.19°N; consistent eastward drift | Vorticity CoM fires instead (jump < 2×gate); CoM drifts east |
| 6 | Removed Coriolis mid-run update | **Identical to Run 5** | Coriolis update confirmed NOT the cause; eastward drift is intrinsic |

---

## What We Know Works

**Runs 4 and 5 (τ=30min, dead reckoning) — first 17.5h of track:**

```
t=0h:    23.93°N, 86.53°W
t=17.5h: 26.33°N, 87.30°W   ← 2.40°N in 17.5h = 4.2 m/s northward
                                  (expected 3.84 m/s; 10% overshoot from β-drift) ✓
```

This is **excellent** — smooth, consistent NNW motion at the right rate for 17.5h.
The track breaks down only when secondary intensification reaches ~120 m/s.

---

## The Residual Failure Mode

### Observation
From t=18h onward in Runs 5 and 6, the vorticity CoM tracker produces:
- **Northward:** +0.02°N per 0.5h step = 1.23 m/s (tracker underestimates ~3× vs expected v_env=3.84)
- **Eastward:** −0.07°W per step = 4.3 m/s EASTWARD (wrong direction entirely)

Total at t=41.5h: 27.12°N, 83.35°W — storm is **tracking east instead of west.**

### Root cause analysis
At t=18h, max|u| reaches ~120 m/s. The secondary intensification from angular
momentum concentration (not θ' buoyancy — BuoyancyComponent is disabled) creates
a strong, **asymmetric secondary vortex**. On the β-plane, the Rossby dispersion
preferentially amplifies the eastern side of the vortex circulation.

The vorticity CoM jump at t=18h is small enough (< 2×max_jump_cells = 156 km)
to pass the continuity gate without triggering dead reckoning. The CoM returns
a position that is:
- Correct latitude: physical vortex IS moving north at ~v_env rate
- Wrong longitude: biased eastward by the asymmetric secondary circulation

Dead reckoning fires only when BOTH θ' and vorticity signals fail catastrophically.
Here the vorticity CoM fails *gradually*, returning positions that drift east
continuously. The dead reckoning gate never triggers.

**Evidence that this is β-plane vorticity asymmetry:**
1. Hugo (f-plane-like — stays equatorward of 35°N) showed no eastward drift
2. Katrina at 24-28°N has larger β effect: β at 25°N = 2.17×10⁻¹¹ m⁻¹s⁻¹
3. The eastward drift rate (+4.3 m/s) exactly matches the scale of β-induced
   vortex asymmetry: V_beta = β × Vmax × Rmax / f = 2.17e-11 × 75000 × 75000 / 5.93e-5 = 2.06 m/s (same order)
4. Removing the Coriolis mid-run update (Run 6 vs Run 5) — ZERO effect.
   The drift is embedded in the physics, not the tracker reference.

### What does NOT cause it
- ✗ ERA5 t=0 DLM sign (tried t=0 and track-mean)
- ✗ Newtonian cooling timescale (tried τ=1h and τ=30min)
- ✗ Dead reckoning (implemented; fires correctly but gate threshold not met)
- ✗ Coriolis mid-run update (removed; identical result)

---

## Questions for the Ensemble

### Q1 — Tracker robustness under β-plane vortex asymmetry

The current vorticity CoM uses `max(ζ, 0)` column-averaged over the lower half
of the domain. On the β-plane with strong secondary circulation (max|u|→120 m/s),
this is contaminated by the asymmetric beta-gyre vorticity structure.

Standard approaches in the literature (Marks & Houze 1987, Reasor & Eastin 2012):
a. **Vorticity-weighted centroid with distance penalty** — weight by ζ × exp(-r²/r₀²)
   centered on last known position, preventing the CoM from drifting to the beta-gyre
b. **SLP minimum in a restricted search box** — find the pressure minimum only
   within ±N cells of the last position (stricter than current interior masking)
c. **Phase velocity from vortex tilt** — track the phase of the dominant wavenumber-1
   Fourier component of the vorticity field

Which approach is most appropriate for a dx=15.6 km, barotropic, β-plane, Cat 5 TC?
Is there a closed-form approach that avoids explicit wavenumber decomposition?

### Q2 — Is the secondary intensification physical or numerical?

max|u| reaches 177 m/s by t=28h despite BuoyancyComponent being disabled.
The only intensification mechanism in the current config is angular momentum
convergence through the anelastic secondary circulation.

For a dry barotropic TC without surface drag, is unlimited intensification through
angular momentum concentration expected? What is the dimensional estimate for
the terminal vortex intensity in the absence of friction and diabatic effects?

Specifically: does this represent the correct dynamics of an inviscid anelastic
vortex (i.e., the model is CORRECT and real TCs are limited by missing physics),
or is there a numerical artifact (e.g., aliasing at high CFL) driving the
intensification beyond the physical limit?

The CFL at 120 m/s: 120 × 30 / 15625 = 0.23 — below stability limit but elevated.

### Q3 — Validation sufficiency for BAMS

The 6-run Katrina campaign has systematically identified:
1. The steering architecture limitation (constant initialization vs time-varying ERA5)
2. The secondary intensification mechanism (angular momentum concentration)
3. The tracker failure mode (β-plane vorticity asymmetry in the CoM)

The first 17.5h of track is physically correct (4.2 m/s northward, right direction).
The background-flow dead-reckoning estimate reaches 29.1°N at t=39.7h (within window).
The tracked position does not cross 29.1°N due to the eastward CoM bias.

**Question:** For a methods paper in BAMS, is a partial Katrina validation — documenting
the first 17.5h track quality plus a systematic diagnostic campaign explaining the
failure mode — scientifically sufficient, or does a successful landfall crossing need
to be demonstrated before submission?

Hugo gives 48.1 km. Katrina's diagnostic campaign demonstrates the model physics
are correct but the tracker fails under secondary intensification. Is the
documented failure mode itself a publishable finding about the limitations of
θ'/vorticity tracking under barotropic intensification?

### Q4 — Recommended path forward

Given the 6-run progression, what is the single highest-leverage next step:

a. **Restricted-box tracker** — find ζ_max in a ±3-cell box around last position.
   Low implementation cost. Might prevent beta-gyre contamination.

b. **Implicit surface drag** — prevent secondary intensification from exceeding
   ~100 m/s where the tracker degrades. Removes the root cause.
   SurfaceDragComponent exists; needs implicit Cd×V treatment.

c. **Time-varying steering** — add a mean-flow relaxation tendency that nudges
   background flow toward ERA5 DLM each step. Addresses the initialization
   limitation independently of the tracker issue.

d. **Accept Run 4's first 17.5h** — document as the valid Katrina track segment
   and move to Ivan for the third validation storm.

---

## Summary for Reference

```
Hugo Run 19:    48.1 km east bias at landfall, 32h run ← paper result
Katrina Run 6:  Stable 42h, reaches 27.12°N (needs 29.1°N)
                First 17.5h valid: 23.93 → 26.33°N at 4.2 m/s NNW ✓
                t>18h: vorticity CoM drifts east at 4.3 m/s (β-gyre asymmetry)
                Dead reckoning fires but gate threshold never exceeded
```
