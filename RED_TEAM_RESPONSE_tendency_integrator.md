# Response to red-team review — `tendency.py` and `integrator.py`

Thank you — this is a sharp review, and your central instinct is correct: the
dynamical core has a stability problem that is currently managed with added
dissipation rather than fixed at the source. I ran a set of controlled ablations to
locate the mechanism precisely, and the data both confirms your root diagnosis
(point 1) and overturns one of your prescriptions (point 2). Details and numbers
below.

Framing point up front, same as the `poisson.py` round: **the production
configuration is barotropic** (`BuoyancyComponent` omitted, θ′≡0 at init). Verified
active components in `build_production_config`:

- FAST (PRE_PROJECTION): **none**
- PROJECTION: `anelastic_projection`
- SLOW: `advection, coriolis, surface_drag, hyper_diffusion, newtonian_cooling,
  helmholtz_divergence_damping`

(`pressure_gradient` and `buoyancy` are both `None`; the 2nd-order
`HorizontalDiffusionComponent` is **not** active — only hyperdiffusion;
`IntensityCapComponent` is applied in the run loop, not the tendency config.)

---

## Scorecard

| # | Claim | Verdict |
|---|---|---|
| 1 | Centered advection is unstable, needs upwind | **Confirmed by measurement — the genuine root cause** |
| 2 | Delete the Helmholtz damper; it suppresses the secondary circulation and the projection is enough | **Backwards** — it is vorticity-neutral *and* the single load-bearing stabilizer |
| 3 | Geostrophic-background subtraction is an "illusion" | Declined — standard technique, identical to your own proposed fix |
| 4a | Intensity cap masks a real instability | Conceded — but it is downstream of point 1 |
| 4b | Newtonian cooling pins θ′; −88 K implies broken thermodynamics | Reframed — θ′ is a passive tracer in production; legitimate for the buoyancy study |
| 5 | Lagged pressure drops RK3 to 1st order | Declined — targets a dormant path |

---

## The mechanism, measured

I isolated a Cat-4 Holland vortex (Vmax 64 m/s) on the production grid (128², dz from
NZ=32), barotropic, no intensity cap, and varied one thing at a time.

**(i) The instability is nonlinear aliasing of centered advection — a spatial-scheme
problem, not a time-scheme one.** With all dissipation off:

| Integration of identical centered advection (no dissipation) | Result |
|---|---|
| Forward-Euler / Strang (production time scheme) | blew up at step 150 |
| SSP-RK3, same advection operator | blew up at step 152 |
| Forward-Euler + hyperdiffusion (production coeff.) | blew up at step 163 |

FE and RK3 fail identically, with the same grid-scale energy growth. RK3's
conditional stability covers *linear* advection; it does nothing for the nonlinear
energy cascade into the 2Δx mode. So this is your point 1, confirmed: centered
advection without adequate dissipation is the root cause.

**(ii) The dissipation that actually controls it is the Helmholtz divergence damper,
not hyperdiffusion.** Ablating the full production stack (300 steps each):

| Variant | Result |
|---|---|
| full production stack | survived (max\|u\|=54.6) |
| − hyperdiffusion | **survived (max\|u\|=54.6, identical)** |
| − helmholtz divergence | **BLEW UP at step 151** |
| − surface drag | survived (max\|u\|=58.8) |
| − newtonian cooling | survived (max\|u\|=54.6) |
| − hyperdiffusion − helmholtz | BLEW UP at step 144 |

Removing hyperdiffusion changes the solution by nothing on this timescale; removing
the Helmholtz damper blows the run up. The reason is physical and is exactly the
observation in your point 2: the nonlinear aliasing of `u·∇u` pumps energy into
**divergent** grid-scale modes, and the Helmholtz damper removes that divergent part
aggressively (ε=0.5 ≈ 75%/step) while leaving the rotational vortex untouched
(measured Δζ = 0.07%). The projection enforces `∇·(ρ̄u)=0` once per stage, but the
aliasing regenerates divergence faster than a single projection removes it; the
damper supplies the missing per-step divergent dissipation.

---

## Point-by-point

### 1 — Centered advection. Confirmed; this is the root cause.

`_cx/_cy/_cz` are pure 2nd-order centered differences with no numerical dissipation,
and the ablations above show this is what drives the instability the rest of the
stack is compensating for. The right fix is upstream, in the advection operator —
either an odd-order upwind scheme (your suggestion; dissipation built in) or an
energy-conserving/skew-symmetric centered form with adequate explicit dissipation.
The one place I'd still soften "non-negotiable 5th-order upwind": for a methods paper
the explicit-dissipation route can be *more* transparent to document. But that is a
style choice; your core point stands and the data backs it.

### 2 — Helmholtz divergence damping. The prescription is backwards.

Two measured facts:

- It does **not** damp the secondary circulation / "the engine." It removes only the
  irrotational part of the horizontal wind (`∇ψ`, `∇²ψ=D`); since `∇×∇ψ=0` it adds
  ~zero vorticity. On the actual vortex, one ε=0.5 half-step changes `max|ζ|` by a
  relative **7.4e-4 (0.07%)**, with a 0.05 m/s wind-speed change against a 32 m/s
  vortex. In a barotropic run there is also no buoyancy-driven overturning for it to
  suppress.
- It is the **single load-bearing stabilizer**: removing it blows the run up at step
  151, while removing hyperdiffusion, drag, or Newtonian cooling individually does
  not. "Delete this; the projection is the only divergence constraint you need" is
  exactly inverted — the projection alone is demonstrably *not* enough.

Your sub-observation is correct and important — the advection *is* generating
spurious divergence, and this damper exists to remove it. That's precisely why the
real fix is point 1: a less aliasing-prone advection operator would generate less
spurious divergence and could let us reduce ε or remove the damper entirely. Until
then, it stays. (The "ε=1.0 worse than ε=0.5" you'd have seen historically was an
FD-vs-spectral inconsistency, since fixed to a fully spectral projection — not
"shutting down the updraft.")

### 3 — Geostrophic background. Declined.

Applying Coriolis to the perturbation `(u−u_env)` is the standard idealized-TC
treatment and is mathematically identical to your "Option 1": an implied constant PGF
that exactly balances Coriolis on the uniform background. A uniform background has no
spatial shear and manufactures none. Without the subtraction the whole airmass
inertially oscillates (T≈26 h) — the bug this fixed. The only real subtlety is
time-varying steering, where the reference and the in-state background must advance in
lockstep; `CoriolisComponent.set_env` documents exactly that hazard. Solved,
conventional, not an illusion.

### 4a — Intensity cap. Conceded — but downstream of point 1.

Agreed that an arbitrary speed cap is not publishable as an intensity control, and
that a 372 m/s runaway signals a real failure. Note from the ablations that the cap
is **not** required for short-run stability of a 64 m/s vortex (these 300-step runs
survive without it); the runaway it addresses appears at higher amplitude / longer
integration and is the same aliasing energy source as point 1. Fix the advection and
the cap should become unnecessary. Until then: track (the paper's subject) is far
less sensitive to peak intensity than intensity itself, and we'll say so explicitly
rather than lean on the cap.

### 4b — Newtonian cooling. Reframed.

In production, buoyancy is off, so θ′ never feeds back into `w`; it is generated by
the advection term `−w·∂θ̄/∂z` and used as a tracer for storm-center tracking. The
ablation confirms it is not load-bearing for stability (removing it: survives,
max\|u\| unchanged). So in production it is a *tracer bound*, and should be labeled as
such, not as thermodynamics. For the buoyancy-enabled study your deeper point is
legitimate: an updraft that keeps accelerating while cooling tens of K indicates a
buoyancy/stratification balance to audit, not a θ′ to relax.

### 5 — Lagged pressure. Declined (dormant path).

The lagged `projection_potential` is consumed only by `PressureGradientComponent`,
which the `OperatorConfig.__post_init__` guard forbids alongside `AnelasticProjection`
(double-counting). Production has no fast components; the buoyancy study has
`[buoyancy]`, which does not read φ. So the lagged-pressure mechanism never drives a
projection-based run — the projection (the real pressure correction) is solved fresh
every substage (`integrator.py:494`). As a separate observation: because FE and RK3
give identical stability here, the RK3 machinery isn't buying anything for barotropic
runs — a simplification opportunity, not a bug.

---

## Summary

Your root diagnosis (point 1) is correct and now demonstrated: centered advection's
nonlinear aliasing is the instability, and neither the time scheme nor hyperdiffusion
is what holds the model together — the **Helmholtz divergence damper** is the
keystone, because the aliasing energy goes into divergent modes. That overturns your
point 2 prescription (the damper is load-bearing and vorticity-neutral, not a
suppressor of the engine), while vindicating the observation underneath it. The
honest path forward is to fix the advection operator (point 1); doing so should let us
shed the intensity cap (4a) and possibly reduce the divergence damping — i.e. delete
the band-aids by removing the wound, which is the spirit of your review. Thank you for
the pass; it materially sharpened where the real work is.
