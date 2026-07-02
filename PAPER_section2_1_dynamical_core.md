# Draft section — §2.1 Dynamical core (+ sidebar)

*Drop-in for §2.1 per [PAPER_outline.md], with the sidebar the outline flags as the BAMS-style
home for core internals. Drafted 2026-07-01 **from the code** (solver/equation_set.py,
integrator.py, tendency.py, poisson.py, grid/staggering.py, production_config.py, run_storm.py)
— every number and claim below was read from the implementation, not recalled. Author notes at
the end. Main text ≈ 420 words; sidebar ≈ 580 words (outline budgets §2 at 1200–1500 words
excluding sidebar).*

---

## 2.1 Dynamical core

Oracle V8 is a three-dimensional, nonhydrostatic dynamical core solving the Lipps and Hemler
(1982) anelastic equations — the standard soundproof system for deep convection, in which
acoustic modes are filtered by the mass-continuity constraint ∇·(ρ̄u) = 0 about a hydrostatic,
stably stratified base state. The prognostic variables are the three velocity components and the
potential-temperature perturbation θ′; pressure enters as the potential of a projection that
enforces the constraint (see sidebar). The domain is doubly periodic in the horizontal with a
rigid surface and lid, on a β-plane whose reference Coriolis parameter is set from each storm's
initialization latitude. The operating grid spacing is Δx = 15.6 km with 32 levels over a 20-km
depth, a 30-s time step, and per-storm domains of 4000–8000 km chosen by a fixed geometric rule
(Section 2.2) so that the storm's circulation never approaches the domain edges.

One configuration choice must be stated prominently, because the results of Section 4 are read
against it. For the track experiments in this paper, the thermodynamic pathway is switched off:
θ′ is zeroed at initialization and its buoyancy forcing is not applied, so temperature
perturbations are dynamically passive and the track dynamics reduce to vortex self-advection,
the β-effect, imposed environmental steering, and boundary-layer friction. In practice, the
configuration is barotropic — deliberately. Track at these scales is a steering-plus-β-drift
problem; running the core in its reduced configuration makes the β-gyre characterization of
Section 4.1 clean, and directly comparable to the idealized barotropic β-drift literature
(Chan and Williams 1987; Fiorino and Elsberry 1989) and to operational barotropic track models
of the BAM/VICBAR class. What distinguishes Oracle from those models is that the reduction is a
*configuration*, not an architecture: the same core runs fully thermodynamically active, and the
anelastic system's foundational small-perturbation assumption has been validated in that mode
(sidebar). The dry, buoyancy-off track configuration is the honest minimum that the two results
of this paper require.

The remaining ingredients are the model's dissipative and bounding operators — scale-selective
∇⁴ hyperdiffusion, a vorticity-preserving divergence damper, bulk-aerodynamic surface drag,
Newtonian cooling, and an intensity ceiling — inventoried in the sidebar. All coefficients are
fixed across every storm in this paper; none is adjusted per storm or per landfall (Section 2.2).

---

### SIDEBAR — Inside the dynamical core

**Equations and base state.** The LH82 anelastic system prognoses (u, v, w, θ′) about a dry base
state with constant buoyancy frequency N = 0.01 s⁻¹ (θ̄ = 300 K · exp(N²z/g); ρ̄ from discrete
hydrostatic integration of the Exner function). Buoyancy, when active, enters the vertical
momentum equation as b = g·θ′/θ̄. The equation set is implemented as a swappable abstraction —
the pseudo-incompressible system of Durran (1989) is the designed alternative — so that the
model can interrogate its own foundational approximation. That approximation, |θ′| ≪ θ̄, was
tested directly in a buoyancy-enabled eyewall configuration with prescribed annular heating: the
terms LH82 neglects grow linearly at ≈0.85 × θ′/θ̄ and remain at the few-percent level for
realistic eyewall amplitudes (θ′/θ̄ ≈ 4%), with stable integration to θ′/θ̄ ≈ 9% and no sign of
qualitative breakdown [ref: LH82 validation study].

**Grid and time stepping.** Vertical staggering is the Lorenz grid (u, v, θ′, and the projection
potential on level centers; w on level interfaces), with w = 0 at the rigid surface and lid.
Horizontal boundary conditions are periodic. Time integration is the three-stage Runge–Kutta of
Wicker and Skamarock (2002), with slow tendencies (advection, Coriolis, drag, diffusion,
damping, cooling) applied as Strang-split half-steps around the buoyancy-driven fast stages; in
the barotropic track configuration the fast stages are dormant and each step reduces to the
split slow dynamics. Advection is second-order centered in the advective form. The anelastic
constraint is re-enforced after every velocity-modifying substep by a projection: the
variable-coefficient Poisson problem ∇·(ρ̄∇φ) = ∇·(ρ̄u*) is solved by 2-D FFT in the horizontal
and a tridiagonal (Thomas) solve in each wavenumber column, with Neumann conditions at surface
and lid and the gauge pinned at the mean mode. The solvability (compatibility) residual of the
singular mean mode is logged at every solve and sits at machine precision across full
multi-day integrations — a running certificate that the projection is well posed, not merely
stable.

**Stabilization and bounding operators (storm-agnostic values in parentheses).**
*Hyperdiffusion:* biharmonic ∇⁴ on momentum (ν₄ = 3 × 10¹¹ m⁴ s⁻¹), scale-selective by design —
its damping time is ~20 minutes at the 2Δx scale but ~days at the vortex scale, so grid noise is
removed without eroding the storm. *Divergence damping:* the divergent part of the horizontal
flow is isolated by a Helmholtz decomposition (an FFT Poisson solve) and relaxed (ε = 0.5 per
half-step); because the correction is a pure gradient it adds exactly zero vorticity, leaving
the balanced circulation untouched. *Surface drag:* bulk-aerodynamic (C_d = 1.5 × 10⁻³),
decaying linearly to zero at H_bl = 1 km. *Newtonian cooling:* θ′ relaxation (τ = 30 min)
bounding adiabatic temperature anomalies, the standard dry-model surrogate for radiative
equilibration (Emanuel 1986). *Intensity ceiling:* perturbation winds are relaxed toward a
70 m s⁻¹ cap (τ = 300 s), a guard against numerical intensity runaway (Section 3.2, chapter 1).
The β-plane's f(y) is tapered smoothly back to f₀ over the outer 20% of the domain at the north
and south boundaries so that f is continuous across the periodic seam; the interior 60% is an
exact β-plane, and the domain rule of Section 2.2 keeps every storm inside it. Coriolis, drag,
and the intensity cap all act on the *departure* of the wind from the imposed environmental
steering flow, which is treated as a maintained geostrophic background.

---

### Author notes (delete before submission)

- **Provenance:** drafted entirely from the implementation — equation_set.py (LH82 class,
  constraint, buoyancy form), poisson.py (solver algorithm, BCs, gauge/compatibility handling),
  integrator.py (WS-RK3, Strang split; the reduced no-fast-component path), tendency.py
  (advection scheme + all operator docstrings and defaults), production_config.py (all
  coefficients), grid/staggering.py (Lorenz), run_storm.py:128 (θ′ zeroed at init). Bless the
  prose, but the facts are read, not recalled.
- **Numbers cross-check:** Δx = 15,625 m; nz = 32, Lz = 20 km (Δz = 625 m); Δt = 30 s;
  ν₄ = 3.0e11 (stability limit dx⁴/64Δt ≈ 3.1e12 — an order below); ε = 0.5; C_d = 1.5e-3;
  H_bl = 1000 m; τ_cool = 1800 s (production value — the component's own default is 3600 s;
  cite 1800); cap 70 m s⁻¹ / τ = 300 s; N = 0.01 s⁻¹; β-taper 20% per boundary → 60% exact-β
  interior; domains 4000–8000 km (nx 256–512 at fixed dx).
- **The barotropic paragraph is the reviewer-critical one.** The external review (Jun 2026)
  flagged "must justify ~barotropic vs BAM/VICBAR" — the main-text paragraph is that
  justification: reduction-as-configuration-not-architecture, plus clean comparability to the
  barotropic β-drift literature. If you want it softer or harder, that paragraph is the dial.
- **[ref: LH82 validation study]** — the eyewall assumption-validation sidebar sentence cites
  LH82_SMALL_PERTURBATION_FINDINGS.md (Phases 1–3, committed 63a74f4). Decide whether this
  becomes (a) a "data availability / supplementary" pointer, (b) a short appendix, or (c) an
  uncited "not shown." The numbers quoted (0.85×, 4%, 9%) match the findings doc.
- **Not duplicated here (belongs elsewhere):** Holland vortex + taper + the single calibrated
  parameter (§2.2); ERA5 DLM steering + lockstep (§2.3 — remember the 850–300 vs 850–200 hPa
  docstring reconciliation flagged in the cheatsheet); along/cross verification (§2.4).
- **Candidate trim if §2 runs long:** the equation-set-abstraction sentence ("swappable…
  Durran 1989") can drop to a citation; the compatibility-residual sentence can compress. Both
  are there because they serve the paper's epistemics theme (§3), not because §2.1 needs them.
- **Citations to verify vol/pages:** Lipps & Hemler (1982, JAS 39, 2192–2210 — from code
  docstring); Durran (1989, JAS 46, 1453–1461); Wicker & Skamarock (2002, MWR 130, 2088–2097);
  Lorenz (1960, Tellus 12, 364–373) if the Lorenz grid gets a citation; Emanuel (1986) for the
  Newtonian-cooling convention (the docstring cites it; confirm the intended paper is the JAS
  air–sea interaction one); Chan & Williams (1987); Fiorino & Elsberry (1989).
- **Sidebar title** is a placeholder; BAMS sidebars usually get a catchier name — e.g. "What's
  inside a 15-km hurricane core" — your call.
