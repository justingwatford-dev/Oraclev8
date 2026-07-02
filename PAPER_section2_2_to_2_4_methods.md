# Draft sections — §2.2–2.4 (initialization, steering, verification)

*Drop-ins for §2.2, §2.3, §2.4 per [PAPER_outline.md]. Drafted 2026-07-01 from the code
(hurdat2.py, vortex_init.py, production_config.py, era5_steering.py, landfall_verify.py,
run_storm.py) and the verified docs. §2.2 is the resolution target for the Results subsection's
"storm-agnostic configuration (Section X)" cross-ref (#3 in the outline map). Author notes at
the end. Combined ≈ 960 words; with §2.1's ≈ 420 the §2 total lands ≈ 1380, inside the outline's
1200–1500 budget (sidebar excluded).*

---

## 2.2 Vortex initialization and the single structural parameter

Every storm is initialized directly from the HURDAT2 best-track file at the chosen
initialization time. Position, maximum wind, central pressure, and the Coriolis parameter (set
from the initialization latitude) are read from the t = 0 fix by the loader — no initialization
value is hand-entered — and the verification reference is simply every subsequent best-track fix
(Section 2.4). The file itself, and our reading of it, were verified against the current and
prior HURDAT2 editions; Section 3.2 recounts what that verification caught.

The initial vortex is a Holland (1980) gradient-wind-balanced cyclone,
V(r) = V_max (R_max/r)^B exp[1 − (R_max/r)^B], with the observed maximum wind, a shape parameter
B = 1.5 frozen across all storms, and R_max fixed at 75 km for every storm — a documented
resolution floor (≈5Δx at the operating grid spacing), not observed inner-core structure. The
profile decays exponentially with height from its surface reference; pressure follows from
gradient-wind balance, integrated inward from the environmental radius. (The associated
warm-core temperature perturbation is then zeroed, per the barotropic configuration of
Section 2.1.) Five projection-only iterations remove residual initialization divergence before
the integration clock starts.

Left untreated, the Holland profile decays so slowly with radius (r^{−B/2}) that the circulation
fills any domain, and the storm's β-drift inflates accordingly (Section 3.2, chapter 3). The
outer wind is therefore tapered to zero by a cosine ramp beginning at the *taper-onset radius*
and completing at the environmental radius R_env = 500 km. The taper-onset radius is the
configuration's **single structural free parameter** — and the sweeps show the two nominal knobs
are degenerate (only the product of R_env and the onset fraction matters), so it is genuinely
one parameter, not two. It was calibrated in the quiescent-environment testbed of Section 4.1:
sweeping the onset radius moves the simulated β-drift magnitude and aim together, and an onset
radius of 200 km (≈2.7 R_max) brings the drift magnitude into the canonical published range
while minimizing the poleward aim bias. That value was locked before any production storm run
and shared by every storm in this paper. We emphasize the resulting property on which the
"blind" in blind track skill rests: **no parameter anywhere in the system has ever been adjusted
against a landfall.** The calibration target is idealized-vortex physics from the literature,
and no landfall-fitting pathway exists in the code.

Domain size and run length are likewise derived, not chosen: a fixed geometric rule takes each
storm's initialization and threshold latitudes and returns the smallest standard domain
(4000–8000 km; 256–512 grid points at fixed Δx) that keeps the storm's full circulation, with
margin, inside the exact-β interior of Section 2.1's tapered β-plane; integrations run 10 h past
the observed landfall time.

## 2.3 Environmental steering

The environment is imposed as a horizontally uniform background flow relaxed toward the ERA5
deep-layer mean (DLM) — the mass-weighted 850/700/500/300-hPa average, the standard
tropical-cyclone steering layer — computed on the raw ERA5 grid and averaged over a 3–7°
annulus centered on the *model* storm, so that the inner core's own circulation does not steer
it. The annulus is sampled at the model's position and interpolated in time: the storm feels the
environment where it actually is, not where the observed storm was. Every 30 minutes of
simulation the background flow relaxes toward the locally sampled DLM with a 3-h timescale; each
increment is applied simultaneously to the flow carried in the model state and to the
perturbation references of the Coriolis, drag, and intensity-cap operators (the
perturbation-relative convention of Section 2.1), so the forcing operators and the flow they
reference advance in lockstep. Because the increment is horizontally uniform it carries no
divergence and no vorticity: the projection, and the vortex itself, are untouched by
construction. Time-varying steering is architecture here, not refinement: in a configuration
accident that became an informative A/B test, steering frozen at its initial value failed in
*opposite* directions on different storms (+182 km east on Hugo, −192 km west on Ivan), which no
steering-independent mechanism can produce. One honest limitation is deferred to Section 5: the
annulus average is not a storm-removed environmental field, and at these radii it retains some
storm-induced flow.

## 2.4 Track verification

All track scoring passes through one shared verification module, against the HURDAT2 fixes, in
two layers. The first is the headline same-latitude threshold comparison: the model's and the
observed track's first northward crossings of a threshold latitude (the observed landfall
latitude) are interpolated, and the crossing-time difference and the longitude offset at
crossing separate timing from cross-track displacement. The second is the full along/cross
decomposition evaluated at every observed fix: the along-track component is the model's
displacement ahead of (or behind) the observed storm along its instantaneous direction of
motion, and the cross-track component is the displacement to the right of that motion —
approximately eastward for the poleward-moving storms considered here. Two retired metrics
motivate this convention. The legacy landfall-point comparison — model threshold crossing versus
the observed landfall point — folds the observed storm's remaining along-track travel into
apparent cross-track error and overstated one re-scored run by half. And any single cross-track
scalar can manufacture error for a storm that is not moving due north at threshold: Michael's
poleward along-track overshoot against an observed track recurving northeast projects onto the
same-latitude axis as a spurious *westward* miss, which the full decomposition correctly
renders as timing (Sections 3.3, 4.2). Section 4 therefore reports both layers for all six
storms.

---

### Author notes (delete before submission)

- **§2.2 is cross-ref target #3**: the Results subsection's "storm-agnostic configuration
  (Section X)" resolves to §2.2 — the bolded no-landfall-adjustment sentence is the anchor.
- **Numbers/facts provenance:** B = 1.5 and R_max = 75 km frozen in hurdat2.storm_init defaults +
  production_config (RMAX_RUN_M, "5×dx"); resolution-floor honesty per red-team item 18.
  Taper-onset 200 km / R_env 500 km / degeneracy / "locked before the runs" = cheatsheet β-DRIFT
  CALIBRATION (gate-beta-renv/-taper: (400, 0.5) ≡ (650, 0.30) at taper-start ≈195–200 km);
  200/75 ≈ 2.7. Domain rule = production_config.choose_domain (STANDARD_NX 256–512,
  N_STEPS_MARGIN_H = 10). Pre-balance = N_PREBAL = 5. DLM layer/weights + 3–7° annulus =
  era5_steering.py (_LEVELS_HPA + get_dlm defaults). τ_steer = 10800 s; update cadence
  DIAG_EVERY·DT = 30 min. Frozen-steering A/B (+182 Hugo reval-1 / −192 Ivan run-1) =
  cheatsheet storm history + HURDAT2_VERIFICATION.md (Ivan run-1 section). "Overstated one
  re-scored run by half" = Hugo re-val 2: +137 km landfall-point vs +91 km same-latitude
  (landfall_verify.py docstring — note those were vs the pipeline track pre-correction; if you'd
  rather quote a HURDAT2-scored pair, swap in the +128.3 legacy vs +102.7 same-lat pair from the
  corrected Hugo run in HURDAT2_VERIFICATION.md).
- **H_wind (wind-profile scale height):** the exponential vertical decay constant in
  vortex_init.wind_vertical_structure — value not quoted in the draft ("decays exponentially
  with height"); insert the default if you want the number in print. [check vortex_init.py]
- **"Canonical published range"** in §2.2 deliberately avoids re-quoting numbers — §4.1 owns the
  1–3 m s⁻¹ band and its citations; the testbed acceptance band was 1.5–2.5 m s⁻¹ (gate-beta) if
  a reviewer asks for the calibration criterion precisely.
- **§2.3 annulus caveat** is red-team hold item 2 (annulus ≠ storm-removed field; DLM
  sensitivity test outstanding) — the §5 draft must actually pick this up.
- **DLM layer reconciliation:** code = 850–300 hPa (four levels, layer-thickness weights);
  storm_data docstring says 850–200. The draft states 850–300 per the implementation — fix the
  docstring or add an erratum note so the paper and code agree.
- **Michael metric sentence** — verify the recurve direction wording against the six-storm
  note ("obs recurves NE at landfall") and the track plot (the note's caveat asks for an
  eyeball check of Michael's steering before this goes to print).
