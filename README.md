# Oracle V8

Oracle V8 is a research-grade, from-scratch **anelastic tropical-cyclone dynamical core**
(Lipps–Hemler 1982) plus a HURDAT2/ERA5-driven verification workspace. It is built to study
*why* an idealized TC model tracks the way it does — beta-drift mechanics, track-error structure,
and (in progress) the eyewall thermodynamics where the small-perturbation assumption is tested.

The model runs from a **single storm-agnostic configuration**: the only per-storm inputs are a
HURDAT2 id and an initialization time; the domain, run length, and physics are derived/fixed.
It runs on CPU (NumPy) or GPU (CuPy) through one backend shim.

> Status: active research code, not a packaged product. Interfaces move, and some abstraction
> layers are intentionally thin (see *Honesty notes* below). Reproducibility of the headline
> runs is the priority, not API stability.

## Where things stand

The project currently rests on **two independent results**, plus an active solver-hardening and
buoyancy-study effort:

### 1. Storm-agnostic landfall reproduction + a verification methodology

One shared configuration reproduces the US landfalls of multiple historical storms from
first-principles dynamics. The current validation set is six storms; same-latitude cross-track:

| Storm   | Same-latitude cross-track | Type     |
| ------- | ------------------------- | -------- |
| Hugo    | +102.7 km (E)             | direct   |
| Katrina | +114.3 km (E)             | direct   |
| Ivan    | +67.6 km (E)              | recurver |
| Fran    | +23.0 km (E)              | direct   |
| Michael | −75.0 km (W)              | direct   |
| Laura   | −38.7 km (W)              | direct   |

**Honest caveat:** the tight *eastward* cluster seen in the original calibration storms
(Hugo/Katrina/Ivan) does **not** generalize out of sample — Michael and Laura miss to the
**west**. So the eastward bias is a property of the calibration set, not a universal model bias.
A separate strand of the work is the *verification methodology* itself: a HURDAT2 audit that
caught init/“observed-track” provenance errors (recalled-from-memory values vs. read-from-file),
and same-latitude along/cross-track scoring via `landfall_verify`.

### 2. A characterized model property: the β-gyre over-rotates poleward

Idealized β-drift experiments show the model's **β-gyre equilibrates ~20° too far poleward**
(NNW instead of the canonical NW). The effect is **structural and intensity-independent** —
flat across a 3× range of vortex strength, ruling out a swirl-shear origin — and the drift
*speed* saturates at the correct β-drift magnitude (~2.5 m/s) while its *direction* precesses.
An f-plane null (zero drift with β off) confirms it is genuine β-drift, not a numerical artifact.

![β-gyre m=1 asymmetry precessing](gyre_precession.png)

*β-gyre m=1 vorticity asymmetry (azimuthal mean removed, re-centered), t = 12–60 h. The black
arrow is the swirl-removed steering flow, which matches the simulated β-drift.*

**What is NOT claimed:** earlier notes tied this poleward gyre bias *quantitatively* to the
real-storm cross-track error (the “landfall bridge”). That causal claim has been **retired** —
a predicted along/cross projection did not hold up. The β-gyre over-rotation stands on its own as
a characterized model property; the link to landfall error is an open question, not a result.

### 3. (Active) Solver hardening + a buoyancy/eyewall study

The solver core recently went through a thorough, file-by-file **external red-team review**. That
pass and the follow-on work produced, among other things:

- A fix to the anelastic Poisson **zero-mode** (solve it instead of discarding it; negligible for
  the barotropic track runs, ~12% of `max|w|` once buoyancy is on).
- A **single source of truth** for the discrete divergence operator, so the validation harness
  measures the exact constraint the solver enforces (it previously disagreed ~48×).
- An honest downgrade of the staggering abstraction to what it actually implements (Lorenz).
- A 5th-order **upwind advection** option (Wicker–Skamarock 2002): the root cause of the
  barotropic instability the old config managed with a divergence damper was the centered
  advection scheme's lack of dissipation.
- A prescribed **eyewall-heating driver** for the (otherwise inert) dry buoyancy study. With a
  diabatic driver present, upwind advection + a light/zero divergence damper gives both stability
  *and* a physical secondary circulation — toward the original goal of probing where the LH82
  small-perturbation assumption (|θ′| ≪ θ̄) breaks.

These changes live on the `red-team-fixes-and-buoyancy-driver` branch / recent history. The
point-by-point review responses are in `RED_TEAM_RESPONSE_*.md`.

## Repository layout

- `oracle_v8/solver/` — the dynamical core: `equation_set.py` (LH82 anelastic), `tendency.py`
  (advection, Coriolis, buoyancy, diabatic heating, divergence damping, projection components),
  `poisson.py` (variable-coefficient anelastic Poisson solver, FFT + batched Thomas),
  `integrator.py` (split-explicit RK3), `operator_config.py` (composable run configuration).
- `oracle_v8/grid/staggering.py` — Lorenz vertical staggering (the live seam).
- `oracle_v8/vortex_init.py` — Holland (1980) gradient-wind-balanced vortex; θ′ derived from the
  wind via hydrostatic balance.
- `oracle_v8/run_storm.py` — unified storm runner (per-storm input = HURDAT2 id + init time).
- `oracle_v8/production_config.py` — the single storm-agnostic physics/config factory.
- `oracle_v8/run_translation_test.py` — idealized β-drift / translation harnesses.
- `oracle_v8/reanalyze_gyre.py` — tagged β-gyre snapshot reanalysis and figure generation.
- `oracle_v8/landfall_verify.py` — same-latitude along/cross-track verification.
- `oracle_v8/hurdat2.py`, `hurdat2.txt` — HURDAT2 loader and best-track data.
- `oracle_v8/validation/` — analytical-solution tests, manufactured-solution convergence, and the
  shared test harness.
- `oracle_v8/legacy/` — deprecated per-storm drivers, retained for provenance.
- `CHEATSHEET.md` — the detailed running research log (build state, audit notes, next steps).

Generated data products (ERA5 steering `.nc`, gyre snapshot `.npz`, caches, ad-hoc output) are
ignored by `.gitignore`; the agnostic validation logs are committed as evidence.

## Setup

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

For GPU runs, install the CuPy build matching the local CUDA version (the model auto-selects
NumPy if CuPy/CUDA is unavailable). Force CPU with `ORACLE_GPU=0`.

## Quick checks

```powershell
# validation tests (analytical solutions, manufactured-solution convergence, component checks)
python -m oracle_v8.validation.tests.test_poisson_solver
python -m oracle_v8.validation.tests.test_buoyancy_component

# import / syntax sanity over the package
python -B -c "import ast, pathlib; files=list(pathlib.Path('oracle_v8').rglob('*.py')); [ast.parse(p.read_text(encoding='utf-8')) for p in files]; print('parsed', len(files), 'files ok')"
```

## Running a storm

Unified runs read HURDAT2 from `oracle_v8/hurdat2.txt` and require an ERA5 DLM steering file
(`<storm>_era5_steering.nc`, a local/generated product — regenerate with
`oracle_v8/era5_steering.py`). Missing ERA5 aborts rather than silently falling back.

```powershell
python -m oracle_v8.run_storm Ivan
python -m oracle_v8.run_storm Hugo
python -m oracle_v8.run_storm Katrina
```

β-gyre reanalysis from saved snapshots:

```powershell
python oracle_v8/reanalyze_gyre.py --tag v64
```

## Honesty notes

This is research code under active review, and a few things are deliberately partial:

- The `EquationSet`, `GridStaggering`, and `PseudoIncompressible` abstractions are thin —
  several methods are placeholders, and the production path computes some stencils inline. The
  docstrings now say so rather than overstating a “drop-in swap.”
- The production config is **barotropic** (buoyancy off, θ′ a passive tracer); the buoyancy
  pathway is exercised only in the study runs.
- A running research log with full provenance — including where earlier claims were walked
  back — lives in `CHEATSHEET.md`. When this README and the cheat sheet disagree, the cheat
  sheet is newer.
