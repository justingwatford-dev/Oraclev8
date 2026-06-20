# Oracle V8

Oracle V8 is an experimental tropical-cyclone modeling and verification workspace focused on
track-error characterization, beta-drift mechanics, and HURDAT2/ERA5-driven storm reruns.

The current codebase contains the V8 solver components, a unified storm runner for Hugo, Katrina,
and Ivan, beta-drift translation harnesses, verification utilities, and working notes for the
publication path. The project is in cleanup/publishability mode: the mechanism study is strong,
and the three-storm validation has now been regenerated under one storm-agnostic configuration.

## Key Result

Oracle reproduces the US landfalls of Hugo (1989), Katrina (2005), and Ivan (2004) from
first-principles physics under a **single storm-agnostic configuration** — the only per-storm
inputs are a HURDAT2 id and an initialization time.

All three share a systematic **eastward cross-track bias**. On the clean same-latitude metric the
two direct storms cluster tightly (the recurving Ivan is the understood outlier — the model misses
its zonal stall, so it crosses the scoring latitude early):

| Storm   | Same-latitude cross-track | Timing | Type     |
| ------- | ------------------------- | ------ | -------- |
| Hugo    | +102.7 km                 | −2.3 h | direct   |
| Katrina | +114.3 km                 | −2.6 h | direct   |
| Ivan    | +67.6 km                  | −8.1 h | recurver |

Idealized β-drift experiments locate the mechanism: the model's **β-gyre equilibrates ~20° too far
poleward** (NNW instead of canonical NW). The bias is *structural and intensity-independent* —
confirmed flat across a 3× range of vortex strength, ruling out a swirl-shear origin. A poleward
rotation of the few-m/s β-drift accumulated over a multi-day track is consistent in magnitude with
the ~+100 km eastward cross-track seen in the two direct storms, tying the idealized mechanism to the
real-storm error.

![β-gyre m=1 asymmetry precessing](gyre_precession.png)

*β-gyre m=1 vorticity asymmetry (azimuthal mean removed, re-centered) at t = 12–60 h. The asymmetry
intensifies as the storm drifts; the black arrow is the swirl-removed steering flow, which matches
the simulated β-drift. The gyre amplitude saturates while its orientation holds poleward of the
canonical NW — the structural source of the eastward track bias.*

## Repository Layout

- `oracle_v8/solver/` - anelastic equation-set, operator config, Poisson/projection, and RK3 machinery.
- `oracle_v8/run_storm.py` - unified storm runner; per-storm inputs are HURDAT2 id and init time.
- `oracle_v8/production_config.py` - single production physics/config factory shared by every storm.
- `oracle_v8/hurdat2.py` and `oracle_v8/hurdat2.txt` - HURDAT2 loader and packaged best-track data.
- `oracle_v8/legacy/` - deprecated per-storm drivers (`run_hugo.py`, `run_katrina.py`, `run_ivan.py`); superseded by `run_storm.py`, retained for provenance and gated behind `ORACLE_ALLOW_LEGACY=1`.
- `oracle_v8/run_translation_test.py` - beta-drift and translation experiment harness.
- `oracle_v8/reanalyze_gyre.py` - tagged gyre snapshot reanalysis and figure generation.
- `oracle_v8/landfall_verify.py` - shared same-latitude and along/cross-track verification.
- `CHEATSHEET.md` - current research state, audit notes, and next-step map.
- `PAPER_track_error_characterization.md` - manuscript-oriented notes.

Generated data products are intentionally not committed: ERA5 steering `.nc` files, gyre snapshot
`.npz` files, caches, and ad-hoc output are ignored by `.gitignore`. The three agnostic validation
logs are committed as evidence for the current storm-agnostic claim.

## Setup

Create and activate a virtual environment, then install dependencies:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

For GPU runs, install the CuPy package matching the local CUDA version. See
`oracle_v8/requirements.txt` for the optional CuPy lines.

## Quick Checks

Run the lightweight test suite:

```powershell
python -m pytest oracle_v8
```

Run a syntax/import sanity check:

```powershell
python -B -c "import ast, pathlib; files=list(pathlib.Path('oracle_v8').rglob('*.py')); [ast.parse(p.read_text(encoding='utf-8')) for p in files]; print('parsed ok files', len(files))"
```

## Data Notes

Unified storm runs read HURDAT2 data from `oracle_v8/hurdat2.txt`. Storm drivers also require ERA5
steering files in `oracle_v8/`:

- `hugo_era5_steering.nc`
- `katrina_era5_steering.nc`
- `ivan_era5_steering.nc`

Those files are ignored because they are generated/local data products. Use `oracle_v8/era5_steering.py`
or the local ERA5 workflow to regenerate them.

Run the unified driver from the repository root:

```powershell
python -m oracle_v8.run_storm Hugo
python -m oracle_v8.run_storm Katrina
python -m oracle_v8.run_storm Ivan
```

Gyre snapshots are also generated artifacts. To reanalyze tagged snapshots:

```powershell
python oracle_v8/reanalyze_gyre.py --tag v64
```

## Current Status

The cleanup audit fixed script provenance around taper-start settings, ERA5-required storm drivers,
translation-harness `f_ref` metadata, tagged gyre snapshot selection, and per-storm config drift.
`run_storm.py` revalidates Hugo, Katrina, and Ivan under one storm-agnostic configuration (see
**Key Result** above for the same-latitude cross-track headline). The landfall-fix cross-track
cluster — Hugo +110.2 km, Katrina +124.6 km, Ivan +126.3 km (mean ~+120 km) — survives the unified
config, while the same-latitude metric is the cleaner cross-track diagnostic for the recurving Ivan
case. The next publication-critical work is a baseline comparison and expansion beyond the current
three-storm set.
