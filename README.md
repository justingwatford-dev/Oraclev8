# Oracle V8

Oracle V8 is an experimental tropical-cyclone modeling and verification workspace focused on
track-error characterization, beta-drift mechanics, and HURDAT2/ERA5-driven storm reruns.

The current codebase contains the V8 solver components, a unified storm runner for Hugo, Katrina,
and Ivan, beta-drift translation harnesses, verification utilities, and working notes for the
publication path. The project is in cleanup/publishability mode: the mechanism study is strong,
and the three-storm validation has now been regenerated under one storm-agnostic configuration.

## Repository Layout

- `oracle_v8/solver/` - anelastic equation-set, operator config, Poisson/projection, and RK3 machinery.
- `oracle_v8/run_storm.py` - unified storm runner; per-storm inputs are HURDAT2 id and init time.
- `oracle_v8/production_config.py` - single production physics/config factory shared by every storm.
- `oracle_v8/hurdat2.py` and `oracle_v8/hurdat2.txt` - HURDAT2 loader and packaged best-track data.
- `oracle_v8/run_hugo.py`, `oracle_v8/run_katrina.py`, `oracle_v8/run_ivan.py` - legacy storm-specific drivers.
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
python -B -c "import ast, pathlib; files=list(pathlib.Path('oracle_v8').rglob('*.py'))+[pathlib.Path('vortex_init.py')]; [ast.parse(p.read_text(encoding='utf-8')) for p in files]; print('parsed ok files', len(files))"
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
`run_storm.py` revalidates Hugo, Katrina, and Ivan under one storm-agnostic configuration:

| Storm | Same-latitude cross | Landfall-fix cross | Timing | Type |
| --- | ---: | ---: | ---: | --- |
| Hugo | +102.7 km | +110.2 km | -2.3 h | direct |
| Katrina | +114.3 km | +124.6 km | -2.6 h | direct |
| Ivan | +67.6 km | +126.3 km | -8.1 h | recurver |

The landfall-fix cluster survives the unified config, while the same-latitude metric is the cleaner
cross-track diagnostic for the recurving Ivan case. The next publication-critical work is a baseline
comparison and expansion beyond the current three-storm set.
