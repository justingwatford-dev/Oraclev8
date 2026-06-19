# Oracle V8

Oracle V8 is an experimental tropical-cyclone modeling and verification workspace focused on
track-error characterization, beta-drift mechanics, and HURDAT2/ERA5-driven storm reruns.

The current codebase contains the V8 solver components, storm drivers for Hugo, Katrina, and Ivan,
beta-drift translation harnesses, verification utilities, and working notes for the publication
path. The project is in cleanup/publishability mode: the mechanism study is strong, but the
three-storm landfall-fix cluster should be treated as provisional until fresh logs are regenerated
from the cleaned scripts.

## Repository Layout

- `oracle_v8/solver/` - anelastic equation-set, operator config, Poisson/projection, and RK3 machinery.
- `oracle_v8/run_hugo.py`, `oracle_v8/run_katrina.py`, `oracle_v8/run_ivan.py` - storm rerun drivers.
- `oracle_v8/run_translation_test.py` - beta-drift and translation experiment harness.
- `oracle_v8/reanalyze_gyre.py` - tagged gyre snapshot reanalysis and figure generation.
- `oracle_v8/landfall_verify.py` - shared same-latitude and along/cross-track verification.
- `CHEATSHEET.md` - current research state, audit notes, and next-step map.
- `PAPER_track_error_characterization.md` - manuscript-oriented notes.

Generated data products are intentionally not committed: ERA5 steering `.nc` files, gyre snapshot
`.npz` files, logs, caches, and ad-hoc output are ignored by `.gitignore`.

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

Storm drivers require ERA5 steering files in `oracle_v8/`:

- `hugo_era5_steering.nc`
- `katrina_era5_steering.nc`
- `ivan_era5_steering.nc`

Those files are ignored because they are generated/local data products. Use `oracle_v8/era5_steering.py`
or the local ERA5 workflow to regenerate them.

Gyre snapshots are also generated artifacts. To reanalyze tagged snapshots:

```powershell
python oracle_v8/reanalyze_gyre.py --tag v64
```

## Current Status

The cleanup audit fixed script provenance around taper-start settings, ERA5-required storm drivers,
translation-harness `f_ref` metadata, and tagged gyre snapshot selection. The next publication-critical
work is to regenerate checked-in storm logs from the cleaned scripts, add a baseline comparison, and
expand beyond the current three-storm set.
