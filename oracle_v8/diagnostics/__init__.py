"""Oracle V8 diagnostics package.

The vortex diagnostics currently live in the legacy sibling module
``oracle_v8/diagnostics.py``.  This package wrapper re-exports that module so
``from oracle_v8 import diagnostics as dg`` works even though the diagnostics
package directory also exists.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

_legacy_path = Path(__file__).resolve().parents[1] / "diagnostics.py"
_spec = importlib.util.spec_from_file_location(
    "_oracle_v8_vortex_diagnostics",
    _legacy_path,
)
if _spec is None or _spec.loader is None:
    raise ImportError(f"Could not load vortex diagnostics from {_legacy_path}")

_vortex = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_vortex)

low_level_vmax = _vortex.low_level_vmax
vmax_3d = _vortex.vmax_3d
ventilation_flow = _vortex.ventilation_flow
vorticity_center = _vortex.vorticity_center
center_separation_km = _vortex.center_separation_km
perturbation_ke = _vortex.perturbation_ke
enstrophy = _vortex.enstrophy
compute_diagnostics = _vortex.compute_diagnostics
format_diag_line = _vortex.format_diag_line

__all__ = [
    "low_level_vmax",
    "vmax_3d",
    "ventilation_flow",
    "vorticity_center",
    "center_separation_km",
    "perturbation_ke",
    "enstrophy",
    "compute_diagnostics",
    "format_diag_line",
]
