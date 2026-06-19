"""
Oracle V8 Array Backend
========================

Provides a unified `xp` alias that resolves to CuPy when a GPU is
available and NumPy otherwise.  Every compute file does:

    from oracle_v8.backend import xp

and then uses `xp.zeros`, `xp.fft.fft2`, etc. throughout.  Switching
devices requires only toggling this module — no changes to physics code.

GPU detection is automatic at import time.  To force CPU mode (e.g. for
debugging or on a machine without CUDA):

    import os
    os.environ["ORACLE_GPU"] = "0"    # set before importing oracle_v8

Helper functions
----------------
to_numpy(arr)       — convert any device array to numpy (for I/O, plots,
                       Python-level assertions)
asarray(arr)        — move array to the active device (no-op if already there)
wrap_base(base)     — return a copy of a base-state object with z, rho0,
                       theta0 on the compute device; used by RK3Integrator
                       so every tendency computation sees device arrays

V7 precedent
------------
V7 achieved a 23× speedup on the RTX 5070 using this pattern.  V8's
hot path (BuoyancyComponent, AdvectionComponent, Poisson FFT+Thomas)
maps directly to CuPy array operations and should see comparable gains.
"""

from __future__ import annotations

import os as _os

_USE_GPU = _os.environ.get("ORACLE_GPU", "1") != "0"

GPU_AVAILABLE: bool = False

if _USE_GPU:
    try:
        import cupy as xp          # type: ignore[import]
        # Importing CuPy can succeed even when the CUDA runtime/compiler DLLs
        # are missing.  Force one tiny operation now so backend selection fails
        # early and falls back to NumPy instead of crashing deep in a run.
        xp.arange(1).sum()
        GPU_AVAILABLE = True
    except Exception:
        import numpy as xp         # type: ignore[assignment]
else:
    import numpy as xp             # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Device helpers
# ---------------------------------------------------------------------------

def to_numpy(arr):
    """
    Return `arr` as a NumPy array.

    If GPU_AVAILABLE and arr is a CuPy array, copies it to CPU first.
    If arr is already NumPy (or GPU not available), returns it unchanged.
    Used for diagnostics, file I/O, and test assertions that need
    Python-level values.
    """
    if GPU_AVAILABLE:
        try:
            import cupy  # type: ignore[import]
            if isinstance(arr, cupy.ndarray):
                return cupy.asnumpy(arr)
        except ImportError:
            pass
    return arr


def asarray(arr, dtype=None):
    """
    Move `arr` to the active compute device.

    Wraps xp.asarray.  No-op if arr is already on the correct device.
    If dtype is specified, also casts.
    """
    if dtype is not None:
        return xp.asarray(arr, dtype=dtype)
    return xp.asarray(arr)


def wrap_base(base):
    """
    Return a device-resident copy of a base-state object.

    Moves the canonical array attributes (z, rho0, theta0) to the
    compute device.  Non-array attributes (integration_scheme, etc.)
    are preserved as-is.

    Called once at RK3Integrator construction so tendency components
    always receive device arrays from `base`.

    Parameters
    ----------
    base : any object with z, rho0, theta0 array attributes

    Returns
    -------
    A lightweight object with the same interface as `base` but with
    array fields on the active device.
    """
    class _DeviceBase:
        pass

    dev = _DeviceBase()

    # Move the three canonical array fields to device
    for attr in ("z", "rho0", "theta0"):
        if hasattr(base, attr):
            setattr(dev, attr, asarray(getattr(base, attr)))

    # Preserve all other non-private attributes — both non-callable values
    # AND callable methods (e.g. is_load_bearing_ready).
    # wrap_base creates a plain object that does not inherit from the original
    # class, so class-level methods become invisible unless explicitly forwarded.
    # LH82AnelasticEquationSet.base_state_compatibility() calls
    # base.is_load_bearing_ready() and returns (False, ...) if the method is
    # absent — silently degrading the compatibility guard in any code path
    # that passes a wrapped base to the equation set.
    # Callables are forwarded as-is (no device transfer needed for Python methods).
    # Fix: Opus 4.6 red-team audit C8, May 2026.
    for attr in dir(base):
        if attr.startswith("_") or attr in ("z", "rho0", "theta0"):
            continue
        val = getattr(base, attr, None)
        if val is not None:
            setattr(dev, attr, val)   # scalars and callables forwarded as-is

    return dev
