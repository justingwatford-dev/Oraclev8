"""
Oracle V8 physical constants — single source of truth.
=======================================================

Centralizes the dry-air thermodynamic constants and standard gravity that
were previously duplicated as inline literals across the solver, vortex init,
storm tracker, and base-state builders. Import the named constant rather than
hardcoding the literal:

    from oracle_v8.constants import GRAVITY, C_P, R_D, P_REF

Keep this module a leaf (no internal imports) so anything can depend on it
without circular-import risk.

These values are deliberately the *exact* literals the codebase has always
used, not their more precise SI counterparts (e.g. GRAVITY is 9.81, not
9.80665). Changing any value here changes physics results bit-for-bit, so
treat edits as a model change, not a refactor.

Sources:
    - C_P, R_D: AMS Glossary of Meteorology ("specific heat", "gas constant")
    - GRAVITY, P_REF: WMO standard atmosphere
"""

GRAVITY = 9.81        # m/s², standard gravity
R_D = 287.04          # J/(kg·K), gas constant for dry air
C_P = 1004.5          # J/(kg·K), specific heat of dry air at constant pressure
P_REF = 100_000.0     # Pa, reference pressure for the Exner/Π definition
