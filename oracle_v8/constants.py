"""
Oracle V8 physical constants — single source of truth.

Import the named constant rather than hardcoding the literal:

    from oracle_v8.constants import GRAVITY

Keep this module a leaf (no internal imports) so anything can depend on it
without circular-import risk.

NOTE: the dry thermodynamic constants used by the base-state construction
(c_p ≈ 1004.5, R_d ≈ 287.04, p_ref = 100000) are NOT yet centralized here —
they currently live inline in production_config.build_base_state and in
validation/base_states.py. Migrating them is a separate, larger cleanup
(those values are embedded in a base-state recipe that is copy-pasted across
several test files); centralizing gravity first because it is the most widely
duplicated and the one the red-team review flagged.
"""

from __future__ import annotations

GRAVITY = 9.81  # m/s², standard gravity
