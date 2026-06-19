"""
Oracle V8 — Storm Initialization Data
=======================================
HURDAT2 best-track initialization parameters for Hugo (1989),
Katrina (2005), and Ivan (2004).

Vmax and position from:
    NOAA NHC HURDAT2 Atlantic Hurricane Database (hurdat2-1851-2024.txt)
    Available: https://www.nhc.noaa.gov/data/

Rmax from:
    RAMMB/CIRA Extended Best Track Dataset (EBTRK)
    Demuth, DeMaria, and Knaff (2006)
    Available: https://rammb2.cira.colostate.edu/research/tropical-cyclones/
               tc_extended_best_track_dataset/

DLM steering wind (u_env, v_env):
    Estimated from storm motion vector at initialization time.
    850–200 hPa deep-layer mean wind ≈ storm motion for mature TCs.
    For production paper runs, replace with ERA5 850–200 hPa mean
    from ECMWF Copernicus Climate Data Store.

Ensemble review decision (May 2026):
    Option C — gradient-wind-balanced Holland vortex + V7-compatible
    DLM steering environment. B = 1.5 frozen for all three storms.
    Parameters MUST NOT be retuned after seeing landfall errors.

HURDAT2 format reminder:
    Columns: date, time, record_id, status, lat, lon, Vmax(kt), P_min(mb),
             R34_NE, R34_SE, R34_SW, R34_NW, R50_..., R64_..., Rmax(nm)

Verify all values below against the actual HURDAT2 file before
submitting to BAMS. Values here are rounded to available precision.

HURDAT2 VERIFICATION (June 2026, hurdat2-1851-2025-02272026.txt; values
identical in the 2023 edition, so these are stable, not recent revisions):
    *** The init parameters below do NOT match HURDAT2 for any storm. ***
    HUGO    t=0: HURDAT2 27.2N 73.4W 100kt 950mb  (here: 27.7/74.5/110/942)
            Storm ID is AL111989, NOT AL131989 (fixed below).
            The position here matches HURDAT2 interpolated to ~21/0400Z;
            110 kt is the 21/1200Z value; 942 mb matches no fix.
    KATRINA t=0: HURDAT2 24.8N 85.9W 100kt 941mb  (here: 24.0/86.3/145/906)
            145 kt is the 28/1200Z value; 906 mb matches no fix
            (Katrina's minimum was 902 at 28/1800Z).
            HURDAT2 landfall: 29.3N 89.6W at 29/1110Z = t+35.17h.
    IVAN    t=0: HURDAT2 23.0N 86.0W 125kt 930mb  (here: 22.9/85.0/130/928)
            Longitude looks like an 85.0/86.0 transcription slip.
            HURDAT2 landfall: 30.2N 87.9W at 16/0650Z = t+42.83h.
    Init values are left UNCHANGED pending ensemble review (changing them
    changes every physics run).  The obs_track lists added below ARE the
    verified HURDAT2 fixes and are safe to use for verification metrics.
"""

from __future__ import annotations
import math


# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------

def kt_to_ms(kt: float) -> float:
    """Knots → m/s  (1 kt = 0.51444 m/s)."""
    return kt * 0.51444


def nm_to_m(nm: float) -> float:
    """Nautical miles → metres  (1 nm = 1852 m)."""
    return nm * 1852.0


def coriolis(lat_deg: float) -> float:
    """f = 2Ω sin(lat) at given latitude."""
    return 2.0 * 7.2921e-5 * math.sin(math.radians(lat_deg))


# ---------------------------------------------------------------------------
# Hugo — AL131989
# ---------------------------------------------------------------------------
# Initialization: September 21, 1989, 00 UTC
# ~28 h before landfall near Sullivan's Island, SC (Sep 22, ~0400 UTC)
# Storm was re-intensifying after Puerto Rico weakening, Cat 4.
#
# HURDAT2:  19890921, 0000, , HU, 27.7N, 74.5W, 110, 942
# Rmax:     Extended Best Track ~17 nm (Hugo had a compact eye)
# Landfall: 32.8°N, 79.8°W (Sullivan's Island, SC)  Cat 4, ~120 kt
# Steering: NNW at ~13 m/s inferred from 6-h storm positions

HUGO = dict(
    name           = "Hugo",
    year           = 1989,
    basin          = "AL",
    hurdat2_id     = "AL111989",   # was AL131989 — wrong; AL13 is an Oct TD

    # Initialization time
    init_date      = "1989-09-21 00Z",

    # Storm position at t=0
    lat0_deg       = 27.7,         # °N
    lon0_deg       = -74.5,        # °W

    # Holland profile parameters
    Vmax_kt        = 110,          # HURDAT2 best-track
    Vmax_ms        = kt_to_ms(110),  # 56.6 m/s
    P_min_mb       = 942,          # HURDAT2
    Rmax_nm        = 17,           # Extended Best Track (~compact eye)
    Rmax_m         = nm_to_m(17),  # 31 484 m
    B              = 1.5,          # frozen (Five/ensemble)

    # Coriolis at initialization latitude
    f              = coriolis(27.7),   # ~6.59e-5 s⁻¹

    # DLM steering (850–200 hPa deep-layer mean)
    # Hugo moving NNW at ~13 m/s at this time.
    # Decomposition: 340° true → u = -13·sin(20°), v = 13·cos(20°)
    u_env_ms       = -4.5,         # m/s  (westward component)
    v_env_ms       = +12.2,        # m/s  (northward component)

    # Observed landfall for error computation
    landfall_lat   = 32.8,         # °N
    landfall_lon   = -79.8,        # °W
    landfall_time_h = 28.0,        # hours after initialization

    # Verification (landfall_verify.py)
    threshold_lat  = 32.5,         # crossing latitude used by the runs
    # HURDAT2 best-track fixes from init time (t_h, latN, lon signed, W neg).
    # Source: hurdat2-1851-2025-02272026.txt, AL111989.  NOTE: differs from
    # the obs positions previously printed by the ERA5 pipeline — those
    # tracked ~0.3-0.6° north of HURDAT2 mid-track (provenance unknown).
    obs_track      = [
        ( 0.0,  27.2, -73.4),
        ( 6.0,  28.0, -74.9),
        (12.0,  29.0, -76.1),
        (18.0,  30.2, -77.5),
        (24.0,  31.7, -78.8),
        (28.0,  32.8, -79.8),      # landfall record, 22/0400Z
        (30.0,  33.5, -80.3),
    ],
)

# ---------------------------------------------------------------------------
# Katrina — AL122005
# ---------------------------------------------------------------------------
# Initialization: August 28, 2005, 00 UTC
# ~35 h before Gulf landfall near Buras, LA (Aug 29, ~1100 UTC)
# Katrina was near peak Cat 5 intensity (peak: 150 kt at 1800Z Aug 28)
#
# HURDAT2:  20050828, 0000, , HU, 24.0N, 86.3W, 145, 906
# Rmax:     Extended Best Track ~10 nm at peak; ~25 nm at 00Z Aug 28
#           (before the final ERC; inner eye still dominant)
# Landfall: 29.1°N, 89.6°W (near Buras, LA)  Cat 3, ~110 kt
# Steering: NNW at ~10 m/s

KATRINA = dict(
    name           = "Katrina",
    year           = 2005,
    basin          = "AL",
    hurdat2_id     = "AL122005",

    init_date      = "2005-08-28 00Z",

    lat0_deg       = 24.0,
    lon0_deg       = -86.3,

    Vmax_kt        = 145,
    Vmax_ms        = kt_to_ms(145),  # 74.6 m/s
    P_min_mb       = 906,
    Rmax_nm        = 25,           # ~25 nm pre-ERC (conservative)
    Rmax_m         = nm_to_m(25),  # 46 300 m
    B              = 1.5,

    f              = coriolis(24.0),   # ~5.94e-5 s⁻¹

    # Katrina moving NNW at ~10 m/s.
    # Decomposition: 340° true → u = -10·sin(20°), v = 10·cos(20°)
    u_env_ms       = -3.4,
    v_env_ms       = +9.4,

    landfall_lat   = 29.1,
    landfall_lon   = -89.6,
    landfall_time_h = 35.0,        # HURDAT2 L record: 29.3N 89.6W t+35.17h

    # Verification (landfall_verify.py)
    threshold_lat  = 29.1,         # crossing latitude used by the runs
    # HURDAT2 best-track fixes from init time (t_h, latN, lon signed, W neg).
    # Source: hurdat2-1851-2025-02272026.txt, AL122005.
    obs_track      = [
        ( 0.0,   24.8, -85.9),
        ( 6.0,   25.2, -86.7),
        (12.0,   25.7, -87.7),
        (18.0,   26.3, -88.6),
        (24.0,   27.2, -89.2),
        (30.0,   28.2, -89.6),
        (35.17,  29.3, -89.6),     # landfall record 1, 29/1110Z (Buras, LA)
        (36.0,   29.5, -89.6),
        (38.75,  30.2, -89.6),     # landfall record 2, 29/1445Z (LA/MS)
        (42.0,   31.1, -89.6),
    ],
)

# ---------------------------------------------------------------------------
# Ivan — AL092004
# ---------------------------------------------------------------------------
# Initialization: September 14, 2004, 12 UTC
# ~42 h before Alabama landfall near Gulf Shores (Sep 16, ~0650 UTC)
# Ivan was intensifying again after Cuba passage; Cat 4 in Gulf.
#
# HURDAT2:  20040914, 1200, , HU, 22.9N, 85.0W, 130, 928
# Rmax:     Extended Best Track ~25 nm (Ivan had a larger eye than Hugo)
# Landfall: 30.3°N, 87.7°W (Gulf Shores, AL)  Cat 3, ~105 kt
# Steering: NNW at ~8 m/s

IVAN = dict(
    name           = "Ivan",
    year           = 2004,
    basin          = "AL",
    hurdat2_id     = "AL092004",

    init_date      = "2004-09-14 12Z",

    lat0_deg       = 22.9,
    lon0_deg       = -85.0,

    Vmax_kt        = 130,
    Vmax_ms        = kt_to_ms(130),  # 66.9 m/s
    P_min_mb       = 928,
    Rmax_nm        = 25,
    Rmax_m         = nm_to_m(25),  # 46 300 m
    B              = 1.5,

    f              = coriolis(22.9),   # ~5.63e-5 s⁻¹

    # Ivan moving NNW at ~8 m/s.
    u_env_ms       = -2.7,
    v_env_ms       = +7.5,

    landfall_lat   = 30.3,
    landfall_lon   = -87.7,
    landfall_time_h = 42.0,        # HURDAT2 L record: 30.2N 87.9W t+42.83h

    # Verification (landfall_verify.py)
    threshold_lat  = 30.0,         # provisional — set before first Ivan run
    # HURDAT2 best-track fixes from init time (t_h, latN, lon signed, W neg).
    # Source: hurdat2-1851-2025-02272026.txt, AL092004.
    obs_track      = [
        ( 0.0,   23.0, -86.0),
        ( 6.0,   23.7, -86.5),
        (12.0,   24.7, -87.0),
        (18.0,   25.6, -87.4),
        (24.0,   26.7, -87.9),
        (30.0,   27.9, -88.2),
        (36.0,   28.9, -88.2),
        (42.0,   30.0, -87.9),
        (42.83,  30.2, -87.9),     # landfall record, 16/0650Z (Gulf Shores)
        (48.0,   31.4, -87.7),
    ],
)

# ---------------------------------------------------------------------------
# Convenience list
# ---------------------------------------------------------------------------

ALL_STORMS = [HUGO, KATRINA, IVAN]


# ---------------------------------------------------------------------------
# V8.6 — initialization source toggle (HURDAT2 verification, June 2026)
# ---------------------------------------------------------------------------
# The literal dicts above carry the LEGACY init values used by every run
# through V8.5.x.  HURDAT2 verification (HURDAT2_VERIFICATION.md) showed
# they do not match the best track for any storm; the corrected values
# live here and are applied when INIT_SOURCE = "hurdat2".
#
#   INIT_SOURCE = "hurdat2"  → V8.6+ production (HURDAT2-true inits)
#   INIT_SOURCE = "legacy"   → bit-reproduces all V8.5.x and earlier runs
#
# This is a DATA CORRECTION, not parameter tuning: every value below is
# read directly from hurdat2-1851-2025-02272026.txt at the init time, and
# the fallback steering (u_env_ms / v_env_ms, used only when ERA5 is
# unavailable) is the centred-difference best-track motion about t=0.
# Registered predictions for the V8.6 re-runs are in
# HURDAT2_VERIFICATION.md — read them BEFORE looking at new results.
#
# Notes:
#   * Katrina's true t=0 Vmax is 100 kt (51.4 m/s) — BELOW the 70 m/s
#     intensity cap, which therefore no longer clamps at init.
#   * f is recomputed from the true lat0.
#   * landfall_* fields update to the HURDAT2 landfall records.
#   * Rmax (EBTRK) is NOT re-verified here — separate dataset; keep the
#     existing values until checked against EBTRK directly.

INIT_SOURCE = "hurdat2"

_HURDAT2_INITS = {
    "Hugo": dict(                          # AL111989, 1989-09-21 00Z
        lat0_deg       = 27.2,
        lon0_deg       = -73.4,
        Vmax_kt        = 100,
        Vmax_ms        = kt_to_ms(100),    # 51.4 m/s
        P_min_mb       = 950,
        f              = coriolis(27.2),
        # Best-track motion, centred 20/18Z-21/06Z: WNW ~7.6 m/s
        # (legacy claimed NNW 13 m/s — also not supported by the track)
        u_env_ms       = -6.2,
        v_env_ms       = +4.4,
        landfall_lat   = 32.8,             # L record 22/0400Z = t+28.0h
        landfall_lon   = -79.8,
        landfall_time_h = 28.0,
    ),
    "Katrina": dict(                       # AL122005, 2005-08-28 00Z
        lat0_deg       = 24.8,
        lon0_deg       = -85.9,
        Vmax_kt        = 100,
        Vmax_ms        = kt_to_ms(100),    # 51.4 m/s (cap now inactive at init)
        P_min_mb       = 941,
        f              = coriolis(24.8),
        # Best-track motion, centred 27/18Z-28/06Z: WNW ~3.7 m/s
        u_env_ms       = -3.3,
        v_env_ms       = +1.8,
        landfall_lat   = 29.3,             # L record 29/1110Z (Buras, LA)
        landfall_lon   = -89.6,
        landfall_time_h = 35.17,
    ),
    "Ivan": dict(                          # AL092004, 2004-09-14 12Z
        lat0_deg       = 23.0,
        lon0_deg       = -86.0,             # legacy -85.0 was a transcription slip
        Vmax_kt        = 125,
        Vmax_ms        = kt_to_ms(125),    # 64.3 m/s
        P_min_mb       = 930,
        f              = coriolis(23.0),
        # Best-track motion, centred 14/06Z-14/18Z: NNW ~4.7 m/s
        u_env_ms       = -2.1,
        v_env_ms       = +3.4,
        landfall_lat   = 30.2,             # L record 16/0650Z (Gulf Shores)
        landfall_lon   = -87.9,
        landfall_time_h = 42.83,
    ),
}

if INIT_SOURCE == "hurdat2":
    for _s in ALL_STORMS:
        _s.update(_HURDAT2_INITS[_s["name"]])
for _s in ALL_STORMS:
    _s["init_source"] = INIT_SOURCE


def print_summary() -> None:
    """Print initialization recipe for all three storms."""
    print("Oracle V8 — Storm initialization summary")
    print(f"INIT_SOURCE = {INIT_SOURCE!r}")
    print("=" * 68)
    fmt = "{name:<10}  {init_date:<22}  {lat0_deg:>6.1f}°N  {lon0_deg:>7.1f}°W  "
    fmt2 = "Vmax={Vmax_ms:.0f}m/s  Rmax={Rmax_nm:.0f}nm  f={f:.2e}"
    for s in ALL_STORMS:
        print((fmt + fmt2).format(**s))
    print("=" * 68)
    print("B = 1.5 frozen for all storms (ensemble decision May 2026).")
    print("Verify Vmax/Rmax against HURDAT2/EBTRK before BAMS submission.")


if __name__ == "__main__":
    print_summary()
