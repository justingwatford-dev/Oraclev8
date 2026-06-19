"""
Oracle V8 — ERA5 Deep Layer Mean Steering Winds
================================================
Replaces the constant u_env / v_env in run_hugo.py with time- and
position-varying DLM winds from ERA5 reanalysis.

What is the DLM?
----------------
The Deep Layer Mean (850–300 hPa mass-weighted average) is the standard
NHC steering layer.  Replacing the fixed constant with ERA5 DLM is
expected to reduce the 325 km east bias in Run 16 by capturing Hugo's
actual recurvature toward the NW as it approached the SC coast.

Setup (one-time)
----------------
1. Create a free account at https://cds.climate.copernicus.eu
2. Install the API client:
       conda activate woe_env
       pip install cdsapi
3. Create ~/.cdsapirc with your UID and API key:
       url: https://cds.climate.copernicus.eu/api/v2
       key: <UID>:<API-key>
4. Run the download helper once:
       python -m oracle_v8.era5_steering --download
   This saves hugo_era5_steering.nc to the oracle_v8/ directory (~10 MB).

Usage in run_hugo.py
--------------------
    from oracle_v8.era5_steering import ERA5Steering

    steering = ERA5Steering.load('oracle_v8/hugo_era5_steering.nc')

    # Inside the main loop, after tracking:
    t_hours = step * dt / 3600.0
    u_env, v_env = steering.get_dlm(t_hours, lat_c, lon_c)
    config.coriolis.u_env = u_env   # works if CoriolisComponent is not frozen
    config.coriolis.v_env = v_env

ERA5 request spec
-----------------
  Dataset : reanalysis-era5-pressure-levels
  Variable: u_component_of_wind, v_component_of_wind
  Levels  : 300, 500, 700, 850 hPa
  Dates   : 1989-09-21 to 1989-09-22
  Times   : 00/03/06/09/12/15/18/21 UTC  (3-hourly)
  Area    : N=45, W=-92, S=15, E=-62  (covers the 2000×2000 km domain)
  Format  : NetCDF
"""

from __future__ import annotations

import argparse
import os
import sys
import numpy as np

# ── optional heavy imports (not needed at import time) ────────────────────────
try:
    import netCDF4 as nc
    _HAS_NC4 = True
except ImportError:
    _HAS_NC4 = False

try:
    from scipy.interpolate import RegularGridInterpolator
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


# ─── constants ────────────────────────────────────────────────────────────────

# ---------------------------------------------------------------------------
# Per-storm configuration
# ---------------------------------------------------------------------------
# Each entry contains everything ERA5Steering needs that differs by storm:
#   init_*      : initialisation datetime (Hugo: Sep 21 1989 00Z, etc.)
#   area        : CDS API bounding box [N, W, S, E] in degrees (W/E as °E)
#   days        : calendar days to request (covers init + 2 days)
#   nc_path     : where to save/load the NetCDF file
#   obs_track   : {t_hours: (lat, lon_W)} HURDAT2 best-track fixes.
#                 Used ONLY by print_summary (sanity table), the --levels
#                 CLI, and get_track_mean_dlm (frozen-steering mode,
#                 retired since V8.4.0).  The time-varying lockstep
#                 pathway samples get_dlm at the MODEL position and
#                 never touches obs_track.
#
# Katrina ERA5 area centred at 24°N, 86.3°W, 2000km domain →
#   [40, -105, 12, -70] gives comfortable margin.

STORM_CONFIGS = {
    "hugo": dict(
        init_year  = 1989,
        init_month = 9,
        init_day   = 21,
        init_hour  = 0,
        area       = [45, -92, 15, -62],   # [N, W, S, E] in °N / °E
        year_str   = "1989",
        month_str  = "09",
        days       = ["21", "22"],
        nc_path    = os.path.join(os.path.dirname(__file__),
                                  "hugo_era5_steering.nc"),
        # HURDAT2 best-track fixes (AL111989, hurdat2-1851-2025 ed.).
        # ⚠ The previous values here were NOT HURDAT2 (synthetic track,
        # ~0.3-0.6° north of best track mid-run) — see HURDAT2_VERIFICATION.md.
        obs_track  = {
             0: (27.2, 73.4),
             6: (28.0, 74.9),
            12: (29.0, 76.1),
            18: (30.2, 77.5),
            24: (31.7, 78.8),
            28: (32.8, 79.8),   # landfall record, 22/0400Z
        },
    ),
    "katrina": dict(
        init_year  = 2005,
        init_month = 8,
        init_day   = 28,
        init_hour  = 0,
        area       = [40, -105, 12, -70],  # Gulf domain, 2000km margin
        year_str   = "2005",
        month_str  = "08",
        days       = ["28", "29", "30"],   # covers 32h run + buffer
        nc_path    = os.path.join(os.path.dirname(__file__),
                                  "katrina_era5_steering.nc"),
        # HURDAT2 best-track fixes (AL122005, hurdat2-1851-2025 ed.).
        # ⚠ The previous values here were NOT HURDAT2: they were a straight
        # line from the (incorrect) init to the landfall point — every fix
        # matches linear interpolation exactly.  See HURDAT2_VERIFICATION.md.
        obs_track  = {
             0:     (24.8, 85.9),
             6:     (25.2, 86.7),
            12:     (25.7, 87.7),
            18:     (26.3, 88.6),
            24:     (27.2, 89.2),
            30:     (28.2, 89.6),
            35.17:  (29.3, 89.6),   # landfall record 1, 29/1110Z (Buras, LA)
        },
    ),
    "ivan": dict(
        init_year  = 2004,
        init_month = 9,
        init_day   = 14,
        init_hour  = 12,                   # ⚠ 12Z init — only storm not at 00Z
        area       = [42, -105, 12, -70],  # Gulf + recurvature headroom to 42°N
        year_str   = "2004",
        month_str  = "09",
        days       = ["14", "15", "16"],   # t = −12h … +57h relative to init
        nc_path    = os.path.join(os.path.dirname(__file__),
                                  "ivan_era5_steering.nc"),
        # HURDAT2 best-track fixes (AL092004, hurdat2-1851-2025 ed.).
        obs_track  = {
             0:     (23.0, 86.0),
             6:     (23.7, 86.5),
            12:     (24.7, 87.0),
            18:     (25.6, 87.4),
            24:     (26.7, 87.9),
            30:     (27.9, 88.2),
            36:     (28.9, 88.2),
            42:     (30.0, 87.9),
            42.83:  (30.2, 87.9),   # landfall record, 16/0650Z (Gulf Shores)
        },
    ),
}

_DEFAULT_STORM = "hugo"

# ---------------------------------------------------------------------------
# DLM pressure-level weights (shared across all storms)
# ---------------------------------------------------------------------------
_LEVELS_HPA = np.array([850., 700., 500., 300.])
_DLM_WEIGHTS = np.array([
    (850. - 700.) / 2.,
    (850. - 700.) / 2. + (700. - 500.) / 2.,
    (700. - 500.) / 2. + (500. - 300.) / 2.,
    (500. - 300.) / 2.,
], dtype=float)
_DLM_WEIGHTS /= _DLM_WEIGHTS.sum()   # [0.136, 0.318, 0.364, 0.182]


# ─── downloader ──────────────────────────────────────────────────────────────

def download_era5(output_path: str = None, storm: str = _DEFAULT_STORM) -> None:
    """
    Download ERA5 pressure-level winds for the given storm via the CDS API.
    Requires a configured ~/.cdsapirc (see module docstring).

    Parameters
    ----------
    output_path : override save path (default: storm-specific path in STORM_CONFIGS)
    storm       : 'hugo' or 'katrina' (default: 'hugo')
    """
    cfg = STORM_CONFIGS[storm.lower()]
    if output_path is None:
        output_path = cfg["nc_path"]
    try:
        import cdsapi
    except ImportError:
        sys.exit(
            "cdsapi not found.  Install with:\n"
            "  conda activate woe_env && pip install cdsapi"
        )

    print(f"Downloading ERA5 steering winds for {storm.capitalize()} "
          f"→ {output_path}")
    print("  (this may take a few minutes depending on CDS queue depth)")

    c = cdsapi.Client()
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'variable'    : ['u_component_of_wind', 'v_component_of_wind'],
            'pressure_level': [str(int(p)) for p in _LEVELS_HPA],
            'year'  : cfg["year_str"],
            'month' : cfg["month_str"],
            'day'   : cfg["days"],
            'time'  : [f'{h:02d}:00' for h in range(0, 24, 3)],
            'area'  : cfg["area"],
            'format': 'netcdf',
        },
        output_path,
    )
    print(f"  Saved → {output_path}")


# ─── loader / interpolator ───────────────────────────────────────────────────

class ERA5Steering:
    """
    Load and interpolate ERA5 DLM steering winds.

    Attributes
    ----------
    times_h   : 1-D array, hours since Hugo init (t=0 = 1989-09-21 00Z)
    lats      : 1-D array, latitude grid (°N, ascending)
    lons      : 1-D array, longitude grid (°W, ascending → stored as negative °)
    u_dlm     : 3-D array (time, lat, lon), DLM zonal wind (m/s)
    v_dlm     : 3-D array (time, lat, lon), DLM meridional wind (m/s)
    """

    def __init__(self, times_h, lats, lons, u_dlm, v_dlm):
        if not _HAS_SCIPY:
            raise ImportError(
                "scipy required for ERA5Steering interpolation.\n"
                "  conda activate woe_env && conda install scipy"
            )
        self.times_h = times_h
        self.lats    = lats
        self.lons    = lons   # in °W (positive values, e.g. 74.5 not -74.5)
        self.u_dlm   = u_dlm
        self.v_dlm   = v_dlm

        # Build bilinear interpolators (time, lat, lon)
        # scipy wants axes in ascending order
        self._u_itp = RegularGridInterpolator(
            (times_h, lats, lons), u_dlm, method='linear',
            bounds_error=False, fill_value=None,
        )
        self._v_itp = RegularGridInterpolator(
            (times_h, lats, lons), v_dlm, method='linear',
            bounds_error=False, fill_value=None,
        )

    # ── constructor ──────────────────────────────────────────────────────────

    @classmethod
    def load(cls, nc_path: str = None,
             storm: str = _DEFAULT_STORM) -> 'ERA5Steering':
        """Load from a pre-downloaded NetCDF file.

        Parameters
        ----------
        nc_path : explicit path override (default: storm-specific path)
        storm   : 'hugo' or 'katrina' (default: 'hugo')
        """
        cfg = STORM_CONFIGS[storm.lower()]
        if nc_path is None:
            nc_path = cfg["nc_path"]

        if not _HAS_NC4:
            raise ImportError(
                "netCDF4 not found.  Install with:\n"
                "  conda activate woe_env && conda install netcdf4"
            )
        if not os.path.exists(nc_path):
            raise FileNotFoundError(
                f"ERA5 file not found: {nc_path}\n"
                f"Run: python -m oracle_v8.era5_steering --download "
                f"--storm {storm}"
            )

        import netCDF4 as _nc
        with _nc.Dataset(nc_path) as ds:
            raw_time = ds.variables['valid_time'][:]
            from datetime import datetime, timezone
            init_epoch = datetime(
                cfg["init_year"], cfg["init_month"], cfg["init_day"],
                cfg["init_hour"], tzinfo=timezone.utc
            ).timestamp()
            times_h = (np.array(raw_time, dtype=float) - init_epoch) / 3600.0

            lat_raw = np.array(ds.variables['latitude'][:],  dtype=float)
            lon_raw = np.array(ds.variables['longitude'][:], dtype=float)

            lon_raw = np.where(lon_raw < 0, lon_raw + 360.0, lon_raw)
            lon_w   = 360.0 - lon_raw

            if lat_raw[0] > lat_raw[-1]:
                lat_raw  = lat_raw[::-1]
                flip_lat = True
            else:
                flip_lat = False

            sort_idx = np.argsort(lon_w)
            lon_w    = lon_w[sort_idx]

            u_raw = np.array(ds.variables['u'][:], dtype=float)
            v_raw = np.array(ds.variables['v'][:], dtype=float)

            if flip_lat:
                u_raw = u_raw[:, :, ::-1, :]
                v_raw = v_raw[:, :, ::-1, :]
            u_raw = u_raw[:, :, :, sort_idx]
            v_raw = v_raw[:, :, :, sort_idx]

            n_lev = len(_LEVELS_HPA)
            u_dlm = np.zeros((len(times_h), len(lat_raw), len(lon_w)))
            v_dlm = np.zeros_like(u_dlm)
            for k in range(n_lev):
                u_dlm += _DLM_WEIGHTS[k] * u_raw[:, k, :, :]
                v_dlm += _DLM_WEIGHTS[k] * v_raw[:, k, :, :]

        print(
            f"ERA5Steering loaded [{storm.capitalize()}]: "
            f"{len(times_h)} times "
            f"({times_h[0]:+.1f}h to {times_h[-1]:+.1f}h relative to init), "
            f"lat {lat_raw[0]:.1f}–{lat_raw[-1]:.1f}°N, "
            f"lon {lon_w[0]:.1f}–{lon_w[-1]:.1f}°W"
        )

        obj = cls(times_h, lat_raw, lon_w, u_dlm, v_dlm)
        obj._u_levels = u_raw
        obj._v_levels = v_raw
        obj._storm    = storm.lower()
        obj._cfg      = cfg
        obj._level_itps_u = [
            RegularGridInterpolator(
                (times_h, lat_raw, lon_w), u_raw[:, k, :, :],
                method='linear', bounds_error=False, fill_value=None)
            for k in range(n_lev)
        ]
        obj._level_itps_v = [
            RegularGridInterpolator(
                (times_h, lat_raw, lon_w), v_raw[:, k, :, :],
                method='linear', bounds_error=False, fill_value=None)
            for k in range(n_lev)
        ]
        return obj

    # ── interpolation ─────────────────────────────────────────────────────────

    def get_dlm(self, t_hours: float, lat: float, lon_w: float,
                inner_deg: float = 3.0, outer_deg: float = 7.0):
        """
        Return DLM steering (u_env, v_env) as a spatial ring-average.

        Sampling ERA5 AT the TC centre gives Hugo's own circulation winds
        (u ≈ −30 m/s westward at the south side of the vortex), not the
        large-scale environmental steering.  We instead average the DLM
        over an annulus centred on the TC that excludes the inner core.

        The ring average (default 3–7° radius) captures the synoptic flow
        while the cyclonic circulation largely cancels in the average.
        ERA5 resolution ≈ 0.25° so the ring contains ~100–200 sample points.

        Parameters
        ----------
        t_hours   : hours since Hugo init
        lat       : storm latitude (°N)
        lon_w     : storm longitude (°W, positive OR negative — both accepted)
        inner_deg : inner radius of the averaging annulus in degrees (default 3°)
        outer_deg : outer radius of the averaging annulus in degrees (default 7°)

        Returns
        -------
        (u_env, v_env) in m/s  — ERA5 convention (u eastward+, v northward+)
        """
        lon_w = abs(lon_w)   # accept both -74.5 and 74.5 for 74.5°W
        # Build sample grid at 1° spacing across the outer box
        lats_s = np.arange(lat  - outer_deg, lat  + outer_deg + 0.5, 1.0)
        lons_s = np.arange(lon_w - outer_deg, lon_w + outer_deg + 0.5, 1.0)

        pts = []
        for la in lats_s:
            for lo in lons_s:
                r = np.sqrt((la - lat) ** 2 + (lo - lon_w) ** 2)
                if inner_deg <= r <= outer_deg:
                    pts.append([t_hours, la, lo])

        if not pts:
            raise ValueError(
                f"No sample points in ring [{inner_deg}°–{outer_deg}°] "
                f"around ({lat:.2f}°N, {lon_w:.2f}°W)"
            )

        pts    = np.array(pts)
        u_mean = float(self._u_itp(pts).mean())
        v_mean = float(self._v_itp(pts).mean())
        return u_mean, v_mean

    def get_track_mean_dlm(self):
        """
        Compute the time-average DLM over the storm's full observed track.

        For storms like Katrina whose steering evolves dramatically over the run
        (v_env: +1.3 → +6.1 m/s over 35h), the t=0 DLM alone underrepresents
        the steering.  This method averages the DLM over all obs_track positions
        to give a better constant representative for the initial background flow.

        Returns
        -------
        (u_mean, v_mean) : m/s
        """
        u_vals, v_vals = [], []
        for t_h, (lat, lon) in sorted(self._cfg["obs_track"].items()):
            u, v = self.get_dlm(t_h, lat, lon)
            u_vals.append(u)
            v_vals.append(v)
        u_mean = float(np.mean(u_vals))
        v_mean = float(np.mean(v_vals))
        print(f"  ERA5 track-mean DLM: u={u_mean:+.2f} m/s  v={v_mean:+.2f} m/s"
              f"  (avg over {len(u_vals)} obs track points)")
        return u_mean, v_mean

    def get_dlm_per_level(self, t_hours: float, lat: float, lon_w: float,
                          inner_deg: float = 3.0, outer_deg: float = 7.0):
        """
        Diagnostic: DLM ring-average broken down by pressure level.
        Prints a table for verifying the steering at each level.
        """
        lats_s = np.arange(lat  - outer_deg, lat  + outer_deg + 0.5, 1.0)
        lons_s = np.arange(lon_w - outer_deg, lon_w + outer_deg + 0.5, 1.0)
        pts = []
        for la in lats_s:
            for lo in lons_s:
                r = np.sqrt((la - lat) ** 2 + (lo - lon_w) ** 2)
                if inner_deg <= r <= outer_deg:
                    pts.append([t_hours, la, lo])
        pts = np.array(pts)

        print(f"\nPer-level ring-avg winds (r={inner_deg}°–{outer_deg}°) "
              f"at t={t_hours:.1f}h, lat={lat:.2f}°N, lon={lon_w:.2f}°W "
              f"({len(pts)} sample points):")
        print(f"  {'Level':>8}  {'weight':>7}  {'u (m/s)':>9}  {'v (m/s)':>9}")
        print(f"  {'-----':>8}  {'------':>7}  {'-------':>9}  {'-------':>9}")
        u_sum = v_sum = 0.0
        for k, (lev, wt) in enumerate(zip(_LEVELS_HPA, _DLM_WEIGHTS)):
            u_k = float(self._level_itps_u[k](pts).mean())
            v_k = float(self._level_itps_v[k](pts).mean())
            u_sum += wt * u_k
            v_sum += wt * v_k
            print(f"  {lev:>5.0f} hPa  {wt:>7.3f}  {u_k:>+9.2f}  {v_k:>+9.2f}")
        print(f"  {'DLM':>8}  {'':>7}  {u_sum:>+9.2f}  {v_sum:>+9.2f}")
        return u_sum, v_sum

    def print_summary(self, t_start=0.0, t_end=35.0, dt_h=6.0):
        """Print DLM winds at observed track positions for sanity-check."""
        obs_track = self._cfg["obs_track"]
        storm_name = self._storm.capitalize()
        print(f"\nERA5 DLM steering at {storm_name} HURDAT2 best-track positions:")
        print(f"  {'t(h)':>5}  {'lat':>6}  {'lon':>6}  {'u_env':>7}  {'v_env':>7}")
        print(f"  {'----':>5}  {'---':>6}  {'---':>6}  {'------':>7}  {'------':>7}")
        for t, (lat, lon) in sorted(obs_track.items()):
            if t_start <= t <= t_end:
                u, v = self.get_dlm(t, lat, lon)
                print(f"  {t:>5.1f}  {lat:>6.2f}  {lon:>6.2f}  {u:>+7.2f}  {v:>+7.2f}")


# ─── integration helpers for run_hugo.py ─────────────────────────────────────

def update_coriolis_steering(config, u_env: float, v_env: float) -> None:
    """
    Update CoriolisComponent steering winds in-place.

    Tries direct attribute assignment first (works if CoriolisComponent is a
    non-frozen dataclass).  Falls back to object.__setattr__ for frozen cases.

    Call this each output step (or each model step if you want sub-step updates).

    Parameters
    ----------
    config  : OperatorConfig instance (from solver.operator_config)
    u_env   : new zonal steering wind (m/s, ERA5 convention: eastward positive)
    v_env   : new meridional steering wind (m/s, northward positive)
    """
    try:
        config.coriolis.u_env = u_env
        config.coriolis.v_env = v_env
    except (AttributeError, TypeError):
        # Frozen dataclass fallback
        object.__setattr__(config.coriolis, 'u_env', u_env)
        object.__setattr__(config.coriolis, 'v_env', v_env)


# ─── CLI entry point ─────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Oracle V8 ERA5 steering winds utility'
    )
    parser.add_argument(
        '--download', action='store_true',
        help='Download ERA5 data for Hugo 1989 (requires ~/.cdsapirc)'
    )
    parser.add_argument(
        '--storm', default=_DEFAULT_STORM,
        choices=list(STORM_CONFIGS.keys()),
        help=f'Storm to operate on (default: {_DEFAULT_STORM})'
    )
    parser.add_argument(
        '--check', action='store_true',
        help='Print DLM summary at observed track positions'
    )
    parser.add_argument(
        '--levels', action='store_true',
        help='Print per-level ring-avg winds at t=0 init position'
    )
    parser.add_argument(
        '--path', default=None,
        help='Override NetCDF file path'
    )
    args = parser.parse_args()

    if args.download:
        download_era5(output_path=args.path, storm=args.storm)

    if args.check:
        s = ERA5Steering.load(nc_path=args.path, storm=args.storm)
        s.print_summary()

    if args.levels:
        s    = ERA5Steering.load(nc_path=args.path, storm=args.storm)
        cfg  = STORM_CONFIGS[args.storm]
        lat0 = list(cfg["obs_track"].values())[0][0]
        lon0 = list(cfg["obs_track"].values())[0][1]
        s.get_dlm_per_level(0.0, lat0, lon0)

    if not (args.download or args.check or args.levels):
        parser.print_help()
