"""
Oracle V8 — Landfall Verification (shared)
==========================================
Along/cross-track landfall verification used by run_hugo / run_katrina /
run_ivan.  Replaces the per-script footer that compared the model's
threshold-latitude crossing against the LANDFALL POINT — a comparison that
folds remaining observed along-track travel into "cross-track bias"
(Hugo re-validation 2 read +137 km east when the same-latitude error was
+91 km vs the pipeline track).

Three layers:
  1. Same-latitude threshold comparison (the headline): the model's and the
     observed track's crossings of `lat_thresh`, separating timing
     (along-track) from cross-track error.
  2. Error decomposition at each observed fix: total / along-track
     (+ = ahead of obs along its motion) / cross-track (+ = right of
     observed motion; ≈ east for a NNW-moving storm).  This is the Run-12
     paper-figure machinery.
  3. The legacy landfall-point metric, labelled, for continuity with the
     Runs 12–15 logs.

The observed track comes from the storm dict's `obs_track` key
(HURDAT2 fixes — see storm_data.py for sourcing) as
(t_h, lat_degN, lon_deg_signed) tuples, west longitude NEGATIVE.
"""
from __future__ import annotations

import numpy as np

KM_PER_DEG_LAT = 111.32


def _km_offset(lat, lon, lat_ref, lon_ref):
    """(east_km, north_km) of (lat, lon) relative to (lat_ref, lon_ref)."""
    east  = (lon - lon_ref) * KM_PER_DEG_LAT * np.cos(np.radians(lat_ref))
    north = (lat - lat_ref) * KM_PER_DEG_LAT
    return east, north


def _interp_crossing(ts, lats, lons, lat_thresh):
    """First northward crossing of lat_thresh → (t, lon) or (None, None)."""
    for i in range(1, len(lats)):
        if lats[i] >= lat_thresh and lats[i - 1] < lat_thresh:
            frac = (lat_thresh - lats[i - 1]) / (lats[i] - lats[i - 1])
            return (ts[i - 1]  + frac * (ts[i]  - ts[i - 1]),
                    lons[i - 1] + frac * (lons[i] - lons[i - 1]))
    return None, None


def landfall_report(track_t, track_lat, track_lon, s,
                    obs_track=None, lat_thresh=None):
    """Print the along/cross-track landfall verification.

    Parameters
    ----------
    track_t, track_lat, track_lon : model track (lon signed, west negative)
    s          : storm dict (storm_data.py)
    obs_track  : [(t_h, lat, lon_signed), ...]; defaults to s["obs_track"]
    lat_thresh : crossing latitude; defaults to s["threshold_lat"],
                 then s["landfall_lat"]
    """
    if obs_track is None:
        obs_track = s.get("obs_track")
    if lat_thresh is None:
        lat_thresh = s.get("threshold_lat", s["landfall_lat"])

    print(f"\n{'=' * 68}")
    print(f"LANDFALL VERIFICATION - threshold {lat_thresh} degN")
    print("=" * 68)

    # -- model threshold crossing (needed by every layer) --------------------
    m_t, m_lon = _interp_crossing(track_t, track_lat, track_lon, lat_thresh)
    if m_t is None:
        print(f"V8: storm did not cross {lat_thresh} degN in run window")
        print("=" * 68)
        return

    if not obs_track:
        print("  (no obs_track in storm dict - legacy metric only)")
    else:
        obs_t   = [p[0] for p in obs_track]
        obs_lat = [p[1] for p in obs_track]
        obs_lon = [p[2] for p in obs_track]

        # -- 1. Same-latitude threshold crossing -----------------------------
        o_t, o_lon = _interp_crossing(obs_t, obs_lat, obs_lon, lat_thresh)
        if o_t is None:
            print(f"  (obs track never crosses {lat_thresh} degN - check obs_track)")
        else:
            dt_h  = m_t - o_t
            dx_km = ((m_lon - o_lon)
                     * KM_PER_DEG_LAT * np.cos(np.radians(lat_thresh)))
            print(f"  Obs crossed:  {lat_thresh} degN at {abs(o_lon):.2f} degW  "
                  f"t+{o_t:.1f}h")
            print(f"  V8  crossed:  {lat_thresh} degN at {abs(m_lon):.2f} degW  "
                  f"t+{m_t:.1f}h")
            print(f"  Timing:       {dt_h:+.1f} h  "
                  f"({'early' if dt_h < 0 else 'late'})")
            print(f"  Cross-track:  {dx_km:+.1f} km  "
                  f"({'east' if dx_km > 0 else 'west'})  [same-latitude]")

        # -- 2. Along/cross decomposition at observed fixes -------------------
        print(f"\n  Track-error decomposition at observed fixes")
        print(f"  (along + = ahead of obs;  "
              f"cross + = right of obs motion, approx E here)")
        print(f"  {'t(h)':>5}  {'obs (N/W)':>13}  {'V8 (N/W)':>13}  "
              f"{'total':>7}  {'along':>7}  {'cross':>7}")
        t_arr = np.asarray(track_t)
        for i, (t_o, la_o, lo_o) in enumerate(obs_track):
            if t_o < t_arr[0] or t_o > t_arr[-1]:
                continue   # outside the model run window
            la_m = float(np.interp(t_o, track_t, track_lat))
            lo_m = float(np.interp(t_o, track_t, track_lon))
            e_e, e_n = _km_offset(la_m, lo_m, la_o, lo_o)
            total = float(np.hypot(e_e, e_n))
            # Observed motion direction — centred difference where possible.
            j0, j1 = max(i - 1, 0), min(i + 1, len(obs_track) - 1)
            m_e, m_n = _km_offset(obs_lat[j1], obs_lon[j1],
                                  obs_lat[j0], obs_lon[j0])
            norm = float(np.hypot(m_e, m_n))
            if norm == 0.0:
                along = cross = float("nan")
            else:
                ue, un = m_e / norm, m_n / norm     # along-track unit vector
                along  = e_e * ue + e_n * un
                cross  = e_e * un - e_n * ue        # right-of-motion (un, -ue)
            print(f"  {t_o:>5.1f}  {la_o:>6.2f}/{abs(lo_o):>6.2f}  "
                  f"{la_m:>6.2f}/{abs(lo_m):>6.2f}  "
                  f"{total:>7.1f}  {along:>+7.1f}  {cross:>+7.1f}")

    # -- 3. Legacy landfall-point metric (Runs 12-15 continuity) -------------
    err_legacy = ((m_lon - s["landfall_lon"])
                  * KM_PER_DEG_LAT * np.cos(np.radians(lat_thresh)))
    print(f"\n  Legacy metric - vs landfall point {s['landfall_lat']} degN "
          f"{abs(s['landfall_lon']):.1f} degW t+{s['landfall_time_h']:.0f}h "
          f"(conflates along+cross):")
    print(f"    {err_legacy:+.1f} km  ({'east' if err_legacy > 0 else 'west'})")
    print("=" * 68)
