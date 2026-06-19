"""Re-score a run log with landfall_verify against the HURDAT2 obs_track.

Usage:
    python test_footer.py <run_log.txt> <STORM>
    python -m oracle_v8.test_footer <run_log.txt> <STORM>

STORM is one of HUGO, KATRINA, or IVAN. The script parses the printed
(t, lat, lon) track table from an old run log and feeds it through the
storm's HURDAT2 obs_track.
"""
from __future__ import annotations

import re
import sys

try:
    from oracle_v8.landfall_verify import landfall_report
    from oracle_v8 import storm_data
except ImportError:  # Allow running from inside oracle_v8/.
    from landfall_verify import landfall_report
    import storm_data


ROW = re.compile(r"^\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s+[-+]?\d")


def _parse_track(logfile: str):
    track_t, track_lat, track_lon = [], [], []
    in_track = False

    with open(logfile, encoding="utf-8", errors="replace") as f:
        for line in f:
            if "Time(h)" in line:
                in_track = True
                continue
            if not in_track:
                continue
            m = ROW.match(line)
            if not m:
                continue
            track_t.append(float(m.group(1)))
            track_lat.append(float(m.group(2)))
            track_lon.append(-float(m.group(3)))  # logs print abs lon; W negative

    return track_t, track_lat, track_lon


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("Usage: python test_footer.py <run_log.txt> <STORM>")
        print("       STORM: HUGO | KATRINA | IVAN")
        return 2

    logfile, storm_name = argv[0], argv[1].upper()
    if not hasattr(storm_data, storm_name):
        print(f"Unknown storm {storm_name!r}; expected HUGO, KATRINA, or IVAN")
        return 2

    track_t, track_lat, track_lon = _parse_track(logfile)
    if not track_t:
        print(f"No track rows parsed from {logfile!r}")
        return 1

    print(
        f"Parsed {len(track_t)} track points from {logfile} "
        f"(t={track_t[0]}..{track_t[-1]}h)"
    )
    landfall_report(track_t, track_lat, track_lon, getattr(storm_data, storm_name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
