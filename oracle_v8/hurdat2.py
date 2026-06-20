"""
Oracle V8 — HURDAT2 loader
==========================
Reads the NHC HURDAT2 Atlantic best-track file and returns storm initialization
+ verification data DIRECTLY from the data file, eliminating the hand-transcription
that produced the Hugo storm-ID (AL13 vs AL11), Katrina 145-vs-100 kt, and Ivan
longitude (85 vs 86) errors.

With one storm-agnostic OperatorConfig, the ONLY per-storm inputs are:
    (hurdat2_id, init_datetime)
Everything physical at t=0 — lat, lon, Vmax, MSLP, Coriolis f, the fallback
steering (centred-difference best-track motion), the verification track, and the
landfall record — is read from here.  Domain size and run length are computed
downstream from the init→landfall latitude span, not hand-set.

File:  hurdat2-1851-2025-*.txt  (2024+ edition: 21 columns incl. trailing Rmax;
older editions have 20 and Rmax reads as absent — handled).

Format (comma-separated):
  header:  AL092004,            IVAN,     61,
  fix:     20040914, 1200,  , HU, 22.9N,  85.0W, 125,  928, <12 wind-radii>, <Rmax>
           [0]date   [1]time [2]id [3]status [4]lat [5]lon [6]Vmax(kt) [7]MSLP(mb)
           ... [8..19] 34/50/64-kt wind radii ... [20] Rmax(nm)
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
import math

OMEGA = 7.2921e-5
R_EARTH = 6_371_000.0


def kt_to_ms(kt: float) -> float:
    return kt * 0.51444


def nm_to_m(nm: float) -> float:
    return nm * 1852.0


def coriolis(lat_deg: float) -> float:
    return 2.0 * OMEGA * math.sin(math.radians(lat_deg))


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Fix:
    dt: datetime
    lat: float                 # °N
    lon: float                 # ° signed, W negative
    vmax_kt: float
    mslp_mb: float | None
    status: str                # HU / TS / TD / EX / ...
    record_id: str             # '', 'L' (landfall), 'I', etc.
    rmax_nm: float | None


@dataclass
class StormTrack:
    storm_id: str
    name: str
    year: int
    fixes: list[Fix]

    def fix_at(self, dt: datetime) -> Fix:
        """Exact fix at dt, else linear interpolation between bracketing fixes."""
        for f in self.fixes:
            if f.dt == dt:
                return f
        before = [f for f in self.fixes if f.dt <= dt]
        after = [f for f in self.fixes if f.dt >= dt]
        if not before or not after:
            raise ValueError(f"{dt:%Y-%m-%d %HZ} is outside the {self.storm_id} track range "
                             f"({self.fixes[0].dt:%Y-%m-%d %HZ} – {self.fixes[-1].dt:%Y-%m-%d %HZ})")
        a, b = before[-1], after[0]
        if a.dt == b.dt:
            return a
        w = (dt - a.dt).total_seconds() / (b.dt - a.dt).total_seconds()
        return Fix(dt, a.lat + w * (b.lat - a.lat), a.lon + w * (b.lon - a.lon),
                   a.vmax_kt + w * (b.vmax_kt - a.vmax_kt), None, a.status, "", None)

    def landfall_fix(self, on_or_after: datetime) -> Fix | None:
        """First US-landfall ('L') record at or after the given time."""
        for f in self.fixes:
            if f.dt >= on_or_after and f.record_id == "L":
                return f
        return None


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def _parse_latlon(tok: str) -> float:
    tok = tok.strip()
    val = float(tok[:-1])
    return -val if tok[-1] in ("S", "W") else val


def _missing(tok: str) -> float | None:
    try:
        v = float(tok)
    except (TypeError, ValueError):
        return None
    return None if v == -999 else v


def parse_hurdat2(path: str | Path) -> dict[str, StormTrack]:
    """Parse the whole file into {storm_id: StormTrack}."""
    with open(path, "r", encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    storms: dict[str, StormTrack] = {}
    i = 0
    while i < len(lines):
        hdr = [t.strip() for t in lines[i].split(",")]
        storm_id, name, n = hdr[0], hdr[1], int(hdr[2])
        year = int(storm_id[4:8])
        fixes = []
        j = i + 1
        end = i + 1 + n
        while j < end and j < len(lines):
            c = [t.strip() for t in lines[j].split(",")]
            dt = datetime.strptime(c[0] + c[1], "%Y%m%d%H%M")
            fixes.append(Fix(
                dt=dt,
                lat=_parse_latlon(c[4]),
                lon=_parse_latlon(c[5]),
                vmax_kt=float(c[6]),
                mslp_mb=_missing(c[7]),
                status=c[3],
                record_id=c[2],
                rmax_nm=_missing(c[20]) if len(c) > 20 else None,
            ))
            j += 1
        storms[storm_id] = StormTrack(storm_id, name, year, fixes)
        i = end
    return storms


# ---------------------------------------------------------------------------
# Per-storm init builder  (the storm-agnostic discipline lives here)
# ---------------------------------------------------------------------------

def storm_init(track: StormTrack, init_date: datetime, *,
               B: float = 1.5, rmax_run_m: float = 75_000.0,
               threshold_lat: float | None = None) -> dict:
    """Build the init/verification dict for one storm, all values read from the
    best track at `init_date`.  Mirrors the storm_data.py dict schema so the
    existing run pipeline and landfall_verify consume it unchanged."""
    f0 = track.fix_at(init_date)

    # Fallback steering (used ONLY if ERA5 is unavailable; REQUIRE_ERA5 normally
    # blocks that path).  Centred-difference best-track motion about t=0.
    try:
        pm = track.fix_at(init_date - timedelta(hours=6))
        pp = track.fix_at(init_date + timedelta(hours=6))
        dlat_m = math.radians(pp.lat - pm.lat) * R_EARTH
        dlon_m = math.radians(pp.lon - pm.lon) * R_EARTH * math.cos(math.radians(f0.lat))
        dt_s = 12 * 3600.0
        u_env, v_env = dlon_m / dt_s, dlat_m / dt_s
    except ValueError:
        u_env = v_env = 0.0

    # Verification track: every fix from init forward, as (t_h, lat, lon_signed).
    obs_track = [((fx.dt - init_date).total_seconds() / 3600.0, fx.lat, fx.lon)
                 for fx in track.fixes if fx.dt >= init_date]

    lf = track.landfall_fix(init_date)
    if threshold_lat is None and lf is not None:
        threshold_lat = round(lf.lat, 1)

    rmax_nm = f0.rmax_nm if f0.rmax_nm else rmax_run_m / 1852.0
    return dict(
        name=track.name.title(), year=track.year, basin=track.storm_id[:2],
        hurdat2_id=track.storm_id, init_date=init_date.strftime("%Y-%m-%d %HZ"),
        lat0_deg=f0.lat, lon0_deg=f0.lon,
        Vmax_kt=f0.vmax_kt, Vmax_ms=kt_to_ms(f0.vmax_kt), P_min_mb=f0.mslp_mb,
        Rmax_nm=rmax_nm, Rmax_m=nm_to_m(rmax_nm), B=B, f=coriolis(f0.lat),
        u_env_ms=u_env, v_env_ms=v_env,
        landfall_lat=(lf.lat if lf else None),
        landfall_lon=(lf.lon if lf else None),
        landfall_time_h=((lf.dt - init_date).total_seconds() / 3600.0 if lf else None),
        threshold_lat=threshold_lat,
        obs_track=obs_track,
        init_source="hurdat2_file",
    )


# ---------------------------------------------------------------------------
# Registry — the ONLY hand-set per-storm inputs: id + init start time
# (+ optional threshold override; defaults to the HURDAT2 landfall latitude)
# ---------------------------------------------------------------------------

ATLANTIC_FILE = Path(__file__).with_name("hurdat2.txt")

REGISTRY = {
    "Hugo":    dict(id="AL111989", init="1989-09-21 00Z", threshold_lat=32.5),
    "Katrina": dict(id="AL122005", init="2005-08-28 00Z", threshold_lat=29.1),
    "Ivan":    dict(id="AL092004", init="2004-09-14 12Z", threshold_lat=30.0),
    # add storms here: "Michael": dict(id="AL142018", init="2018-10-09 12Z"),
}


def _parse_init(s: str) -> datetime:
    return datetime.strptime(s.replace("Z", "").strip(), "%Y-%m-%d %H")


def load_storm(name_or_id: str, hurdat2_path: str = ATLANTIC_FILE,
               init_date: str | None = None, **kw) -> dict:
    """Resolve a storm by registry name (e.g. 'Ivan') or raw id (e.g. 'AL092004')
    and return its init/verification dict, read straight from HURDAT2.

      load_storm("Ivan")                         # registry
      load_storm("AL142018", init_date="2018-10-09 12Z", threshold_lat=30.4)
    """
    storms = parse_hurdat2(hurdat2_path)
    if name_or_id in REGISTRY:
        entry = REGISTRY[name_or_id]
        sid = entry["id"]
        init = _parse_init(entry["init"])
        kw.setdefault("threshold_lat", entry.get("threshold_lat"))
    else:
        sid = name_or_id
        if init_date is None:
            raise ValueError(f"{name_or_id} not in REGISTRY — pass init_date='YYYY-MM-DD HHZ'")
        init = _parse_init(init_date)
    if sid not in storms:
        raise KeyError(f"storm id {sid} not found in {hurdat2_path}")
    return storm_init(storms[sid], init, **kw)


# ---------------------------------------------------------------------------
# Validation entry point: load a registry storm and print its init + track so
# you can diff against the hand-coded values in storm_data.py.
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    path = sys.argv[2] if len(sys.argv) > 2 else ATLANTIC_FILE
    name = sys.argv[1] if len(sys.argv) > 1 else "Ivan"
    s = load_storm(name, hurdat2_path=path)
    print(f"{s['name']} ({s['year']})  {s['hurdat2_id']}  init {s['init_date']}")
    print(f"  t=0:  {s['lat0_deg']:.1f}N {abs(s['lon0_deg']):.1f}W  "
          f"{s['Vmax_kt']:.0f}kt  {s['P_min_mb']}mb  f={s['f']:.2e}")
    print(f"  fallback steering (best-track motion): u={s['u_env_ms']:+.1f} v={s['v_env_ms']:+.1f} m/s")
    print(f"  landfall: {s['landfall_lat']}N {abs(s['landfall_lon'])}W  "
          f"t+{s['landfall_time_h']:.2f}h   threshold_lat={s['threshold_lat']}")
    print("  obs_track (diff this against storm_data.py):")
    for t_h, lat, lon in s["obs_track"][:10]:
        print(f"    ({t_h:6.2f}, {lat:6.1f}, {lon:7.1f}),")
