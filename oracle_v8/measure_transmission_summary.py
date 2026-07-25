"""
Oracle V8 — transmission-ratio single source of truth
=====================================================
Review response, 2026-07-24. The reviewer asked whether the aggregate and
the mean transmission ratios "are computed in the same place." Honest
answer: they were not — the guard-clean mean lived in the Stage-3 A/B
append, the 0.42 product in the AM-budget accounting, and Fig. 3's fit in
the figure script. They are now all computed HERE, from one input table,
and the manuscripts cite the values this script prints.

Populations, stated once:
  guard-clean six-storm set = {Katrina, Fran, Michael, Laura}
  flagged (intensity guard) = {Hugo, Ivan}
  expansion set (dominant-axis ratios, own registration) =
      {Charley (cross), Florence (along), Ida (cross)} — all guard-clean

Rules the manuscripts follow:
  - The EXPORTED number is the guard-clean mean of ratios.
  - A storm excluded from an estimate never appears in that estimate's
    validation target.
  - Aggregates (sum-obs / sum-pred) are reported for the record, labeled
    by population.

Run:  python -m oracle_v8.measure_transmission_summary
"""
from __future__ import annotations

# observed |shift| and strong-form |projection| (km), landfall-fix
# decomposition; six-storm values from the Stage-3 A/B record, expansion
# values from oracle_v8/Logs/{Charley,Florence,Ida} (2026-07-21 runs).
SIX = {
    #  storm     obs    strong  guard-clean
    "Hugo":     (31.5,   93.0,  False),
    "Katrina":  (46.4,  108.0,  True),
    "Ivan":     (86.1,  139.0,  False),
    "Fran":     (37.6,   81.0,  True),
    "Michael":  (14.5,   73.0,  True),
    "Laura":    (19.9,   74.0,  True),
}
EXPANSION = {
    "Charley":  (13.2,  59.0,  True),   # cross axis
    "Florence": (46.0,  73.0,  True),   # along axis (zonal mover)
    "Ida":      (23.0,  97.0,  True),   # cross axis
}
# per-storm decomposition predictions (measure_transmission_decomp /
# measure_expansion_score)
PRED = {"Hugo": 0.45, "Katrina": 0.52, "Ivan": 0.67, "Fran": 0.43,
        "Michael": 0.45, "Laura": 0.38,
        "Charley": 0.42, "Florence": 0.56, "Ida": 0.31}


def ratios(d):
    return {s: o / p for s, (o, p, _) in d.items()}


def agg(d, keys):
    return sum(d[s][0] for s in keys) / sum(d[s][1] for s in keys)


def main():
    r6 = ratios(SIX)
    clean6 = [s for s, (_, _, c) in SIX.items() if c]
    flagged = [s for s in SIX if s not in clean6]

    print("=" * 74)
    print("TRANSMISSION RATIOS — every variant, one place")
    print("=" * 74)
    for s, (o, p, c) in {**SIX, **EXPANSION}.items():
        tag = "" if c else "  guard-flagged"
        pr = PRED.get(s)
        print(f"  {s:<9} obs {o:6.1f} / strong {p:6.1f} = {o/p:5.2f}"
              f"   decomp-pred {pr:4.2f}{tag}")
    print("-" * 74)
    m_clean = sum(r6[s] for s in clean6) / len(clean6)
    print(f"guard-clean MEAN of ratios (THE exported number) : {m_clean:.3f}")
    print(f"guard-clean aggregate  sum-obs/sum-strong        : {agg(SIX, clean6):.3f}")
    print(f"all-six aggregate (includes flagged {'+'.join(flagged)})   : {agg(SIX, list(SIX)):.3f}")
    xs = [(SIX[s][1], SIX[s][0]) for s in clean6]
    slope = sum(x * y for x, y in xs) / sum(x * x for x, _ in xs)
    print(f"through-origin fit, guard-clean (Fig. 3 line)    : {slope:.3f}")
    pm = sum(PRED[s] for s in clean6) / len(clean6)
    print(f"decomposition mean, guard-clean predictions      : {pm:.3f}")
    print(f"  -> honest validation: {pm:.2f} predicted vs {m_clean:.2f} observed "
          f"(overshoot {pm - m_clean:+.2f}, driven by Michael)")
    r9c = {**{s: r6[s] for s in clean6},
           **{s: o / p for s, (o, p, _) in EXPANSION.items()}}
    print(f"nine-storm guard-clean mean (dominant-axis, incl. expansion): "
          f"{sum(r9c.values()) / len(r9c):.3f}")
    print("NOTE: the 0.64 x 0.65 ~ 0.42 product is the Katrina-history "
          "accounting;\n      its comparator is Katrina's observed 0.43, "
          "not any population mean.")


if __name__ == "__main__":
    main()
