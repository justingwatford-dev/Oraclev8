"""Shared figure style for the companion papers (2026-07).

Convention (fixed assignment, never cycled — Sol's cross-paper system):
    COMPACT  rust   — compact-taper / control / paper-1 model
    ENVELOPE blue   — Gaussian-envelope / treatment
    OBS      black  — observations / canonical references
    REF      gray   — reference lines, bands, guides
Rust/blue is the canonical CVD-safe complementary axis; guard-violated or
falsified items use OPEN markers, confirmed use FILLED (ledger semantics).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COMPACT  = "#B34A21"   # rust
ENVELOPE = "#2A6FB0"   # blue
OBS      = "#1A1A1A"
REF      = "#9A9A9A"
BAND     = "#2A6FB0"   # canonical-band fill (used at low alpha)

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.22, "grid.linewidth": 0.5,
    "lines.linewidth": 1.6, "figure.dpi": 120, "savefig.dpi": 300,
    "legend.frameon": False, "legend.fontsize": 8.5,
})


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", bbox_inches="tight")
    print(f"  wrote {stem}.png/.pdf")
