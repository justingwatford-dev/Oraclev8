"""
Oracle V8 — vortex diagnostics (V8.4.1)
=======================================

Built to answer the Run 12 along-track-lag question with *measured* quantities
instead of proxies, following Five's "diagnostics first" discipline and the
ensemble's convergence on the Fiorino & Elsberry (1989) ventilation flow.

The load-bearing routine is `ventilation_flow`. It measures the asymmetric
(translation) flow the vortex actually feels at its center, which discriminates
the two surviving hypotheses for the poleward deficit WITHOUT a new campaign:

    measured v_vent ≈ storm v_N (both ~1.4)   → the flow itself is reduced;
                                                 the vortex is faithfully advected
                                                 by a gyre-modified flow → β-gyre /
                                                 reversed-taper mechanism (Gemini).
    measured v_vent ≈ background v_env (~4.5)  → the flow is NOT reduced but the
        but storm only does ~1.4               vortex still lags → numerical
                                                 advection lag (Copilot / Five).

Everything is computed on the perturbation wind where intensity is concerned, so
"intensity" here is vector |V'| at low level with its radius/height — NOT the
max|u|-component the run log has been printing (Five's catch).

All routines take a duck-typed `state` (needs .u, .v, .theta_prime of shape
(nx, ny, nz)) and use the project backend, so they run unchanged on CuPy in the
repo and on NumPy in the test rig.
"""
from __future__ import annotations

from oracle_v8.backend import xp as np
from oracle_v8.backend import to_numpy


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _grids(nx, ny, dx, dy):
    """Cell-center coordinates (m), shapes (nx,1) and (1,ny)."""
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    return x[:, None], y[None, :]


def _radius(nx, ny, dx, dy, cx_m, cy_m):
    """Radius (m) from center (cx_m, cy_m), shape (nx, ny)."""
    x, y = _grids(nx, ny, dx, dy)
    return np.sqrt((x - cx_m) ** 2 + (y - cy_m) ** 2)


def _level_band(nz, base, z_lo_m, z_hi_m):
    """Boolean (nz,) mask for levels with z in [z_lo, z_hi]. Falls back to all."""
    z = getattr(base, "z", None)
    if z is None:
        return np.ones(nz, dtype=bool)
    z = np.asarray(z)[:nz]
    return (z >= z_lo_m) & (z <= z_hi_m)


# ---------------------------------------------------------------------------
# intensity — the real thing, not max|u|
# ---------------------------------------------------------------------------
def low_level_vmax(state, base, cx_m, cy_m, dx, dy,
                   u_env=0.0, v_env=0.0, z_top_m=3000.0):
    """
    Max *perturbation* vector wind |V'| = |(u-u_env, v-v_env)| in the lowest
    z_top_m, with the radius (km) and height (m) at which it occurs.

    Returns dict: vmax_lowlvl, r_vmax_km, z_vmax_m.
    """
    nx, ny, nz = state.u.shape
    up = state.u - u_env
    vp = state.v - v_env
    spd = np.sqrt(up ** 2 + vp ** 2)                      # (nx,ny,nz)

    band = _level_band(nz, base, 0.0, z_top_m)
    if not bool(np.any(band)):
        band = np.zeros(nz, dtype=bool); band[0] = True
    spd_b = np.where(band[None, None, :], spd, -1.0)

    flat = int(to_numpy(np.argmax(spd_b)))
    iz = flat % nz
    iy = (flat // nz) % ny
    ix = flat // (nz * ny)
    vmax = float(to_numpy(spd_b[ix, iy, iz]))

    r = _radius(nx, ny, dx, dy, cx_m, cy_m)
    r_km = float(to_numpy(r[ix, iy])) / 1e3
    z = getattr(base, "z", None)
    z_m = float(to_numpy(np.asarray(z)[iz])) if z is not None else float(iz)
    return {"vmax_lowlvl": vmax, "r_vmax_km": r_km, "z_vmax_m": z_m}


def vmax_3d(state, base, u_env=0.0, v_env=0.0):
    """Max perturbation |V'| anywhere, with its height (m)."""
    nx, ny, nz = state.u.shape
    up = state.u - u_env
    vp = state.v - v_env
    spd = np.sqrt(up ** 2 + vp ** 2)
    flat = int(to_numpy(np.argmax(spd)))
    iz = flat % nz
    vmax = float(to_numpy(spd.reshape(-1)[flat]))
    z = getattr(base, "z", None)
    z_m = float(to_numpy(np.asarray(z)[iz])) if z is not None else float(iz)
    return {"vmax_3d": vmax, "z_vmax3d_m": z_m}


# ---------------------------------------------------------------------------
# ventilation flow — the discriminator (Fiorino & Elsberry 1989)
# ---------------------------------------------------------------------------
def ventilation_flow(state, cx_m, cy_m, dx, dy, base=None,
                     r_vent_m=150e3, z_lo_m=0.0, z_hi_m=12000.0,
                     u_env=None, v_env=None, n_bins=12):
    """
    The asymmetric (translation) flow at the vortex center: the azimuthal mean
    of the *Cartesian* (u, v) over rings out to r_vent_m, vertically averaged
    over the deep layer [z_lo, z_hi].

    A purely symmetric (tangential) vortex contributes ZERO to the Cartesian
    azimuthal mean, so what survives is the environmental + β-gyre flow steering
    the vortex. Ring-binning (rather than a raw disk mean) makes the symmetric
    cancellation robust to a coarse, off-center grid.

    Returns dict:
      u_vent, v_vent           — steering flow felt at the center (m/s)
      u_gyre, v_gyre           — residual vs the imposed uniform background
                                 (= the β-gyre / asymmetric contribution), only
                                 if u_env/v_env supplied. v_gyre < 0 ⇒ the gyre
                                 is dragging the vortex equatorward (Gemini's
                                 prediction for the mid-run taper interaction).
    """
    nx, ny, nz = state.u.shape
    band = _level_band(nz, base, z_lo_m, z_hi_m) if base is not None \
        else np.ones(nz, dtype=bool)
    if not bool(np.any(band)):
        band = np.ones(nz, dtype=bool)

    # deep-layer mean of Cartesian u, v  (nx, ny)
    u_dl = np.mean(state.u[:, :, band], axis=2)
    v_dl = np.mean(state.v[:, :, band], axis=2)

    r = _radius(nx, ny, dx, dy, cx_m, cy_m)               # (nx, ny)
    # azimuthal mean per ring, then average rings inside r_vent
    edges = np.linspace(0.0, r_vent_m, n_bins + 1)
    u_rings, v_rings = [], []
    for b in range(n_bins):
        m = (r >= edges[b]) & (r < edges[b + 1])
        if bool(to_numpy(np.any(m))):
            u_rings.append(np.mean(u_dl[m]))
            v_rings.append(np.mean(v_dl[m]))
    if not u_rings:                                       # center fell off-grid
        m = r < r_vent_m
        u_vent = float(to_numpy(np.mean(u_dl[m])))
        v_vent = float(to_numpy(np.mean(v_dl[m])))
    else:
        u_vent = float(to_numpy(np.mean(np.stack(u_rings))))
        v_vent = float(to_numpy(np.mean(np.stack(v_rings))))

    out = {"u_vent": u_vent, "v_vent": v_vent}
    if u_env is not None and v_env is not None:
        out["u_gyre"] = u_vent - float(u_env)
        out["v_gyre"] = v_vent - float(v_env)
    return out


# ---------------------------------------------------------------------------
# vorticity center + separation from the θ′ center
# ---------------------------------------------------------------------------
def vorticity_center(state, cx_m, cy_m, dx, dy, base=None,
                     window_m=200e3, z_top_m=3000.0):
    """
    Low-level relative-vorticity maximum within a window about (cx_m, cy_m).
    Returns dict: xv_m, yv_m  (the vorticity center, in m).
    """
    nx, ny, nz = state.u.shape
    band = _level_band(nz, base, 0.0, z_top_m) if base is not None \
        else np.ones(nz, dtype=bool)
    u_ll = np.mean(state.u[:, :, band], axis=2)
    v_ll = np.mean(state.v[:, :, band], axis=2)
    dvdx = (np.roll(v_ll, -1, axis=0) - np.roll(v_ll, 1, axis=0)) / (2 * dx)
    dudy = (np.roll(u_ll, -1, axis=1) - np.roll(u_ll, 1, axis=1)) / (2 * dy)
    zeta = dvdx - dudy                                    # (nx, ny)

    r = _radius(nx, ny, dx, dy, cx_m, cy_m)
    zeta_w = np.where(r < window_m, zeta, -np.inf)
    flat = int(to_numpy(np.argmax(zeta_w)))
    iy = flat % ny
    ix = flat // ny
    x, y = _grids(nx, ny, dx, dy)
    return {"xv_m": float(to_numpy(x[ix, 0])),
            "yv_m": float(to_numpy(y[0, iy]))}


def center_separation_km(theta_xy_m, vort_xy_m):
    """Distance (km) between the θ′ center and the vorticity center."""
    dx_ = theta_xy_m[0] - vort_xy_m[0]
    dy_ = theta_xy_m[1] - vort_xy_m[1]
    return float((dx_ ** 2 + dy_ ** 2) ** 0.5) / 1e3


# ---------------------------------------------------------------------------
# integral budgets
# ---------------------------------------------------------------------------
def perturbation_ke(state, u_env=0.0, v_env=0.0):
    """Mean perturbation kinetic energy density 0.5<(u')²+(v')²> (m²/s²)."""
    up = state.u - u_env
    vp = state.v - v_env
    return {"pert_ke": float(to_numpy(0.5 * np.mean(up ** 2 + vp ** 2)))}


def enstrophy(state, dx, dy, base=None, z_top_m=3000.0):
    """Mean low-level enstrophy 0.5<ζ²> (s⁻²)."""
    nx, ny, nz = state.u.shape
    band = _level_band(nz, base, 0.0, z_top_m) if base is not None \
        else np.ones(nz, dtype=bool)
    u_ll = np.mean(state.u[:, :, band], axis=2)
    v_ll = np.mean(state.v[:, :, band], axis=2)
    dvdx = (np.roll(v_ll, -1, axis=0) - np.roll(v_ll, 1, axis=0)) / (2 * dx)
    dudy = (np.roll(u_ll, -1, axis=1) - np.roll(u_ll, 1, axis=1)) / (2 * dy)
    zeta = dvdx - dudy
    return {"enstrophy": float(to_numpy(0.5 * np.mean(zeta ** 2)))}


# ---------------------------------------------------------------------------
# one-call wrapper + log formatter
# ---------------------------------------------------------------------------
def compute_diagnostics(state, base, cx_m, cy_m, dx, dy,
                        u_env=0.0, v_env=0.0, theta_xy_m=None,
                        r_vent_m=120e3):
    """Run all diagnostics; return a flat dict.

    The ventilation flow and Vmax are anchored on the VORTICITY centre, not the
    passed (θ′) centre — the azimuthal-mean cancellation of the symmetric vortex
    requires the rotational centre, and the two can differ by many cells in a
    distorted vortex (Run 13 showed θ-vs-ζ separations up to ~90 km, which
    contaminated an θ-anchored ventilation read with leaked vortex rotation).
    """
    d = {}
    # vorticity centre first — the rotational centre everything else anchors on
    vc = vorticity_center(state, cx_m, cy_m, dx, dy, base=base)
    d.update(vc)
    cxv, cyv = vc["xv_m"], vc["yv_m"]
    d.update(low_level_vmax(state, base, cxv, cyv, dx, dy, u_env, v_env))
    d.update(vmax_3d(state, base, u_env, v_env))
    d.update(ventilation_flow(state, cxv, cyv, dx, dy, base=base,
                              r_vent_m=r_vent_m, u_env=u_env, v_env=v_env))
    if theta_xy_m is not None:
        d["sep_km"] = center_separation_km(theta_xy_m, (cxv, cyv))
    d.update(perturbation_ke(state, u_env, v_env))
    d.update(enstrophy(state, dx, dy, base=base))
    return d


def format_diag_line(d):
    """Compact one-line log: real Vmax, ventilation/gyre flow, sep, KE."""
    s = (f"Vmax'={d['vmax_lowlvl']:5.1f}@{d['r_vmax_km']:4.0f}km/"
         f"{d['z_vmax_m']/1e3:3.1f}km  "
         f"vent=[{d['u_vent']:+.2f},{d['v_vent']:+.2f}]")
    if "v_gyre" in d:
        s += f" gyre=[{d['u_gyre']:+.2f},{d['v_gyre']:+.2f}]"
    if "sep_km" in d:
        s += f"  sep={d['sep_km']:4.0f}km"
    s += f"  KE={d['pert_ke']:.0f}"
    return s
