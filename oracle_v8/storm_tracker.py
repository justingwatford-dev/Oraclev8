"""
Oracle V8 Storm Tracker
=======================
Lightweight replacement for find_storm_centre() in run_hugo.py.

Three-layer robustness (lessons from V7 StormTracker V50.4):
  Layer 1 — Interior masking:    exclude domain edges (12% margin) from
                                  the θ′ pressure proxy search, exactly like
                                  V7's sponge-region pressure mask.
  Layer 2 — Continuity gate:     reject positions more than max_jump_cells
                                  from the previous fix; fall through to
                                  vorticity fallback.
  Layer 3 — Vorticity CoM:       positive-vorticity centre-of-mass over the
                                  lower half of the domain (z-levels 0..nz//2),
                                  analogous to V7's Vorticity Anchor used in
                                  GENESIS and HURRICANE-VORT phases.

Primary signal: θ′ hydrostatic pressure proxy (V8 barotropic config).
  p_sfc = -∫ ρ₀(z) · g·θ′/θ₀(z) dz      (column integral, each (x,y))
  anomaly = p_sfc − ⟨p_sfc⟩

Fallback signal: vertical vorticity CoM.
  ζ = ∂v/∂x − ∂u/∂y  (finite differences)
  weighted by positive ζ, column-averaged over lower nz//2 levels.

Usage
-----
  from oracle_v8.storm_tracker import find_storm_centre

  # Inside run_hugo.py main loop:
  last_pos = None
  for step in range(n_steps):
      ...
      (x_c, y_c), method = find_storm_centre(
          state, base_w, nx, ny, dx,
          last_pos=last_pos,
      )
      last_pos = (x_c, y_c)
      # AXIS: x_c = E–W (→ lon), y_c = N–S (→ lat).  Do NOT swap these.
      lat_c = init_lat + (y_c - Ly/2) / 111_000
      lon_c = init_lon - (x_c - Lx/2) / (111_000 * cos(radians(lat_c)))

Returns
-------
  (x_m, y_m) : storm centre in metres from SW domain corner
  method      : 'theta_prime' | 'vorticity' | 'extrapolated'
                (logged every output step for diagnostics)
"""

import numpy as np
from oracle_v8.backend import xp as _xp, to_numpy

# ─── constants ────────────────────────────────────────────────────────────────
_G = 9.81          # m s⁻²

# Gaussian localisation radius (grid cells), shared by the θ′ proxy weight
# (Layer 1) and the vorticity-CoM fallback (Layer 3) so the two windows can
# never drift apart.  10 cells ≈ 156 km at dx = 15.6 km.
_R0_CELLS = 10.0


# ─── helpers ──────────────────────────────────────────────────────────────────

def _theta_prime_proxy(state, base_w, xp):
    """
    Hydrostatic surface-pressure anomaly from θ′ field.
    Returns 2D numpy array (nx, ny).
    """
    b         = _G * state.theta_prime / base_w.theta0[None, None, :]  # buoyancy
    integrand = base_w.rho0[None, None, :] * b                          # ρ₀·b
    p_sfc     = -xp.trapz(integrand, base_w.z, axis=2)                 # (nx, ny)
    p_anm     = p_sfc - xp.mean(p_sfc)                                 # anomaly
    return to_numpy(p_anm)                                              # to CPU


def _gaussian_weight(p_anm, nx, ny, dx, last_pos, r0_cells=10.0):
    """
    Apply Gaussian distance weight to the θ′ pressure proxy field.

    Multiplies p_anm (nx × ny, negative minimum at TC centre) by
    exp(−r² / 2r₀²) centred on last_pos.  False minima from secondary
    circulation features and beta-gyres at large radii are suppressed;
    the primary eyewall signal near last_pos is preserved.

    r0_cells=10 (≈156 km for dx=15.6 km) gives:
      r=4.8 cells (Rmax=75 km):  weight=0.89   ← eyewall preserved
      r=10  cells (156 km):      weight=0.61
      r=20  cells (312 km):      weight=0.14
      r=32  cells (500 km):      weight=0.006  ← beta-gyres suppressed
    """
    if last_pos is None:
        return p_anm
    ix_last = last_pos[0] / dx
    iy_last = last_pos[1] / dx
    ix_arr  = np.arange(nx, dtype=float).reshape(-1, 1)   # (nx, 1)
    iy_arr  = np.arange(ny, dtype=float).reshape(1, -1)   # (1, ny)
    r2      = (ix_arr - ix_last) ** 2 + (iy_arr - iy_last) ** 2
    return p_anm * np.exp(-r2 / (2.0 * r0_cells ** 2))



def _smooth3(f, xp):
    """Separable 3-point [1,2,1]/4 smooth in x then y (periodic via roll).

    Suppresses single-cell vorticity spikes so the arg-max locks onto the
    broad eye-scale peak rather than a noisy secondary filament.  Applied to
    the Gaussian-windowed field, the periodic wrap at the domain edge is
    already negligible (window weight ≈ 0 there).
    """
    f = (xp.roll(f, 1, axis=0) + 2.0 * f + xp.roll(f, -1, axis=0)) * 0.25
    f = (xp.roll(f, 1, axis=1) + 2.0 * f + xp.roll(f, -1, axis=1)) * 0.25
    return f


def _vorticity_centre(u, v, nx, ny, dx, xp, last_pos=None, r0_cells=_R0_CELLS):
    """
    Locate the storm centre from the lower-tropospheric vorticity field
    (ζ = ∂v/∂x − ∂u/∂y, positive part, averaged over the lower nz//2 levels).

    Two estimators:

    * `last_pos` given  → ARG-MAX of the Gaussian-windowed, lightly-smoothed
      ζ, with a sub-cell centre-of-mass refinement over a 5×5 patch about the
      peak.  The eye is the strongest, most compact cyclonic vorticity peak;
      an arg-max is NOT pulled off it by surrounding asymmetric vorticity the
      way a centre of mass is.  This replaces the V8.3.3 windowed-CoM, which
      lagged the eye south-west in Katrina Run 9 (the 13 h "park" at ~27°N,
      t=24–37h): surface drag adds rear-flank boundary-layer vorticity that
      biased the centroid toward the SW quadrant.  The peak is immune to that.

    * `last_pos=None`   → original GLOBAL centre of mass (Hugo / legacy path,
      unchanged).

    The Gaussian window (same r₀ as the θ′ proxy) still confines the search to
    the eye's neighbourhood so a distant secondary peak can't capture the
    arg-max.

    Returns (x_m, y_m) in metres from SW corner.
    """
    nz_half = u.shape[2] // 2
    u_low   = u[:, :, :nz_half]
    v_low   = v[:, :, :nz_half]

    dvdx = (xp.roll(v_low, -1, axis=0) - xp.roll(v_low, 1, axis=0)) / (2.0 * dx)
    dudy = (xp.roll(u_low, -1, axis=1) - xp.roll(u_low, 1, axis=1)) / (2.0 * dx)
    zeta     = dvdx - dudy                           # (nx, ny, nz_half)
    zeta_pos = xp.maximum(zeta, 0.0)                 # cyclonic only
    zeta_col = xp.mean(zeta_pos, axis=2)             # (nx, ny) column average

    ix_arr = xp.arange(nx, dtype=xp.float64)[:, None]   # (nx,1)
    iy_arr = xp.arange(ny, dtype=xp.float64)[None, :]   # (1,ny)

    # ── Legacy global CoM (no prior fix) ────────────────────────────────────
    if last_pos is None:
        total = float(xp.sum(zeta_col))
        if total < 1e-14:
            return (nx * dx / 2.0), (ny * dx / 2.0)   # ultimate fallback
        xc = float(xp.sum(ix_arr * zeta_col) / total) * dx + 0.5 * dx
        yc = float(xp.sum(iy_arr * zeta_col) / total) * dx + 0.5 * dx
        return xc, yc

    # ── Windowed arg-max + sub-cell refinement (production fallback) ─────────
    r2       = (ix_arr - last_pos[0] / dx) ** 2 + (iy_arr - last_pos[1] / dx) ** 2
    zeta_win = _smooth3(zeta_col, xp) * xp.exp(-r2 / (2.0 * r0_cells ** 2))

    if float(xp.max(zeta_win)) < 1e-14:
        # No usable cyclonic vorticity in the window — hold last position
        return last_pos[0], last_pos[1]

    # Peak cell of the windowed, smoothed field
    idx        = int(xp.argmax(zeta_win))
    ix0, iy0   = divmod(idx, ny)

    # Sub-cell centre of mass over a 5×5 patch about the peak (clipped to grid)
    w        = 2
    i0, i1   = max(0, ix0 - w), min(nx, ix0 + w + 1)
    j0, j1   = max(0, iy0 - w), min(ny, iy0 + w + 1)
    patch    = xp.maximum(zeta_win[i0:i1, j0:j1], 0.0)
    ptot     = float(xp.sum(patch))
    if ptot < 1e-14:
        xc_cell, yc_cell = float(ix0), float(iy0)
    else:
        ii = xp.arange(i0, i1, dtype=xp.float64)[:, None]
        jj = xp.arange(j0, j1, dtype=xp.float64)[None, :]
        xc_cell = float(xp.sum(ii * patch) / ptot)
        yc_cell = float(xp.sum(jj * patch) / ptot)

    return (xc_cell + 0.5) * dx, (yc_cell + 0.5) * dx


def _masked_argmin(field_2d, nx, ny, margin):
    """
    Find argmin of field_2d, excluding the outer `margin` cells on all sides.
    Returns (ix, iy) in grid cells.
    """
    p = field_2d.copy()
    p[:margin,  :]  = p.max()   # edges → maximum so argmin ignores them
    p[-margin:, :]  = p.max()
    p[:,  :margin]  = p.max()
    p[:, -margin:]  = p.max()
    idx    = int(np.argmin(p))
    ix, iy = divmod(idx, ny)
    return ix, iy


def _parabolic_refine(p, ix, iy):
    """Sub-cell refinement of an argmin cell via 1-D parabola fits (V8.6.2).

    The Galilean-control trace showed cell-snapped fixes flicker with ±1.5
    cell amplitude, putting a 23 km noise floor under every track at
    dx=15.6 km.  For a smooth extremum the 3-point parabola vertex recovers
    the sub-cell position to ~0.1 cell.  Offsets are clamped to ±0.5 cell;
    edge cells return unrefined.  Works for minima (θ′ proxy) directly.
    """
    nx_, ny_ = p.shape

    def _off(a, b, c):
        d = a - 2.0 * b + c
        if abs(d) < 1e-30:
            return 0.0
        return float(np.clip(0.5 * (a - c) / d, -0.5, 0.5))

    ox = _off(p[ix - 1, iy], p[ix, iy], p[ix + 1, iy]) \
        if 0 < ix < nx_ - 1 else 0.0
    oy = _off(p[ix, iy - 1], p[ix, iy], p[ix, iy + 1]) \
        if 0 < iy < ny_ - 1 else 0.0
    return ix + ox, iy + oy


# ─── public API ───────────────────────────────────────────────────────────────

def find_storm_centre(
    state,
    base_w,
    nx,
    ny,
    dx,
    last_pos        = None,
    max_jump_cells  = 5,
    interior_margin = 0.12,
    v_env_ms        = 0.0,
    u_env_ms        = 0.0,
    dt_output_s     = 1800.0,
    return_diag     = False,
    last_raw_pos    = None,
    reacq_jump_cells = None,
    subcell         = True,
):
    """
    Locate storm centre with three-layer robustness.

    Parameters
    ----------
    state           : model State with fields (u, v, theta_prime, …)
    base_w          : base-state object with (rho0, theta0, z) as 1-D arrays
    nx, ny, dx      : grid size and spacing (m)
    last_pos        : (x_m, y_m) from previous call, or None on first step
    max_jump_cells  : Layer-2 continuity gate — positions further than this
                      from last_pos trigger vorticity fallback (default 5 cells
                      ≈ 78 km for dx=15.6 km; at 0.5 h output interval this
                      allows up to 156 km/h apparent motion before rejecting)
    interior_margin : Layer-1 edge fraction to exclude (default 0.12 = 12%;
                      same as V7 V50.4 sponge_margin)
    v_env_ms        : northward background flow (m/s) — used for dead reckoning
                      when both θ′ and vorticity signals fail simultaneously
                      (secondary intensification case).  Pass the current
                      CoriolisComponent v_env.  Default 0 = hold last position.
    u_env_ms        : eastward background flow (m/s, positive=east) for dead
                      reckoning.  x-axis = northward in Oracle domain.
    dt_output_s     : diagnostic output interval in seconds (default 1800 = 0.5h)
    return_diag     : if True, also return a diagnostics dict containing
                      'theta_raw_m' — the raw (unwindowed) θ′ centre in metres.
    last_raw_pos    : (x_m, y_m) raw θ′ centre from the PREVIOUS call, or None.
                      Used by Layer 2a re-acquisition to decide whether the eye
                      is moving coherently.  Pass the previous diag's
                      'theta_raw_m'.
    reacq_jump_cells: max raw-θ′ step (cells) still considered a coherent eye
                      track for re-acquisition (default = max_jump_cells).

    Returns
    -------
    (x_m, y_m) : storm centre in metres from SW domain corner
    method     : str — 'theta_prime' | 'theta_reacq' | 'vorticity' | 'extrapolated'
    diag       : dict (only if return_diag=True) — {'theta_raw_m': (x_m, y_m)}
    """
    xp     = _xp
    margin = max(int(interior_margin * nx), 3)

    # ── Layer 1: θ′ proxy with interior masking + Gaussian weight ────────────
    p_anm = _theta_prime_proxy(state, base_w, xp)   # (nx, ny) numpy

    # RAW θ′ centre — same masked argmin, but WITHOUT the Gaussian continuity
    # window.  Diagnostic only: never steers the fix.  Logging it alongside the
    # chosen method shows whether the eye is marching while a windowed estimator
    # lags/parks (Run 9 [vor] park).  _masked_argmin copies internally, so this
    # does not mutate p_anm before the weight is applied below.
    ix_raw, iy_raw = _masked_argmin(p_anm, nx, ny, margin)
    if subcell:
        fx_raw, fy_raw = _parabolic_refine(p_anm, ix_raw, iy_raw)
    else:
        fx_raw, fy_raw = float(ix_raw), float(iy_raw)
    theta_raw_m    = ((fx_raw + 0.5) * dx, (fy_raw + 0.5) * dx)

    # Gaussian distance weight centred on last_pos — suppresses false minima
    # from secondary circulation features and beta-gyres (Gemini, brief 4).
    p_anm = _gaussian_weight(p_anm, nx, ny, dx, last_pos, r0_cells=_R0_CELLS)
    ix, iy = _masked_argmin(p_anm, nx, ny, margin)
    if subcell:
        fx, fy = _parabolic_refine(p_anm, ix, iy)   # sub-cell (V8.6.2);
    else:                                            # continuity gate below
        fx, fy = float(ix), float(iy)                # still uses cell ix/iy
    x_new  = (fx + 0.5) * dx
    y_new  = (fy + 0.5) * dx
    method = 'theta_prime'

    # ── Layer 2: continuity gate ──────────────────────────────────────────────
    if last_pos is not None:
        last_ix = last_pos[0] / dx
        last_iy = last_pos[1] / dx
        dist    = np.sqrt((ix - last_ix) ** 2 + (iy - last_iy) ** 2)

        if dist > max_jump_cells:
            # The continuity gate has rejected the windowed θ′ fix as too large
            # a jump from last_pos.  But last_pos can itself be STALE: once the
            # [vor] fallback plants it below the true eye, the gate then rejects
            # the (correct) θ′ fix every step because it keeps "jumping" back up
            # to the eye — a self-locking trap that parked Katrina Run 10 ~0.6°
            # south of the eye for t=24–42h while the raw θ′ minimum (the eye)
            # marched north the whole time.
            #
            # Layer 2a — RE-ACQUISITION: before falling to the fallback, ask the
            # raw (unwindowed) θ′ minimum whether it is a COHERENT, smoothly
            # moving feature — i.e. it stepped < reacq_jump_cells from its own
            # previous position.  The eye is the globally deepest pressure
            # minimum, so a coherent raw track is the eye; trust it and snap
            # onto it, resetting the anchor and breaking the trap.  A noisy or
            # teleporting raw minimum (a transient false low from the secondary
            # circulation) fails the coherence test and we fall through to the
            # fallback as before — so the gate's protection against genuine
            # outliers is retained.
            reacq = reacq_jump_cells if reacq_jump_cells is not None else max_jump_cells
            reacquired = False
            if last_raw_pos is not None:
                raw_jump = np.sqrt(
                    (ix_raw - last_raw_pos[0] / dx) ** 2
                    + (iy_raw - last_raw_pos[1] / dx) ** 2
                )
                if raw_jump <= reacq:
                    x_new, y_new = theta_raw_m   # trust the eye
                    method = 'theta_reacq'
                    reacquired = True

            if not reacquired:
                # θ′ minimum jumped too far AND raw θ′ is not coherent — fall
                # through to vorticity arg-max, Gaussian-localised about
                # last_pos so it reads the primary eye rather than the centroid
                # of the whole vorticity field.
                x_new, y_new = _vorticity_centre(
                    state.u, state.v, nx, ny, dx, xp,
                    last_pos=last_pos, r0_cells=_R0_CELLS,
                )
                method = 'vorticity'

                # Sanity-check vorticity result against last position too
                vx_cells = (x_new / dx) - last_ix
                vy_cells = (y_new / dx) - last_iy
                vdist    = np.sqrt(vx_cells ** 2 + vy_cells ** 2)
                if vdist > max_jump_cells * 2:
                    # Both signals disagree badly — dead reckoning with
                    # background flow.  Secondary intensification (max|u| >
                    # ~100 m/s) can generate strong secondary vorticity
                    # features that corrupt both the θ′ proxy AND the vorticity
                    # CoM simultaneously.  Holding last_pos (old behaviour)
                    # anchors the tracker at the wrong location for the
                    # remainder of the run.  Instead, extrapolate with the
                    # environmental steering so the tracker stays close to the
                    # actual vortex position.
                    #
                    # AXIS CONVENTION (Five, brief 4 P0 correction, May 2026):
                    #   x_c (last_pos[0]) = east-west axis  → driven by u_env
                    #   y_c (last_pos[1]) = north-south axis → driven by v_env
                    # Previous code had these SWAPPED.  The swap produced
                    # exactly the +0.02°N / -0.07°W per step seen in Runs 5 & 6.
                    if v_env_ms != 0.0 or u_env_ms != 0.0:
                        x_new  = last_pos[0] + u_env_ms * dt_output_s
                        y_new  = last_pos[1] + v_env_ms * dt_output_s
                    else:
                        x_new  = last_pos[0]
                        y_new  = last_pos[1]
                    method = 'extrapolated'

    if return_diag:
        return (x_new, y_new), method, {'theta_raw_m': theta_raw_m}
    return (x_new, y_new), method
