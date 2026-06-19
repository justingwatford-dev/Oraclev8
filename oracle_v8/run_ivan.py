"""
Oracle V8 — Hurricane Ivan (2004) Track Simulation
==================================================
Third storm of the storm-agnostic validation set.  Initialized from
HURDAT2 (AL092004) at September 14, 2004, 12 UTC (~42.8 h before Gulf
Shores, AL landfall) — the longest lead time of the three storms, and
the only one that RECURVES (zonal motion stalls ~t+30-36, then turns
NNE into the coast).  Cloned from run_katrina.py (V8.6); physics config
identical, domain enlarged for Ivan's longer northward run (see Domain).
Registered predictions for this run: HURDAT2_VERIFICATION.md.

Before first run:  python -m oracle_v8.era5_steering --download --storm ivan

Grid resolution note (paper limitation, documented):
    Ivan observed Rmax = 25 nm (46 km); run uses 75 km (5×dx) for stability.
    Second-order advection requires Rmax >= 4-5×dx without explicit
    diffusion.  Using Rmax_run = 75 km (5×dx) as the numerically
    representative value.  Track is dominated by steering + beta drift,
    not Rmax — this does not systematically bias landfall comparison.

dt = 30s: CFL = 60×30/15625 = 0.115 (60s caused advection instability).
"""
from __future__ import annotations

import sys as _sys
if hasattr(_sys.stdout, "reconfigure"):
    _sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import sys
import time

import numpy as np

from oracle_v8.storm_data import IVAN
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import (
    LH82AnelasticEquationSet,
    BuoyancyComponent,
    AdvectionComponent,
    CoriolisComponent,
    SurfaceDragComponent,
    IntensityCapComponent,
    SpongeDampingComponent,
    HyperDiffusionComponent,
    HelmholtzDivergenceDampingComponent,
    NewtonianCoolingComponent,
    AnelasticProjection,
    OperatorConfig,
    RK3Integrator,
)
from oracle_v8.grid.staggering import LorenzStaggering
from oracle_v8.backend import xp, to_numpy, wrap_base
from oracle_v8.storm_tracker import find_storm_centre as _track_centre
from oracle_v8.era5_steering  import ERA5Steering, update_coriolis_steering
from oracle_v8.landfall_verify import landfall_report
try:                                    # ventilation-flow / real-Vmax diagnostics
    from oracle_v8 import diagnostics as _vdiag   # (optional; safe if absent)
    _HAVE_DIAG = True
except Exception:
    _HAVE_DIAG = False

# ---------------------------------------------------------------------------
# Domain
# ---------------------------------------------------------------------------
# IVAN_DOMAIN: Ivan travels the farthest north of the three storms (init 23.0°N
# → threshold 30.0°N ≈ +779 km; landfall 30.2°N ≈ +801 km).  In the standard
# 4000 km BIG_DOMAIN the reversed-β taper zone starts at y>3200 km; Ivan's centre
# at threshold sits at y≈2779 km, so its 500 km outer circulation reaches
# y≈3279 km — INSIDE the taper during the scored approach.  That is the exact
# contamination mode that invalidated Hugo's small-domain result (and we now
# know the far-field/β-drift representation is THE live physics question — do
# not contaminate it on the third storm).  Enlarge to 5000 km keeping
# dx=15.625 km → nx=ny=320: taper zone moves to y>4000 km, outer edge at
# threshold ≈3779 km → ~220 km margin.  320 is a composite FFT size (2^6·5);
# if the Poisson/FFT path insists on powers of two, fall back to IVAN_DOMAIN=
# False (4000/256) and treat the last ~2 h before threshold as taper-grazed.
# Cost ≈ 1.56× Katrina per step (~3 h wall for 52 h).  Cap must stay ON.
IVAN_DOMAIN = True
if IVAN_DOMAIN:
    Lx = Ly  = 5_000_000.0
    nx = ny   = 320
else:                        # standard BIG_DOMAIN (Hugo/Katrina geometry)
    Lx = Ly  = 4_000_000.0
    nx = ny   = 256
Lz        = 20_000.0
nz        = 32
dx        = Lx / nx        # 15 625 m  (unchanged — dz, dt, CFL all identical)
dz        = Lz / nz        # 625 m
z_centers = (np.arange(nz) + 0.5) * dz

theta0_arr = 300.0 * np.exp(0.01**2 * z_centers / 9.81)
Pi     = np.zeros(nz)
Pi[0]  = 1.0 - (9.81 / 1004.5) * z_centers[0] / theta0_arr[0]
for k in range(nz - 1):
    dl = z_centers[k + 1] - z_centers[k]
    Pi[k + 1] = Pi[k] - (9.81 / 1004.5) * (dl / 2.0) * (
        1.0 / theta0_arr[k] + 1.0 / theta0_arr[k + 1]
    )
p0_arr   = 100_000.0 * Pi ** (1004.5 / 287.04)
rho0_arr = p0_arr / (287.04 * theta0_arr * Pi)


class IvanBase:
    z      = z_centers
    rho0   = rho0_arr
    theta0 = theta0_arr


# ---------------------------------------------------------------------------
# Run parameters
# ---------------------------------------------------------------------------
RMAX_RUN_M = 75_000.0   # m  (5×dx — numerically representative)
R_ENV_M    = 500_000.0  # m  (500 km — with WIND_TAPER, the real vortex-size knob)
TAPER_START_M    = 200_000.0
TAPER_START_FRAC = TAPER_START_M / R_ENV_M   # = 0.50 at R_env = 500 km
# WIND_TAPER (V8.5.1): the bare Holland profile is never bounded — winds run to
# ~10-12 m/s at 1500-2000 km, filling the domain, which inflates the β-drift ~2-3×
# (betadrift: 2.17→~0.8 m/s) and is the Run-14 overshoot driver.  R_env was inert
# (it only shaped the passive θ′).  WIND_TAPER=True brings the wind smoothly to 0 by
# R_ENV_M, so R_ENV_M finally sets the vortex size.  ⚠ STRUCTURAL — changes the
# vortex for every storm; Hugo (validated 48.1 km, untapered) MUST be re-run with the
# same setting before this is the paper config.  False = bit-identical to prior runs.
WIND_TAPER = True
DT         = 30.0        # s
N_STEPS    = 6240        # 52 h — obs crosses 30.0°N at t+42.0h (landfall t+42.83);
                         #        ~10 h margin so a late run (Hugo was +1.2 h) cannot
                         #        miss the threshold on the longest lead time
DIAG_EVERY = 60          # every 30 min
N_PREBAL   = 5
# Hyperdiffusion (nabla^4) replaces 2nd-order Laplacian diffusion.
# Scale selectivity: grid scale (~20 min) vs vortex scale (~7 days).
# nu4 stability limit: dx^4/(64*dt) = 15625^4/(64*30) = 3.1e12 m^4/s
NU4        = 3.0e11     # m^4/s  (10x below stability limit)

# --- ERA5 requirement (V8.6.1) ----------------------------------------------
# A missing ERA5 file used to fall back SILENTLY to constant HURDAT2-motion
# steering — a DIFFERENT steering architecture — which burned a 3 h Ivan run
# (Ivan_run_1) on a non-production config.  Production runs must abort instead
# of quietly downgrading.  Set False only for an explicit frozen-steering A/B.
REQUIRE_ERA5 = True

# --- Time-varying steering (V8.4) -------------------------------------------
# Relax the uniform background flow toward the ERA5 DLM sampled at the storm's
# CURRENT position and time, instead of freezing it at a single value.  The
# frozen track-mean treats an init-time environmental sample as if it followed
# the storm; that error grows with lead time as the vortex translates away from
# where the sample was taken (negligible for Hugo's 32 h, dominant for
# Katrina's ~600 km / 48 h).  Relaxation lets the background EVOLVE without
# hard-overwriting it: each update closes a fraction (Δt_update / TAU_STEER) of
# the gap to the local DLM.  Implementation moves the background flow in the
# state AND the Coriolis/drag reference together (see loop) so the geostrophic
# fix never decouples from the actual flow — the bug that caused the Runs 5–6
# eastward drift when only the reference was updated.
TIME_VARYING_STEER = True
TAU_STEER          = 10800.0   # s — relaxation timescale (3 h).  Lower → tracks
                               # the DLM ramp tighter; too low re-introduces the
                               # hard-overwrite shock.  ~2–4 h is the sane range.

# --- Intensity cap (V8.4.2) -------------------------------------------------
# The f-plane translation test showed the vortex is numerically unstable above
# ~90 m/s (init 120 → peak 372–493, NaN-prone; ε and drag were only crutches).
# Katrina's mid-run max|u|~150 rides that runaway, which corrupts the intensity
# diagnostic AND inflates the outer circulation that drives the reversed-taper
# β-gyre.  This caps the *perturbation* wind |V'| to a physical ceiling so we get
# an honest intensity and an uncontaminated ventilation read (and, if the runaway
# was amplifying the lag, possibly a smaller lag).  Perturbation-relative ⇒ the
# steering flow is untouched and the cap is advanced in lockstep with the steering
# (set_env).  Applied as a gentle per-step relaxation in the loop (forward Euler,
# αΔt<0.1 → stable; near-axisymmetric so it injects ~no divergence, cleaned by the
# next projection).  Set VMAX_CAP_MS=None to disable (recovers the uncapped run).
VMAX_CAP_MS        = 70.0      # m/s ceiling on |V'| (None disables)
TAU_CAP            = 300.0     # s — cap relaxation timescale (tighten to ~150 if
                               # the eyewall sits well above the ceiling)


# ---------------------------------------------------------------------------
# Storm tracking
# ---------------------------------------------------------------------------

def find_storm_centre(state, base_w, nx, ny, dx):
    """Legacy one-shot call — use _track_centre with last_pos for robustness."""
    (x_c, y_c), _ = _track_centre(state, base_w, nx, ny, dx, last_pos=None)
    return x_c, y_c


def centre_to_latlon(x_c, y_c, lat0, lon0, Lx, Ly):
    R = 6_371_000.0
    dlat = np.degrees((y_c - Ly / 2.0) / R)
    dlon = np.degrees((x_c - Lx / 2.0) / (R * np.cos(np.radians(lat0))))
    return lat0 + dlat, lon0 + dlon


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    s = IVAN
    print("=" * 68)
    print(f"ORACLE V8 — HURRICANE {s['name'].upper()} ({s['year']}) TRACK")
    print("=" * 68)
    print(f"  Init:         {s['init_date']}")
    print(f"  Position:     {s['lat0_deg']}°N, {s['lon0_deg']}°W")
    print(f"  Vmax:         {s['Vmax_ms']:.1f} m/s ({s['Vmax_kt']} kt)")
    print(f"  Rmax obs:     {s['Rmax_m']/1000:.1f} km (EBTRK/legacy)")
    print(f"  Rmax run:     {RMAX_RUN_M/1000:.0f} km (5×dx, see docstring)")
    print(f"  Wind taper:   {'ON' if WIND_TAPER else 'OFF'}  (R_env={R_ENV_M/1000:.0f} km)")
    print(f"  Taper-start:  {TAPER_START_M/1000:.0f} km (frac {TAPER_START_FRAC:.2f})")
    print(f"  B:            {s['B']} (frozen)")
    print(f"  Steering:     u={s['u_env_ms']:.1f}, v={s['v_env_ms']:.1f} m/s")
    print(f"  dt:           {DT:.0f}s  (CFL={s['Vmax_ms']*DT/dx:.3f})")
    print(f"  Run:          {N_STEPS} steps = {N_STEPS*DT/3600:.0f} h")
    print(f"  Obs landfall: {s['landfall_lat']}°N, {s['landfall_lon']}°W "
          f"t+{s['landfall_time_h']:.2f}h")

    # ---- ERA5 DLM steering — load FIRST so t=0 values seed the vortex --------
    # IMPORTANT: The TC translation comes from the background flow baked into
    # the initial HollandVortex wind field, NOT from CoriolisComponent.u_env.
    # CoriolisComponent.u_env/v_env only controls the Coriolis inertial fix
    # (prevents the 25.7h oscillation loop) — it does not steer the vortex.
    # Therefore ERA5 DLM must be applied at INITIALIZATION, not mid-run.
    # Mid-run CoriolisComponent updates (below) are still valuable because they
    # keep the inertial fix correct as the storm moves to higher f.
    # For true time-varying steering a mean-flow relaxation tendency is needed
    # (V8.3 target).
    print("\nLoading ERA5 DLM steering winds...")
    try:
        era5     = ERA5Steering.load(storm="ivan")
        era5.print_summary()
        USE_ERA5 = True
        if TIME_VARYING_STEER:
            # Start at the TRUE t=0 environment and let it evolve.  The
            # track-mean was a frozen-steering crutch ("t=0 underrepresents
            # because v varies 1.3→6.1"); with relaxation that crutch is
            # retired — we start at truth and ramp toward the local DLM.
            u_env_t0, v_env_t0 = era5.get_dlm(0.0, s["lat0_deg"], abs(s["lon0_deg"]))
            print(f"  ERA5 steering: ACTIVE  (time-varying — init at t=0 DLM "
                  f"u={u_env_t0:+.2f} v={v_env_t0:+.2f}, relax τ={TAU_STEER/3600:.1f} h "
                  f"toward DLM at model position)")
        else:
            u_env_t0, v_env_t0 = era5.get_track_mean_dlm()
            print(f"  ERA5 steering: ACTIVE  (frozen track-mean DLM — RETIRED "
                  f"mode, kept for A/B only; production = time-varying)")
    except FileNotFoundError as exc:
        if REQUIRE_ERA5:
            print(f"\n  FATAL: {exc}")
            print("  REQUIRE_ERA5=True — refusing to run with fallback constant")
            print("  steering (not the production config).  Download the ERA5")
            print("  file, or set REQUIRE_ERA5=False for an explicit A/B run.")
            return 1
        print(f"  WARNING: {exc}\n  Falling back to HURDAT2 constant steering.")
        u_env_t0 = s["u_env_ms"]
        v_env_t0 = s["v_env_ms"]
        era5     = None
        USE_ERA5 = False

    # ---- Init vortex -------------------------------------------------------
    print("\nInitializing vortex...")
    t0 = time.time()
    init = HollandVortexInit(
        Vmax=s["Vmax_ms"], Rmax=RMAX_RUN_M, B=s["B"], f=s["f"],
        R_env=R_ENV_M, wind_taper=WIND_TAPER, taper_start_frac=TAPER_START_FRAC,
        u_env=u_env_t0, v_env=v_env_t0,   # ERA5 t=0 DLM (or HURDAT2 fallback)
    )
    state = init.build_state(nx, ny, nz, Lx, Ly, IvanBase())
    # Barotropic initialization: zero theta' to prevent BuoyancyComponent
    # from firing on the initial warm core.  Even with BuoyancyComponent
    # removed from config, starting clean avoids any transient adjustment.
    from dataclasses import replace as dc_replace
    state = dc_replace(state, theta_prime=xp.zeros_like(state.theta_prime))
    print(f"  {time.time()-t0:.1f}s  {init.summary()}")

    # ---- Pre-balance: remove grid-scale divergence -------------------------
    # The Holland wind field has large divergence when placed on a discrete
    # Cartesian grid (grid-scale truncation error, ~15 m/s correction).
    # Running N_PREBAL projection-only steps removes this divergence before
    # starting the dynamics.  After 1 step the field is exactly divergence-
    # free; subsequent steps verify convergence.  The resulting state is the
    # nearest divergence-free representation of the Holland vortex.
    print(f"\nPre-balancing ({N_PREBAL} projection-only iterations)...")
    prebal_config = OperatorConfig(
        equation_set = LH82AnelasticEquationSet(),
        staggering   = LorenzStaggering(),
        projection   = AnelasticProjection(
                           nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz),
    )
    prebal_integrator = RK3Integrator(config=prebal_config, base=IvanBase())
    for i in range(N_PREBAL):
        state, pdiag = prebal_integrator.step(state, dt=1.0, step_number=i)
        phi_rms = float(xp.sqrt(xp.mean(state.projection_potential**2)))
        print(f"  iter {i+1}: phi_rms={phi_rms:.2f} Pa  max|u|={pdiag.max_u:.3f} m/s")
    print("  Pre-balance complete — initial divergence removed.")

    # ---- Config ------------------------------------------------------------
    config = OperatorConfig(
        equation_set      = LH82AnelasticEquationSet(),
        staggering        = LorenzStaggering(),
        # BuoyancyComponent: still disabled — testing the two beta-plane fixes
        # in isolation first. Add buoyancy back once 30h stability confirmed.
        advection         = AdvectionComponent(
                                nx=nx,ny=ny,nz=nz,Lx=Lx,Ly=Ly,Lz=Lz),
        coriolis          = CoriolisComponent(
                                f=s["f"], mode="beta_plane", Ly=Ly, ny=ny,
                                u_env=u_env_t0, v_env=v_env_t0,
                                periodic_taper=True),    # FIX 1: f(0)=f(Ly)=f₀
        horiz_diffusion   = HyperDiffusionComponent(
                                nu4=NU4, Lx=Lx, Ly=Ly, nx=nx, ny=ny),
        # Helmholtz divergence damping: correct zero-vorticity implementation.
        # Solves nabla^2(psi)=D via 2D horizontal FFT, then removes epsilon
        # fraction of the purely divergent velocity per SLOW half-step.
        # nabla x nabla(psi) = 0 by definition → ZERO spurious spin-up.
        # (Previous DivergenceDampingComponent used nabla(nabla.u) which
        # contains a vorticity coupling term that drove barotropic spin-up.)
        divergence_damping = HelmholtzDivergenceDampingComponent(
                                epsilon=0.5, Lx=Lx, Ly=Ly, nx=nx, ny=ny),
        # Newtonian cooling: suppresses theta' accumulation from AdvectionComponent's
        # base-state lapse-rate term (-w*dtheta0/dz). tau=1800s (30min) bounds
        # theta' at ~5K — tight enough to eliminate false tracker minima from
        # secondary circulation features while preserving the primary warm-core
        # signal. tau=3600s (Run 3) was not aggressive enough: tracker still
        # stalled at 4-cell grid clusters. Angular momentum intensification
        # (max|u|→177 m/s) is unaffected by cooling — it's a barotropic mechanism.
        newtonian_cooling  = NewtonianCoolingComponent(tau=1800.0),  # 30min — tighter θ' bound (~5K)
        # Surface drag (Run 9 / A/B): caps the inviscid runaway intensification
        # (max|u|→177 m/s) by supplying the boundary-layer angular-momentum sink
        # that was missing in Runs 1–8.  Perturbation-relative form — drags
        # (u−u_env, v−v_env), so it removes momentum from the VORTEX without
        # spinning down the maintained steering flow (consistent with Coriolis
        # and Sponge).  e-folding ~46 min at the eyewall surface (stable
        # explicitly).  If the vortex over-damps (collapse rather than cap),
        # back Cd to ~1e-3 or H_bl to ~750 m.  Set to None to recover Run 8.
        surface_drag       = SurfaceDragComponent(
                                Cd=1.5e-3, H_bl=1000.0,
                                u_env=u_env_t0, v_env=v_env_t0),
        projection        = AnelasticProjection(
                                nx=nx,ny=ny,nz=nz,Lx=Lx,Ly=Ly,Lz=Lz),
    )
    integrator = RK3Integrator(config=config, base=IvanBase())
    base_w     = wrap_base(IvanBase())
    # Intensity cap — applied as a per-step relaxation in the loop (not an
    # OperatorConfig slot), reusing the tested component's compute_tendency.
    cap_comp = (IntensityCapComponent(v_cap=VMAX_CAP_MS, tau=TAU_CAP,
                                      u_env=u_env_t0, v_env=v_env_t0)
                if VMAX_CAP_MS is not None else None)
    if cap_comp is not None:
        print(f"  Intensity cap: ACTIVE  |V'| → {VMAX_CAP_MS:.0f} m/s "
              f"(τ={TAU_CAP:.0f}s, perturbation-relative)")
    # ERA5 was already loaded above — USE_ERA5 / era5 are in scope.

    # ---- Run ---------------------------------------------------------------
    print(f"\n  {'Time(h)':>6}  {'Lat(°N)':>7}  {'Lon(°W)':>7}  "
          f"{'φ_min(Pa)':>12}  {'max|u|(m/s)':>12}")
    print(f"  {'------':>6}  {'-------':>7}  {'-------':>7}  "
          f"{'----------':>12}  {'-----------':>12}")

    track_t, track_lat, track_lon, track_phi = [], [], [], []
    t_run   = time.time()
    any_nan = False
    last_pos = None          # for storm_tracker continuity gate
    last_raw_pos = None      # previous raw-θ′ centre, for Layer-2a re-acquisition
    # Live background steering — starts at init value, relaxed toward the local
    # ERA5 DLM each diagnostic step when TIME_VARYING_STEER is on.  Held in
    # lockstep with the background flow baked into `state` and the Coriolis/drag
    # references (see relaxation block below).
    u_env, v_env = u_env_t0, v_env_t0
    alpha_steer  = min(1.0, (DIAG_EVERY * DT) / TAU_STEER)

    for n in range(N_STEPS):
        state, diag = integrator.step(state, dt=DT, step_number=n)

        # Intensity cap: bleed |V'| above the ceiling back toward it, every step,
        # preserving direction and leaving the steering flow untouched.
        if cap_comp is not None:
            _ct = cap_comp.compute_tendency(state, None, None, base_w, DT)
            state = dc_replace(state, u=state.u + DT * _ct.du_dt,
                               v=state.v + DT * _ct.dv_dt)

        if n % DIAG_EVERY == 0 and np.isnan(float(diag.max_u)):
            print(f"\n  ✗ NaN at step {n} (t={n*DT/3600:.2f}h)")
            any_nan = True
            break

        if n % DIAG_EVERY == 0:
            t_h = n * DT / 3600.0
            (x_c, y_c), trk_method, trk_diag = _track_centre(
                state, base_w, nx, ny, dx,
                last_pos     = last_pos,
                v_env_ms     = v_env,           # LIVE background flow for dead reckoning
                u_env_ms     = u_env,           # (advances with time-varying steering)
                dt_output_s  = DIAG_EVERY * DT,
                return_diag  = True,            # also report raw (unwindowed) θ′ centre
                last_raw_pos = last_raw_pos,    # for Layer-2a re-acquisition
            )
            last_pos     = (x_c, y_c)
            last_raw_pos = trk_diag['theta_raw_m']
            lat_c, lon_c = centre_to_latlon(
                x_c, y_c, s["lat0_deg"], s["lon0_deg"], Lx, Ly)
            track_t.append(t_h);  track_lat.append(lat_c)
            track_lon.append(lon_c); track_phi.append(diag.surface_phi_min)

            # ---- Mean-flow relaxation (time-varying steering) --------------
            # Sample the ERA5 DLM at the storm's CURRENT position and time and
            # relax the background toward it.  The OLD rule here was "do NOT
            # update the Coriolis reference" — but that was because the early
            # attempts moved the *reference* while the state's background stayed
            # frozen, decoupling the geostrophic fix from the actual flow and
            # driving the Runs 5–6 eastward torque.  We avoid that by moving
            # BOTH together: shift the uniform background already present in the
            # state by Δ(u,v), and advance the Coriolis/drag reference by the
            # SAME Δ.  A uniform shift is divergence- and vorticity-free, so the
            # projection ignores it and the vortex structure (the perturbation
            # u − u_env) is untouched — only the translation flow changes.
            if USE_ERA5:
                u_tgt, v_tgt = era5.get_dlm(t_h, lat_c, abs(lon_c))   # DLM target
                if TIME_VARYING_STEER and n > 0:
                    du = (u_tgt - u_env) * alpha_steer
                    dv = (v_tgt - v_env) * alpha_steer
                    state = dc_replace(state, u=state.u + du, v=state.v + dv)
                    u_env += du
                    v_env += dv
                    config.coriolis.set_env(u_env, v_env)
                    if config.surface_drag is not None:
                        config.surface_drag.set_env(u_env, v_env)
                    if cap_comp is not None:
                        cap_comp.set_env(u_env, v_env)
                    # Bracket shows APPLIED background (what steers the storm),
                    # with the local DLM target it is relaxing toward.
                    steer_tag = f" [{u_env:+.1f},{v_env:+.1f}→{u_tgt:+.1f},{v_tgt:+.1f}]"
                else:
                    steer_tag = f" [{u_tgt:+.1f},{v_tgt:+.1f}]"
            else:
                steer_tag = ""

            tag = " [spin-up]" if n < 120 else ""
            _TRK_TAG = {'theta_prime': 'the', 'theta_reacq': 'rcq',
                        'vorticity': 'vor', 'extrapolated': 'ext'}
            trk_tag = f" [{_TRK_TAG.get(trk_method, trk_method[:3])}]"

            # Raw (unwindowed) θ′ centre — diagnostic only.  If raw-θ′ marches
            # north while the chosen fix (e.g. [vor]) parks, the eye is moving
            # and the estimator is lagging — not a physical stall.
            xr, yr = trk_diag['theta_raw_m']
            lat_raw, lon_raw = centre_to_latlon(
                xr, yr, s["lat0_deg"], s["lon0_deg"], Lx, Ly)
            raw_tag = f"  θraw={lat_raw:5.2f}N/{abs(lon_raw):5.2f}W"

            print(f"  {t_h:>6.1f}  {lat_c:>7.2f}  {abs(lon_c):>7.2f}  "
                  f"{diag.surface_phi_min:>12.1f}  {diag.max_u:>12.4e}"
                  f"{steer_tag}{trk_tag}{tag}{raw_tag}")

            # ---- Measured diagnostics (Five "diagnostics first") -------------
            # Real low-level Vmax (vector |V'|, not max|u|-component), and the
            # ventilation flow at the centre.  gyre = vent − background is the
            # DISCRIMINATOR: gyre v ≈ (storm vN − v_env) and vent vN ≈ storm vN
            # ⇒ the flow itself is reduced (β-gyre / reversed taper, Gemini);
            # vent vN ≈ v_env while the storm still lags ⇒ numerical advection
            # lag (Copilot/Five).  Defensive: never breaks the run.
            if _HAVE_DIAG and n > 0:
                try:
                    _dd = _vdiag.compute_diagnostics(
                        state, IvanBase(), x_c, y_c, dx, dx,
                        u_env=u_env, v_env=v_env,
                        theta_xy_m=trk_diag.get("theta_raw_m", (x_c, y_c)))
                    print(f"          {_vdiag.format_diag_line(_dd)}")
                except Exception as _e:
                    pass

    wall = time.time() - t_run
    print(f"\nWall time: {wall:.1f}s")

    if any_nan:
        print("Run FAILED"); return 1

    # ---- Landfall: along/cross-track verification ---------------------------
    # obs_track + threshold_lat come from the storm dict (storm_data.py,
    # HURDAT2-sourced).  See oracle_v8/landfall_verify.py for the metric.
    # threshold_lat = 30.0°N for Ivan (PROVISIONAL — set in storm_data
    # before this first run; obs crosses it at t+42.0h, 87.9°W).
    landfall_report(track_t, track_lat, track_lon, s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
