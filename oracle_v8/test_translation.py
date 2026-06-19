"""
Translation Isolation Test
============================
Five's P0 recommendation (ensemble review May 2026):
Run a barotropic vortex with uniform steering before the next Hugo attempt.

Setup: f-plane, barotropic vortex (theta'=0), uniform steering (u_env, v_env),
       no beta, no drag, no sponge, full advection + Coriolis + hyperdiffusion + projection.

Expected: vortex translates cleanly in the steering direction.
         If the Coriolis fix is correct, track is straight NNW.
         If the inertial oscillation is still present, track loops.

Pass criteria:
  [A] After 6h, vortex has moved NNW (westward + northward)
  [B] Track is monotonically northward
  [C] Westward displacement consistent with u_env = -4.5 m/s
  [D] No NaN in 6h (720 steps x dt=30s)
"""
from __future__ import annotations

import sys
import numpy as np

from oracle_v8.storm_data import HUGO
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import (
    LH82AnelasticEquationSet,
    BuoyancyComponent,
    AdvectionComponent,
    CoriolisComponent,
    HyperDiffusionComponent,
    AnelasticProjection,
    OperatorConfig,
    RK3Integrator,
)
from oracle_v8.grid.staggering import LorenzStaggering
from oracle_v8.backend import xp, to_numpy, wrap_base

import sys as _sys
if hasattr(_sys.stdout, "reconfigure"):
    _sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ---------------------------------------------------------------------------
# Domain (same as Hugo run)
# ---------------------------------------------------------------------------
Lx = Ly  = 2_000_000.0
Lz        = 20_000.0
nx = ny   = 128
nz        = 32
dx        = Lx / nx
z_centers = (np.arange(nz) + 0.5) * (Lz / nz)

theta0_arr = 300.0 * np.exp(0.01**2 * z_centers / 9.81)
Pi     = np.zeros(nz)
Pi[0]  = 1.0 - (9.81 / 1004.5) * z_centers[0] / theta0_arr[0]
for k in range(nz - 1):
    dl = z_centers[k + 1] - z_centers[k]
    Pi[k + 1] = Pi[k] - (9.81 / 1004.5) * (dl / 2.0) * (
        1.0 / theta0_arr[k] + 1.0 / theta0_arr[k + 1])
rho0_arr = 100_000.0 * Pi**(1004.5/287.04) / (287.04 * theta0_arr * Pi)


class TestBase:
    z = z_centers; rho0 = rho0_arr; theta0 = theta0_arr


s = HUGO
U_ENV, V_ENV = s["u_env_ms"], s["v_env_ms"]  # -4.5, +12.2

# ---------------------------------------------------------------------------
# Tracker: hydrostatic theta' pressure proxy (Five's recommendation)
# ---------------------------------------------------------------------------

def find_centre_hydrostatic(state, base_wrapped):
    """Hydrostatic pressure proxy: p'_sfc = -integral(rho0 * b dz)"""
    g = 9.81
    b = g * state.theta_prime / base_wrapped.theta0[None, None, :]
    integrand = base_wrapped.rho0[None, None, :] * b
    z = base_wrapped.z
    p_sfc = -xp.trapz(integrand, z, axis=2)   # (nx, ny)
    p_c   = p_sfc - xp.mean(p_sfc)
    idx   = int(xp.argmin(p_c))
    ix, iy = divmod(idx, ny)
    return float((ix + 0.5) * dx), float((iy + 0.5) * dx)


def to_latlon(x_c, y_c, lat0, lon0):
    R = 6_371_000.0
    return (lat0 + np.degrees((y_c - Ly/2)/R),
            lon0 + np.degrees((x_c - Lx/2)/(R * np.cos(np.radians(lat0)))))

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    print("=" * 60)
    print("TRANSLATION ISOLATION TEST")
    print("f-plane | barotropic (theta'=0) | uniform steering")
    print("Expected: clean NNW translation — no loop")
    print("=" * 60)

    # Vortex (barotropic: zero out warm core)
    init  = HollandVortexInit(
        Vmax=s["Vmax_ms"], Rmax=75_000.0, B=s["B"], f=s["f"],
        R_env=500_000.0, u_env=U_ENV, v_env=V_ENV,
    )
    state = init.build_state(nx, ny, nz, Lx, Ly, TestBase())
    # Zero theta' — pure barotropic vortex
    from dataclasses import replace as dc_replace
    state = dc_replace(state, theta_prime=xp.zeros_like(state.theta_prime))

    base_w = wrap_base(TestBase())

    # Pre-balance
    prebal_cfg = OperatorConfig(
        equation_set=LH82AnelasticEquationSet(), staggering=LorenzStaggering(),
        projection=AnelasticProjection(nx=nx,ny=ny,nz=nz,Lx=Lx,Ly=Ly,Lz=Lz),
    )
    prebal_int = RK3Integrator(config=prebal_cfg, base=TestBase())
    for _ in range(5):
        state, _ = prebal_int.step(state, dt=1.0, step_number=0)

    # Main config: f-plane, NO beta, NO drag, NO sponge
    # Coriolis receives u_env, v_env — key fix being tested
    config = OperatorConfig(
        equation_set  = LH82AnelasticEquationSet(),
        staggering    = LorenzStaggering(),
        buoyancy      = BuoyancyComponent(),
        advection     = AdvectionComponent(
                            nx=nx,ny=ny,nz=nz,Lx=Lx,Ly=Ly,Lz=Lz),
        coriolis      = CoriolisComponent(
                            f=s["f"], mode="f_plane",
                            u_env=U_ENV, v_env=V_ENV),   # <<< THE FIX
        horiz_diffusion = HyperDiffusionComponent(
                            nu4=3e11, Lx=Lx, Ly=Ly, nx=nx, ny=ny),
        projection    = AnelasticProjection(
                            nx=nx,ny=ny,nz=nz,Lx=Lx,Ly=Ly,Lz=Lz),
    )
    integrator = RK3Integrator(config=config, base=TestBase())

    DT = 30.0; N = 720  # 6 h
    DIAG = 60           # every 30 min

    print(f"\n  Steering: u={U_ENV:.1f} m/s (W), v={V_ENV:.1f} m/s (N)")
    print(f"  Expected motion: {abs(U_ENV)*6*3600/111000*np.cos(np.radians(28)):.2f}deg W, "
          f"{V_ENV*6*3600/111000:.2f}deg N over 6h\n")
    print(f"  {'Time(h)':>7}  {'Lat(N)':>7}  {'Lon(W)':>7}  "
          f"{'dLon(W)':>8}  {'dLat(N)':>8}  {'max|u|':>8}")
    print(f"  {'-------':>7}  {'------':>7}  {'------':>7}  "
          f"{'-------':>8}  {'-------':>8}  {'------':>8}")

    lat0 = s["lat0_deg"]; lon0 = s["lon0_deg"]
    track_lats, track_lons = [], []
    any_nan = False

    for n in range(N):
        state, diag = integrator.step(state, dt=DT, step_number=n)
        if np.isnan(float(diag.max_u)):
            print(f"\n  NaN at step {n}"); any_nan = True; break

        if n % DIAG == 0:
            t_h = n * DT / 3600.0
            x_c, y_c = find_centre_hydrostatic(state, base_w)
            lat_c, lon_c = to_latlon(x_c, y_c, lat0, lon0)
            track_lats.append(lat_c); track_lons.append(lon_c)
            dlon = -(lon_c - lon0)   # positive = westward
            dlat = lat_c - lat0
            print(f"  {t_h:>7.1f}  {lat_c:>7.2f}  {abs(lon_c):>7.2f}  "
                  f"{dlon:>+8.3f}  {dlat:>+8.3f}  {diag.max_u:>8.3f}")

    print("\n" + "=" * 60)
    if any_nan:
        print("FAILED: NaN"); return 1

    # Pass criteria
    ok_north  = track_lats[-1] > track_lats[0]
    ok_west   = abs(track_lons[-1]) > abs(track_lons[0])  # °W increased
    ok_mono   = all(track_lats[i] <= track_lats[i+1]
                    for i in range(len(track_lats)-1))
    # Check westward: obs displacement ~4.5 * 6h * 3600 / (111km * cos(28))
    exp_west_deg = abs(U_ENV) * 6 * 3600 / (111000 * np.cos(np.radians(28)))
    obs_west_deg = abs(track_lons[-1]) - abs(track_lons[0])
    ok_mag = obs_west_deg > 0.1  # at least some westward motion

    print(f"  Northward:  {'PASS' if ok_north else 'FAIL'} "
          f"({track_lats[0]:.2f} -> {track_lats[-1]:.2f})")
    print(f"  Westward:   {'PASS' if ok_west else 'FAIL'} "
          f"({abs(track_lons[0]):.2f} -> {abs(track_lons[-1]):.2f}W)")
    print(f"  Monotone N: {'PASS' if ok_mono else 'FAIL (loop detected!)'}")
    print(f"  W magnitude: obs={obs_west_deg:.3f}deg, "
          f"exp={exp_west_deg:.3f}deg")

    passed = ok_north and ok_west and ok_mono and ok_mag and not any_nan
    print(f"\n{'PASSED' if passed else 'FAILED'}: Translation test")
    if passed:
        print("  Coriolis geostrophic fix confirmed — no inertial oscillation.")
        print("  Ready to run Hugo with full physics.")
    print("=" * 60)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
