"""
LH82 small-perturbation study — Phase 2.

Phase 1 found: eyewall θ′/θ̄ climbs with heating (1% → 4% at Q=1e-2, stable),
then the model BLOWS UP trying to push past ~5% (Q=1.5e-2 @705, 2e-2 @466),
all with the damper OFF (ε=0). Is that blow-up the LH82 small-perturbation
assumption breaking down (PHYSICAL), or a numerical instability of the ε=0 /
strong-forcing setup (NUMERICAL)?

PART A — decisive test: re-run the blow-up cases with the divergence damper
back on (ε=0.1, 0.5) and with half dt. If the blow-up survives every numerical
knob at the SAME θ′/θ̄, it's physical; if damping / smaller dt cures it, numerical.

PART B — neglected-term ratios on the live state (at stable Q=1e-2):
  small parameter  θ′/θ̄            (= LH82 buoyancy relative error)
  density pert.    |ρ′|/ρ̄  = θ′/θ̄  (dry, const-p)
  continuity error max|∇·(ρ′u)| / max|∂(ρ̄w)/∂z|
                   (neglected ρ′ mass-flux divergence vs the RETAINED vertical
                    mass-flux divergence balancing the eyewall convergence)
"""
import os
os.environ.setdefault("ORACLE_GPU", "1")

import numpy as np
from oracle_v8.backend import xp, to_numpy
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import (
    RK3Integrator, BuoyancyComponent, DiabaticHeatingComponent,
    AdvectionComponent, HelmholtzDivergenceDampingComponent)
from oracle_v8.production_config import (
    build_base_state, build_production_config, build_prebal_config, N_PREBAL)

NX = NY = 128
NZ = 32
DX = 15_625.0
LX = LY = NX * DX
LZ = 20_000.0
DZ = LZ / NZ
F = 5.7e-5

zc, rho0_arr, theta0_arr = build_base_state(NZ, LZ)

class Base:
    z = zc; rho0 = rho0_arr; theta0 = theta0_arr

theta0_d = xp.asarray(theta0_arr)[None, None, :]     # θ̄(z)
rho0_d   = xp.asarray(rho0_arr)[None, None, :]        # ρ̄(z)
_kx = 2.0 * xp.pi * xp.fft.fftfreq(NX, d=DX)
_ky = 2.0 * xp.pi * xp.fft.fftfreq(NY, d=DX)
_ikx = 1j * _kx[:, None, None]
_iky = 1j * _ky[None, :, None]


def mx(a):
    return float(to_numpy(xp.max(xp.abs(a))))


def _hdiv(fx, fy):
    """Spectral horizontal divergence ∂fx/∂x + ∂fy/∂y (periodic)."""
    d = _ikx * xp.fft.fft2(fx, axes=(0, 1)) + _iky * xp.fft.fft2(fy, axes=(0, 1))
    return xp.real(xp.fft.ifft2(d, axes=(0, 1)))


def _ddz(f):
    """Vertical derivative, full levels (centered interior, one-sided ends)."""
    g = xp.empty_like(f)
    g[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2 * DZ)
    g[:, :, 0]    = (f[:, :, 1] - f[:, :, 0]) / DZ
    g[:, :, -1]   = (f[:, :, -1] - f[:, :, -2]) / DZ
    return g


def neglected_diag(state):
    """LH82 neglected-term ratios on a live state (full-level, w interpolated)."""
    tp = state.theta_prime
    u, v, w = state.u, state.v, state.w
    w_c = 0.5 * (w[:, :, :-1] + w[:, :, 1:])
    ratio = tp / theta0_d                       # θ′/θ̄ (signed)
    small = float(to_numpy(xp.max(xp.abs(ratio))))
    rho_p = -rho0_d * ratio                      # ρ′ = -ρ̄ θ′/θ̄  (dry, const-p)
    D_neg = _hdiv(rho_p * u, rho_p * v) + _ddz(rho_p * w_c)     # ∇·(ρ′u)
    D_ret = _ddz(rho0_d * w_c)                                   # ∂(ρ̄w)/∂z
    flux_ratio = mx(D_neg) / (mx(D_ret) + 1e-30)
    return small, flux_ratio


def run(Q_max, eps, dt, nsteps, diag=False):
    init = HollandVortexInit(
        Vmax=64.0, Rmax=75_000.0, B=1.5, f=F,
        R_env=500_000.0, wind_taper=True, taper_start_frac=0.40,
        u_env=0.0, v_env=0.0)
    state = init.build_state(NX, NY, NZ, LX, LY, Base())
    pre = RK3Integrator(config=build_prebal_config(NX, NY, NZ, LX, LY, LZ), base=Base())
    for i in range(N_PREBAL):
        state, _ = pre.step(state, dt=1.0, step_number=i)

    cfg = build_production_config(NX, NY, NZ, LX, LY, LZ, F, 0.0, 0.0)
    cfg.advection = AdvectionComponent(nx=NX, ny=NY, nz=NZ, Lx=LX, Ly=LY, Lz=LZ,
                                       scheme="upwind5h")
    cfg.buoyancy = BuoyancyComponent()
    cfg.diabatic_heating = DiabaticHeatingComponent(
        Q_max=Q_max, r_eyewall=75_000.0, width_r=30_000.0,
        z_peak=5_000.0, width_z=3_000.0,
        nx=NX, ny=NY, nz=NZ, Lx=LX, Ly=LY, Lz=LZ)
    cfg.divergence_damping = (None if eps <= 0.0 else
        HelmholtzDivergenceDampingComponent(epsilon=eps, Lx=LX, Ly=LY, nx=NX, ny=NY))
    integ = RK3Integrator(config=cfg, base=Base())

    for n in range(1, nsteps + 1):
        state, _ = integ.step(state, dt=dt, step_number=n)
        if mx(state.u) > 400.0 or not np.isfinite(mx(state.u)):
            return dict(blew=n, tphys=n * dt / 3600.0)
    small, flux_ratio = neglected_diag(state) if diag else (
        float(to_numpy(xp.max(xp.abs(state.theta_prime) / theta0_d))), None)
    return dict(blew=None, max_u=mx(state.u), max_w=mx(state.w),
                small=small, flux_ratio=flux_ratio)


if __name__ == "__main__":
    print(f"Grid {NX}x{NY}x{NZ}  LH82+upwind+heating  xp={xp.__name__}")
    N30 = int(os.environ.get("LH82_STEPS", "1000"))      # dt=30 s
    N15 = 2 * N30                                          # dt=15 s, same phys time

    print("\n=== PART A — is the blow-up physical or numerical? ===")
    print("  (baseline ε=0,dt=30 blew: Q=1.5e-2 @705, Q=2e-2 @466)")
    configs = [
        ("Q=1.5e-2  ε=0.0  dt=30", 1.5e-2, 0.0, 30.0, N30),
        ("Q=1.5e-2  ε=0.1  dt=30", 1.5e-2, 0.1, 30.0, N30),
        ("Q=1.5e-2  ε=0.5  dt=30", 1.5e-2, 0.5, 30.0, N30),
        ("Q=1.5e-2  ε=0.0  dt=15", 1.5e-2, 0.0, 15.0, N15),
        ("Q=2.0e-2  ε=0.5  dt=30", 2.0e-2, 0.5, 30.0, N30),
    ]
    for label, Q, eps, dt, ns in configs:
        r = run(Q, eps, dt, ns)
        if r["blew"]:
            print(f"  {label:26s} BLEW UP @step {r['blew']} (t={r['tphys']:.1f} h)")
        else:
            print(f"  {label:26s} survived  max|u|={r['max_u']:6.1f}  "
                  f"max|w|={r['max_w']:5.2f}  max(θ′/θ̄)={100*r['small']:.2f}%")

    print("\n=== PART B — neglected-term ratios at stable Q=1e-2 ===")
    r = run(1.0e-2, 0.0, 30.0, N30, diag=True)
    if r["blew"]:
        print(f"  Q=1e-2 unexpectedly blew @step {r['blew']}")
    else:
        print(f"  max(θ′/θ̄)  (LH82 buoyancy rel. error)      = {100*r['small']:.2f}%")
        print(f"  max|∇·(ρ′u)| / max|∂(ρ̄w)/∂z|  (continuity) = {100*r['flux_ratio']:.2f}%")
