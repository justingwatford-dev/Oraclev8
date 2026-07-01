"""
LH82 small-perturbation study — Phase 3: tighten the two soft spots before banking.

Phase 2 concluded LH82 holds to θ′/θ̄ ≈ 10% with no visible breakdown, but the
findings doc (caveats 3–4) leaves two loose ends:

PART A — is θ′/θ̄ ≈ 10% reachable with a LIVE circulation?
The 9.74% headline point (Q=2e-2) needed ε=0.5, which crushes the secondary
circulation (max|w| = 0.53 m/s).  The θ′ amplitude is heating-controlled
(θ′_eq ≈ Q·τ_cool), but the DYNAMICS at 10% were sedated.  Re-run Q=2e-2 at
dt=15 with ε=0 (Phase 2 showed halving dt alone cures the ε=0 blow-up at
Q=1.5e-2) and ε=0.1, with the neglected-term diagnostics ON.  If it survives,
the ~10% stability claim rests on a physically alive configuration — and the
continuity ratio there gives a SECOND point on the ratio-vs-θ′ scaling line
(leading-order expectation: ~10%).

PART B — are the Phase 2B ratios at the realistic 4% eyewall grid-converged?
(buoyancy 4.00%, continuity 3.25% at Q=1e-2.)  The heating ring
(width_r = 30 km) is only ~2 cells at dx = 15.6 km.  Repeat Q=1e-2, ε=0 at:
    128x128x32  dx=15.6 km  dt=30   baseline (= Phase 2B)
    128x128x32  dx=15.6 km  dt=15   dt control (isolates dt from resolution)
    256x256x32  dx= 7.8 km  dt=15   horizontal refinement (ring → ~4 cells)
    128x128x64  dz=312.5 m  dt=15   vertical refinement
Converged ≈ both ratios move only ≲10–20% relative across resolutions.

Run on GPU with $env:LH82_STEPS = 1000 (dt=15 runs take 2x steps for the same
8.3 h of physical time).  $env:LH82_PART = "A" | "B" | "AB" (default AB) —
Part B is dominated by the 256² run (~8x the cost of a baseline run).
"""
import os
os.environ.setdefault("ORACLE_GPU", "1")

import sys
if hasattr(sys.stdout, "reconfigure"):     # θ′/θ̄ in prints survives cp1252 pipes
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import time
from types import SimpleNamespace

import numpy as np
from oracle_v8.backend import xp, to_numpy
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import (
    RK3Integrator, BuoyancyComponent, DiabaticHeatingComponent,
    AdvectionComponent, HelmholtzDivergenceDampingComponent)
from oracle_v8.production_config import (
    build_base_state, build_production_config, build_prebal_config, N_PREBAL)

LX = LY = 2_000_000.0        # fixed physical domain (= 128 x 15.625 km)
LZ = 20_000.0
F = 5.7e-5


def mx(a):
    return float(to_numpy(xp.max(xp.abs(a))))


class Grid:
    """Grid-dependent pieces of the Phase 1/2 harnesses, parametrized so the
    same run() works at refined resolutions (same physical domain)."""

    def __init__(self, nx, ny, nz):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dz = LZ / nz
        zc, rho0, theta0 = build_base_state(nz, LZ)
        self.base = SimpleNamespace(z=zc, rho0=rho0, theta0=theta0)
        self.theta0_d = xp.asarray(theta0)[None, None, :]    # θ̄(z)
        self.rho0_d = xp.asarray(rho0)[None, None, :]        # ρ̄(z)
        kx = 2.0 * xp.pi * xp.fft.fftfreq(nx, d=LX / nx)
        ky = 2.0 * xp.pi * xp.fft.fftfreq(ny, d=LY / ny)
        self._ikx = 1j * kx[:, None, None]
        self._iky = 1j * ky[None, :, None]

    def _hdiv(self, fx, fy):
        """Spectral horizontal divergence ∂fx/∂x + ∂fy/∂y (periodic)."""
        d = (self._ikx * xp.fft.fft2(fx, axes=(0, 1))
             + self._iky * xp.fft.fft2(fy, axes=(0, 1)))
        return xp.real(xp.fft.ifft2(d, axes=(0, 1)))

    def _ddz(self, f):
        """Vertical derivative, full levels (centered interior, one-sided ends)."""
        g = xp.empty_like(f)
        g[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2 * self.dz)
        g[:, :, 0]    = (f[:, :, 1] - f[:, :, 0]) / self.dz
        g[:, :, -1]   = (f[:, :, -1] - f[:, :, -2]) / self.dz
        return g

    def neglected_diag(self, state):
        """LH82 neglected-term ratios on a live state (full-level, w interpolated)."""
        tp = state.theta_prime
        u, v, w = state.u, state.v, state.w
        w_c = 0.5 * (w[:, :, :-1] + w[:, :, 1:])
        ratio = tp / self.theta0_d                   # θ′/θ̄ (signed)
        small = float(to_numpy(xp.max(xp.abs(ratio))))
        rho_p = -self.rho0_d * ratio                 # ρ′ = -ρ̄ θ′/θ̄  (dry, const-p)
        D_neg = self._hdiv(rho_p * u, rho_p * v) + self._ddz(rho_p * w_c)  # ∇·(ρ′u)
        D_ret = self._ddz(self.rho0_d * w_c)                               # ∂(ρ̄w)/∂z
        flux_ratio = mx(D_neg) / (mx(D_ret) + 1e-30)
        return small, flux_ratio


def run(g, Q_max, eps, dt, nsteps):
    init = HollandVortexInit(
        Vmax=64.0, Rmax=75_000.0, B=1.5, f=F,
        R_env=500_000.0, wind_taper=True, taper_start_frac=0.40,
        u_env=0.0, v_env=0.0)
    state = init.build_state(g.nx, g.ny, g.nz, LX, LY, g.base)
    pre = RK3Integrator(config=build_prebal_config(g.nx, g.ny, g.nz, LX, LY, LZ),
                        base=g.base)
    for i in range(N_PREBAL):
        state, _ = pre.step(state, dt=1.0, step_number=i)

    cfg = build_production_config(g.nx, g.ny, g.nz, LX, LY, LZ, F, 0.0, 0.0)
    cfg.advection = AdvectionComponent(nx=g.nx, ny=g.ny, nz=g.nz,
                                       Lx=LX, Ly=LY, Lz=LZ, scheme="upwind5h")
    cfg.buoyancy = BuoyancyComponent()
    cfg.diabatic_heating = DiabaticHeatingComponent(
        Q_max=Q_max, r_eyewall=75_000.0, width_r=30_000.0,
        z_peak=5_000.0, width_z=3_000.0,
        nx=g.nx, ny=g.ny, nz=g.nz, Lx=LX, Ly=LY, Lz=LZ)
    cfg.divergence_damping = (None if eps <= 0.0 else
        HelmholtzDivergenceDampingComponent(epsilon=eps, Lx=LX, Ly=LY,
                                            nx=g.nx, ny=g.ny))
    integ = RK3Integrator(config=cfg, base=g.base)

    heartbeat = max(nsteps // 4, 1)
    for n in range(1, nsteps + 1):
        state, _ = integ.step(state, dt=dt, step_number=n)
        if mx(state.u) > 400.0 or not np.isfinite(mx(state.u)):
            return dict(blew=n, tphys=n * dt / 3600.0)
        if n % heartbeat == 0 and n < nsteps:
            # doubles as the equilibration check: θ′/θ̄ should be flat by 50%
            small_now = float(to_numpy(xp.max(xp.abs(state.theta_prime) / g.theta0_d)))
            print(f"      … step {n}/{nsteps}  max|u|={mx(state.u):6.1f}  "
                  f"max(θ′/θ̄)={100*small_now:5.2f}%", flush=True)
    small, flux_ratio = g.neglected_diag(state)
    return dict(blew=None, max_u=mx(state.u), max_w=mx(state.w),
                small=small, flux_ratio=flux_ratio)


def report(label, r, t0):
    mins = (time.time() - t0) / 60.0
    if r["blew"]:
        print(f"  {label:42s} BLEW UP @step {r['blew']} (t={r['tphys']:.1f} h)"
              f"  [{mins:.1f} min]", flush=True)
    else:
        print(f"  {label:42s} max|u|={r['max_u']:6.1f}  max|w|={r['max_w']:5.2f}  "
              f"max(θ′/θ̄)={100*r['small']:5.2f}%  continuity={100*r['flux_ratio']:5.2f}%"
              f"  [{mins:.1f} min]", flush=True)


if __name__ == "__main__":
    N30 = int(os.environ.get("LH82_STEPS", "1000"))      # steps at dt=30 s
    N15 = 2 * N30                                          # dt=15 s, same phys time
    part = os.environ.get("LH82_PART", "AB").upper()
    print(f"LH82 Phase 3  xp={xp.__name__}  N30={N30}  parts={part}", flush=True)

    if "A" in part:
        print("\n=== PART A — θ′/θ̄ ≈ 10% with a LIVE circulation ===")
        print("  (Phase 2 reference: Q=2e-2 ε=0.5 dt=30 → 9.74%, but max|w|=0.53 m/s)")
        gA = Grid(128, 128, 32)
        for label, eps in [("Q=2.0e-2  ε=0.0  dt=15  128x128x32", 0.0),
                           ("Q=2.0e-2  ε=0.1  dt=15  128x128x32", 0.1)]:
            t0 = time.time()
            report(label, run(gA, 2.0e-2, eps, 15.0, N15), t0)

    if "B" in part:
        print("\n=== PART B — grid convergence of the 4% ratios (Q=1e-2, ε=0) ===")
        print("  (Phase 2B reference @128x128x32 dt=30: buoyancy 4.00%, continuity 3.25%)")
        for label, dims, dt, ns in [
                ("Q=1.0e-2  ε=0.0  dt=30  128x128x32 base", (128, 128, 32), 30.0, N30),
                ("Q=1.0e-2  ε=0.0  dt=15  128x128x32 dt/2", (128, 128, 32), 15.0, N15),
                ("Q=1.0e-2  ε=0.0  dt=15  256x256x32 dx/2", (256, 256, 32), 15.0, N15),
                ("Q=1.0e-2  ε=0.0  dt=15  128x128x64 dz/2", (128, 128, 64), 15.0, N15)]:
            t0 = time.time()
            report(label, run(Grid(*dims), 1.0e-2, 0.0, dt, ns), t0)
