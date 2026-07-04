"""
Angular-momentum budget study — two birds, one instrument (2026-07).

BIRD 1 — the envelope-intensification mechanism (Stage 3 P-S4 flag): under the
Gaussian outer envelope, the two storms with strong barotropic
(re)intensification phases (Hugo −14, Ivan −25 m/s) ran much weaker, while the
four decay/steady storms barely changed.  Hypothesis family: dry barotropic
spin-up feeds on angular-momentum import by the drag-driven boundary-layer
inflow; the profile change alters the mid-radius M supply and/or the inflow
itself.  We reproduce the phenomenon on the f-plane (the intensify-ladder's J0
showed spin-up needs no β/steering) and measure the budget directly.

BIRD 2 — the NZ=64 vortex spin-down (LH82 study Phase 3B open flag): ANALYTIC
result, derived 2026-07-03 from SurfaceDragComponent — the historical
α₀ = Cd·|V′|/dz prefactor over-counts when the (1−z/H_bl) profile spans
multiple levels.  Column factor S = Σ_k max(0, 1−z_k/H_bl):
    NZ=32 (dz=625):    S = 0.6875 + 0.0625            = 0.75
    NZ=64 (dz=312.5):  S = 0.84375+0.53125+0.21875    = 1.59375
    → halving dz multiplies the integrated surface drag by 2.125
    → continuum limit S → H_bl/(2dz): DIVERGENT.
The rows below confirm numerically and test the integral-preserving fix
(SurfaceDragComponent(column_normalized=True) — Σα·dz = Cd|V′| on any grid).

ROWS (all f-plane, u=v=0, init Vmax 64, cap 70, production physics values):
    A  compact taper (cos 200→500 km)   NZ=32   historical drag   ← control
    B  gauss r_d=420 km                 NZ=32   historical drag   ← bird 1
    C  gauss r_d=560 km                 NZ=32   historical drag   ← tail dial
    D  compact                          NZ=64   historical drag   ← bird 2 repro
    E  compact                          NZ=64   NORMALIZED drag   ← bird 2 fix
    F  compact                          NZ=32   NORMALIZED drag   ← fix at prod grid

BUDGET INSTRUMENT (every 2 h): low-level Vmax; relative-AM reservoirs
(r < 300/500/800 km); ring imports of relative M at r = 300/500 km (full
column and boundary-layer z<1.5 km branch); BL mass inflow; drag torque and
cap torque inside 500 km; max|w|; center-offset guard.  Comparative
instrumentation — differences between rows carry the story; no claim of a
closed budget (diffusion/damper/projection sinks land in the residual).

Registered predictions: ENVELOPE_INTENSIFICATION.md (written before any GPU
run).  Usage:  $env:AMB_HOURS=48 ; python -m oracle_v8.measure_am_budget
Optional: $env:AMB_ROWS="ABC" to run a subset.
"""
import os
os.environ.setdefault("ORACLE_GPU", "1")

import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import time
from dataclasses import replace as dc_replace
from types import SimpleNamespace

import numpy as np
from oracle_v8.backend import xp, to_numpy, wrap_base
from oracle_v8 import diagnostics as dg
from oracle_v8.vortex_init import HollandVortexInit
from oracle_v8.solver import (
    RK3Integrator, AdvectionComponent, CoriolisComponent, SurfaceDragComponent,
    IntensityCapComponent, HyperDiffusionComponent,
    HelmholtzDivergenceDampingComponent, NewtonianCoolingComponent,
    LH82AnelasticEquationSet, AnelasticProjection, OperatorConfig)
from oracle_v8.grid.staggering import LorenzStaggering
from oracle_v8.production_config import (
    build_base_state, build_prebal_config, N_PREBAL,
    NU4, EPSILON, TAU_COOL, CD, H_BL, VMAX_CAP_MS, TAU_CAP)

LX = LY = 5_000_000.0          # Ivan-grid domain
LZ = 20_000.0
F  = 5.7e-5                    # f-plane (Ivan latitude)
DT = 30.0
DIAG_EVERY = 240               # budget every 2 h


def mx(a):
    return float(to_numpy(xp.max(xp.abs(a))))


class Grid:
    def __init__(self, nx, ny, nz):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dx = LX / nx
        self.dz = LZ / nz
        zc, rho0, theta0 = build_base_state(nz, LZ)
        self.base   = SimpleNamespace(z=zc, rho0=rho0, theta0=theta0)
        self.base_w = wrap_base(self.base)     # device-side base for direct
                                               # component calls (run_storm pattern)
        self.rho_d  = xp.asarray(rho0)[None, None, :]          # ρ̄(z)
        self.zc     = zc
        # polar geometry about the domain center (f-plane: vortex ~stationary)
        x = (xp.arange(nx) + 0.5) * self.dx - LX / 2.0
        y = (xp.arange(ny) + 0.5) * self.dx - LY / 2.0
        X, Y = xp.meshgrid(x, y, indexing="ij")
        self.r    = xp.sqrt(X ** 2 + Y ** 2)                   # (nx, ny)
        r_safe    = xp.where(self.r > 0, self.r, 1.0)
        self.cosp = X / r_safe
        self.sinp = Y / r_safe
        self.kbl  = int(np.sum(zc < 1_500.0))                  # BL levels (z<1.5km)
        self._shell = {R: (xp.abs(self.r - R) < self.dx)       # ring masks
                       for R in (300e3, 500e3)}
        self._disk  = {R: (self.r < R) for R in (300e3, 500e3, 800e3)}

    def budget(self, state, drag_comp, cap_comp):
        """Comparative AM instrumentation (device math, scalars to host)."""
        u, v = state.u, state.v
        vt = -u * self.sinp[:, :, None] + v * self.cosp[:, :, None]
        ur =  u * self.cosp[:, :, None] + v * self.sinp[:, :, None]
        dV   = self.dx * self.dx * self.dz
        Mrel = self.r[:, :, None] * vt                          # r·v_t
        rhoM = self.rho_d * Mrel

        out = {}
        for R, m in self._disk.items():
            out[f"res{int(R/1e3)}"] = float(to_numpy(
                xp.sum(rhoM[m, :]) * dV))                       # kg m²/s (×1)
        for R, sh in self._shell.items():
            flux_col = xp.sum(self.rho_d * ur * Mrel, axis=2) * self.dz
            flux_bl  = xp.sum((self.rho_d * ur * Mrel)[:, :, :self.kbl],
                              axis=2) * self.dz
            mass_bl  = xp.sum((self.rho_d * ur)[:, :, :self.kbl],
                              axis=2) * self.dz
            n = f"{int(R/1e3)}"
            out[f"imp{n}"]   = float(to_numpy(-xp.mean(flux_col[sh]))) * 2*np.pi*R
            out[f"impBL{n}"] = float(to_numpy(-xp.mean(flux_bl[sh])))  * 2*np.pi*R
            out[f"minBL{n}"] = float(to_numpy(-xp.mean(mass_bl[sh])))  * 2*np.pi*R
        # torques inside 500 km (tangential component of each tendency)
        m5 = self._disk[500e3]
        td = drag_comp.compute_tendency(state, None, None, self.base_w, DT)
        dvt = -td.du_dt * self.sinp[:, :, None] + td.dv_dt * self.cosp[:, :, None]
        out["Tdrag"] = float(to_numpy(
            xp.sum((self.rho_d * self.r[:, :, None] * dvt)[m5, :]) * dV))
        if cap_comp is not None:
            tc = cap_comp.compute_tendency(state, None, None, self.base_w, DT)
            cvt = -tc.du_dt * self.sinp[:, :, None] + tc.dv_dt * self.cosp[:, :, None]
            out["Tcap"] = float(to_numpy(
                xp.sum((self.rho_d * self.r[:, :, None] * cvt)[m5, :]) * dV))
        else:
            out["Tcap"] = 0.0
        # BL / above-BL split of the 500-km reservoir (drag lives in the BL)
        out["res500BL"] = float(to_numpy(
            xp.sum(rhoM[self._disk[500e3], :self.kbl]) * dV))
        # intensity — THREE instruments (run-1 lesson: the k=0 metric reads the
        # drag-drained surface level, not the vortex):
        #   vmax   = production instrument (low_level_vmax, max |V'| in z<3km)
        #   v_sfc  = k=0 max speed (the drag layer itself)
        #   max_u  = max |u-component| anywhere (the LH82-harness instrument,
        #            for direct comparability to the Phase-3B 64→48 reading)
        ll = dg.low_level_vmax(state, self.base_w, LX / 2.0, LY / 2.0,
                               self.dx, self.dx)
        out["vmax"]   = ll["vmax_lowlvl"]
        out["z_vmax"] = ll["z_vmax_m"]
        out["r_vmax"] = ll["r_vmax_km"]                        # center guard ≈ Rmax
        out["v_sfc"]  = float(to_numpy(xp.max(
            xp.sqrt(u[:, :, 0] ** 2 + v[:, :, 0] ** 2))))
        out["max_u"]  = max(mx(u), mx(v))
        out["max_w"]  = mx(state.w)
        return out

    def profile(self, state):
        """Final Vmax(z): max horizontal speed per level (decoupling picture)."""
        sp = xp.sqrt(state.u ** 2 + state.v ** 2)
        ks = [k for k in (0, 1, 2, 3, 4, 6, 8, 12, 16) if k < self.nz]
        return [(self.zc[k], float(to_numpy(xp.max(sp[:, :, k])))) for k in ks]


def run_row(label, nz, profile_kw, normalized_drag, hours):
    g = Grid(320, 320, nz)
    init = HollandVortexInit(Vmax=64.0, Rmax=75_000.0, B=1.5, f=F,
                             R_env=500_000.0, u_env=0.0, v_env=0.0,
                             **profile_kw)
    state = init.build_state(g.nx, g.ny, g.nz, LX, LY, g.base)
    state = dc_replace(state, theta_prime=xp.zeros_like(state.theta_prime))
    pre = RK3Integrator(config=build_prebal_config(g.nx, g.ny, g.nz, LX, LY, LZ),
                        base=g.base)
    for i in range(N_PREBAL):
        state, _ = pre.step(state, dt=1.0, step_number=i)

    drag = SurfaceDragComponent(Cd=CD, H_bl=H_BL, u_env=0.0, v_env=0.0,
                                column_normalized=normalized_drag)
    cfg = OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        advection=AdvectionComponent(nx=g.nx, ny=g.ny, nz=g.nz,
                                     Lx=LX, Ly=LY, Lz=LZ),
        coriolis=CoriolisComponent(f=F, mode="f_plane", u_env=0.0, v_env=0.0),
        horiz_diffusion=HyperDiffusionComponent(nu4=NU4, Lx=LX, Ly=LY,
                                                nx=g.nx, ny=g.ny),
        divergence_damping=HelmholtzDivergenceDampingComponent(
            epsilon=EPSILON, Lx=LX, Ly=LY, nx=g.nx, ny=g.ny),
        newtonian_cooling=NewtonianCoolingComponent(tau=TAU_COOL),
        surface_drag=drag,
        projection=AnelasticProjection(nx=g.nx, ny=g.ny, nz=g.nz,
                                       Lx=LX, Ly=LY, Lz=LZ),
    )
    integ = RK3Integrator(config=cfg, base=g.base)
    cap = IntensityCapComponent(v_cap=VMAX_CAP_MS, tau=TAU_CAP,
                                u_env=0.0, v_env=0.0)

    print(f"\n=== {label}  (nz={nz}, drag={'NORMALIZED' if normalized_drag else 'historical'}) ===")
    print(f"  {'t(h)':>5} {'Vmax':>6} {'z_vx':>5} {'v_sfc':>6} {'max|u|':>6} "
          f"{'res500':>10} {'res500BL':>10} {'imp500':>10} {'impBL500':>10} "
          f"{'Tdrag':>10} {'Tcap':>10} {'max|w|':>7} {'r_vx':>5}")
    b0 = g.budget(state, drag, cap)
    hist = [(0.0, b0)]
    _p = lambda t, b: print(
        f"  {t:5.1f} {b['vmax']:6.1f} {b['z_vmax']:5.0f} {b['v_sfc']:6.1f} "
        f"{b['max_u']:6.1f} {b['res500']:10.3e} {b['res500BL']:10.3e} "
        f"{b['imp500']:10.3e} {b['impBL500']:10.3e} {b['Tdrag']:10.3e} "
        f"{b['Tcap']:10.3e} {b['max_w']:7.2f} {b['r_vmax']:5.0f}", flush=True)
    _p(0.0, b0)

    n_steps = int(hours * 3600.0 / DT)
    for n in range(1, n_steps + 1):
        state, _ = integ.step(state, dt=DT, step_number=n)
        ct = cap.compute_tendency(state, None, None, g.base_w, DT)
        state = dc_replace(state, u=state.u + DT * ct.du_dt,
                           v=state.v + DT * ct.dv_dt)
        if mx(state.u) > 400.0 or not np.isfinite(mx(state.u)):
            print(f"  BLEW UP @step {n} (t={n*DT/3600:.1f} h)")
            return label, hist, None
        if n % DIAG_EVERY == 0:
            b = g.budget(state, drag, cap)
            hist.append((n * DT / 3600.0, b))
            _p(n * DT / 3600.0, b)
    vpeak = max(b["vmax"] for _, b in hist)
    vmin  = min(b["vmax"] for _, b in hist)
    final = hist[-1][1]
    prof = g.profile(state)
    print("  Vmax(z) final: " + "  ".join(f"{z/1e3:.1f}km:{s:.0f}" for z, s in prof))
    print(f"  END: Vmax {final['vmax']:.1f} @z={final['z_vmax']:.0f}m  "
          f"(min {vmin:.1f}, peak {vpeak:.1f})  v_sfc {final['v_sfc']:.1f}  "
          f"max|u| {final['max_u']:.1f}  res500 {final['res500']:.3e}")
    return label, hist, final


if __name__ == "__main__":
    hours = float(os.environ.get("AMB_HOURS", "48"))
    which = os.environ.get("AMB_ROWS", "ABCDEF").upper()
    t0 = time.time()
    print(f"AM-budget study  xp={xp.__name__}  {hours:.0f}h  rows={which}")
    print(f"  analytic drag column factor: S(nz=32)=0.75  S(nz=64)=1.59375  "
          f"ratio 2.125 (see module docstring)")

    ROWS = dict(
        A=("A compact cos 200-500km", 32, dict(wind_taper=True,
                                               taper_start_frac=0.40), False),
        B=("B gauss r_d=420km",       32, dict(outer_envelope_m=420e3), False),
        C=("C gauss r_d=560km",       32, dict(outer_envelope_m=560e3), False),
        D=("D compact, nz=64",        64, dict(wind_taper=True,
                                               taper_start_frac=0.40), False),
        E=("E compact, nz=64, drag NORMALIZED", 64,
           dict(wind_taper=True, taper_start_frac=0.40), True),
        F=("F compact, nz=32, drag NORMALIZED", 32,
           dict(wind_taper=True, taper_start_frac=0.40), True),
    )
    finals = {}
    for key in "ABCDEF":
        if key in which and key in ROWS:
            lbl, nz, prof, norm = ROWS[key]
            _, _, fin = run_row(lbl, nz, prof, norm, hours)
            finals[key] = fin

    print("\n" + "=" * 74)
    print("SUMMARY (final; Vmax = production low_level_vmax, z<3km):")
    print(f"  {'row':>38} {'Vmax_end':>8} {'v_sfc':>6} {'max|u|':>7} {'res500':>11}")
    for key in "ABCDEF":
        if key in finals and finals[key] is not None:
            f = finals[key]
            print(f"  {ROWS[key][0]:>38} {f['vmax']:8.1f} {f['v_sfc']:6.1f} "
                  f"{f['max_u']:7.1f} {f['res500']:11.3e}")
        elif key in finals:
            print(f"  {ROWS[key][0]:>38} {'BLEW':>8}")
    print("\nREAD (see ENVELOPE_INTENSIFICATION.md decision rules):")
    print("  bird 1: A vs B/C — spin-up gap + which budget term differs")
    print("  bird 2: D decays, E ≈ F ⇒ drag discretization confirmed+fixed;")
    print("          F vs A = the 0.75→1.00 column-drag effect at nz=32")
    print(f"\nWall time: {time.time()-t0:.0f}s")
