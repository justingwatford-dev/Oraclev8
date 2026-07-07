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
        self.kbl  = int(np.sum(zc < 1_500.0))                  # BL levels (z<1.5km)
        self.set_center(LX / 2.0, LY / 2.0)    # f-plane rows: vortex ~stationary

    def set_center(self, cx_m, cy_m):
        """(Re)build the polar geometry about (cx, cy) — moving-frame support
        for the J2 rows (Run 4); f-plane rows keep the domain center."""
        self._cx, self._cy = float(cx_m), float(cy_m)
        x = (xp.arange(self.nx) + 0.5) * self.dx - self._cx
        y = (xp.arange(self.ny) + 0.5) * self.dx - self._cy
        X, Y = xp.meshgrid(x, y, indexing="ij")
        self.r    = xp.sqrt(X ** 2 + Y ** 2)                   # (nx, ny)
        r_safe    = xp.where(self.r > 0, self.r, 1.0)
        self.cosp = X / r_safe
        self.sinp = Y / r_safe
        self._shell = {R: (xp.abs(self.r - R) < self.dx)       # ring masks
                       for R in (300e3, 500e3)}
        self._disk  = {R: (self.r < R) for R in (300e3, 500e3, 800e3)}

    def budget(self, state, drag_comp, cap_comp, cx=None, cy=None,
               u_env=0.0, v_env=0.0):
        """Comparative AM instrumentation (device math, scalars to host).
        With a background flow, all budget fields use the PERTURBATION wind
        (u−u_env) so the uniform steering cannot leak into the rings through
        the gyre asymmetries."""
        if cx is not None:
            self.set_center(cx, cy)
        u = state.u - u_env
        v = state.v - v_env
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
        ll = dg.low_level_vmax(state, self.base_w, self._cx, self._cy,
                               self.dx, self.dx, u_env=u_env, v_env=v_env)
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


def run_row_j2(label, profile_kw, hours):
    """Run 4 (bird-1 mechanism): the J2 rung — β-plane + ERA5-like steering
    ramp, Ivan structure — with the AM budget in the MOVING frame (vortex
    tracked via vorticity_center, geometry recentered at each budget read).
    Tests the supply-chain-delay hypothesis: does the compact profile's BL
    relative-M import at 300 km rise earlier than the envelope's, leading
    each row's own intensification onset?"""
    g = Grid(320, 320, 32)
    ramp = (-1.1, 3.9, +0.6, 6.6)                 # J2's steering ramp
    V0, B0 = 66.9, 1.5                            # Ivan structure (storm_data)
    init = HollandVortexInit(Vmax=V0, Rmax=75_000.0, B=B0, f=F,
                             R_env=500_000.0, u_env=ramp[0], v_env=ramp[1],
                             **profile_kw)
    state = init.build_state(g.nx, g.ny, g.nz, LX, LY, g.base)
    state = dc_replace(state, theta_prime=xp.zeros_like(state.theta_prime))
    pre = RK3Integrator(config=build_prebal_config(g.nx, g.ny, g.nz, LX, LY, LZ),
                        base=g.base)
    for i in range(N_PREBAL):
        state, _ = pre.step(state, dt=1.0, step_number=i)

    drag = SurfaceDragComponent(Cd=CD, H_bl=H_BL,
                                u_env=ramp[0], v_env=ramp[1])
    cor = CoriolisComponent(f=F, mode="beta_plane", Ly=LY, ny=g.ny,
                            u_env=ramp[0], v_env=ramp[1], periodic_taper=True)
    cfg = OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        advection=AdvectionComponent(nx=g.nx, ny=g.ny, nz=g.nz,
                                     Lx=LX, Ly=LY, Lz=LZ),
        coriolis=cor,
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
                                u_env=ramp[0], v_env=ramp[1])

    print(f"\n=== {label}  (J2: beta + ramp {ramp}, Ivan Vmax={V0:.0f}, "
          f"moving-frame budget) ===")
    print(f"  {'t(h)':>5} {'Vmax':>6} {'y(km)':>7} {'res300':>10} "
          f"{'imp300':>10} {'impBL300':>10} {'minBL300':>10} {'impBL500':>10} "
          f"{'Tdrag':>10} {'Tcap':>10}")
    u_envc, v_envc = ramp[0], ramp[1]
    gx, gy = LX / 2.0, LY / 2.0
    b0 = g.budget(state, drag, cap, gx, gy, u_envc, v_envc)
    hist = [(0.0, b0, gy)]
    _p = lambda t, b, y: print(
        f"  {t:5.1f} {b['vmax']:6.1f} {y/1e3:7.0f} {b['res300']:10.3e} "
        f"{b['imp300']:10.3e} {b['impBL300']:10.3e} {b['minBL300']:10.3e} "
        f"{b['impBL500']:10.3e} {b['Tdrag']:10.3e} {b['Tcap']:10.3e}",
        flush=True)
    _p(0.0, b0, gy)

    n_steps = int(hours * 3600.0 / DT)
    diag_every = 60                                # lockstep cadence (30 min)
    budget_every = int(os.environ.get("AMB_BUDGET_EVERY", "240"))
    alpha_steer = min(1.0, (diag_every * DT) / 10800.0)
    for n in range(1, n_steps + 1):
        state, _ = integ.step(state, dt=DT, step_number=n)
        ct = cap.compute_tendency(state, None, None, g.base_w, DT)
        state = dc_replace(state, u=state.u + DT * ct.du_dt,
                           v=state.v + DT * ct.dv_dt)
        if mx(state.u) > 400.0 or not np.isfinite(mx(state.u)):
            print(f"  BLEW UP @step {n} (t={n*DT/3600:.1f} h)")
            return label, hist, None
        if n % diag_every == 0:
            # lockstep steering relaxation toward the ramp target (run_ivan
            # pattern: shift the state background AND every env reference)
            frac = n / n_steps
            u_tgt = ramp[0] + (ramp[2] - ramp[0]) * frac
            v_tgt = ramp[1] + (ramp[3] - ramp[1]) * frac
            du = (u_tgt - u_envc) * alpha_steer
            dv = (v_tgt - v_envc) * alpha_steer
            state = dc_replace(state, u=state.u + du, v=state.v + dv)
            u_envc += du
            v_envc += dv
            cor.set_env(u_envc, v_envc)
            drag.set_env(u_envc, v_envc)
            cap.set_env(u_envc, v_envc)
            gx += u_envc * (diag_every * DT)      # dead-reckon seed
            gy += v_envc * (diag_every * DT)
        if n % budget_every == 0:
            vc = dg.vorticity_center(state, gx, gy, g.dx, g.dx,
                                     base=g.base_w, window_m=250e3)
            gx, gy = vc["xv_m"], vc["yv_m"]        # recenter the seed
            b = g.budget(state, drag, cap, gx, gy, u_envc, v_envc)
            hist.append((n * DT / 3600.0, b, gy))
            _p(n * DT / 3600.0, b, gy)
    vpeak = max(b["vmax"] for _, b, _ in hist)
    onset = next((t for t, b, _ in hist if b["vmax"] >= 60.0
                  and t > 8.0), None)             # skip the init transient
    final = hist[-1][1]
    final["peak"] = vpeak
    final["onset"] = onset
    final["y_end"] = hist[-1][2] / 1e3
    print(f"  END: Vmax {final['vmax']:.1f}  peak {vpeak:.1f}  "
          f"onset(>60, post-dip) {onset if onset is not None else '—'} h  "
          f"y_end {final['y_end']:.0f} km (taper zone > ~4000)")
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
    J2ROWS = dict(
        G=("G J2 compact taper  (Run-3 onset ~26-28h)", dict(wind_taper=True)),
        H=("H J2 gauss r_d=420km (Run-3 onset ~38-40h)",
           dict(outer_envelope_m=420e3)),
    )
    finals = {}
    for key in "ABCDEF":
        if key in which and key in ROWS:
            lbl, nz, prof, norm = ROWS[key]
            _, _, fin = run_row(lbl, nz, prof, norm, hours)
            finals[key] = fin
    j2h = hours if "AMB_HOURS" in os.environ else 52.0     # J2 runs 52 h
    for key in "GH":
        if key in which:
            lbl, prof = J2ROWS[key]
            _, _, fin = run_row_j2(lbl, prof, j2h)
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
    if any(k in finals for k in "GH"):
        print("\nJ2 SUMMARY (Run 4, moving-frame budget):")
        print(f"  {'row':>44} {'peak':>5} {'end':>5} {'onset>60':>8} {'y_end':>6}")
        for key in "GH":
            if key in finals and finals[key] is not None:
                f = finals[key]
                ons = f"{f['onset']:.0f}h" if f["onset"] is not None else "—"
                print(f"  {J2ROWS[key][0]:>44} {f['peak']:5.1f} {f['vmax']:5.1f} "
                      f"{ons:>8} {f['y_end']:6.0f}")
        print("  READ (registered in ENVELOPE_INTENSIFICATION.md Run 4):")
        print("    P-M2: each row's impBL300 rise LEADS its own Vmax onset, and")
        print("          compact impBL300 >= 2x envelope's in the t=12-24h window;")
        print("    P-M3: minBL300 comparable (within ~30%) while impBL differs")
        print("          -> reservoir-LOCATION mechanism;")
        print("          minBL300 differing >= 2x -> inflow-STRENGTH mechanism.")
    print("\nREAD (see ENVELOPE_INTENSIFICATION.md decision rules):")
    print("  bird 1: A vs B/C — spin-up gap + which budget term differs")
    print("  bird 2: D decays, E ≈ F ⇒ drag discretization confirmed+fixed;")
    print("          F vs A = the 0.75→1.00 column-drag effect at nz=32")
    print(f"\nWall time: {time.time()-t0:.0f}s")
