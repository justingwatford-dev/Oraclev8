"""
Oracle V8 RK3 Split-Explicit Integrator
=========================================

Implements the Wicker-Skamarock (2002) three-stage Runge-Kutta scheme
for the LH82 anelastic dynamical core.

Stage coefficients: dt/3, dt/2, dt (the "3-2-1" scheme):

    Stage 1: y* = y^n + (dt/3) * F(y^n)
    Stage 2: y** = y^n + (dt/2) * F(y*)
    Stage 3: y^{n+1} = y^n + dt * F(y**)

The defining property of WS-RK3: every stage advances FROM y^n (the
anchor state at the start of the timestep), not from the previous
stage result. Tendencies are EVALUATED at the previous stage result
but APPLIED from y^n. This gives third-order accuracy for the
time-averaged tendency.

At each stage, the step loop is:
    1. Evaluate PRE_PROJECTION tendencies (buoyancy, pressure gradient)
       at the current stage state.
    2. Advance all prognostic variables from y^n by alpha*dt * tendency.
    3. Enforce ∇·(ρ̄u) = 0 via AnelasticProjection.

SLOW-stage components (advection, Coriolis, surface drag, sponge
damping) enter via Strang operator splitting at the full timestep:
apply SLOW for dt/2 before the three fast stages, then dt/2 after.
These are currently stubs; the framework is in place for when they
are implemented.

θ′ evolution: advection is a SLOW component; until AdvectionComponent
is implemented, θ′ is carried unchanged through each step. Calls to
compute_tendency on the buoyancy component will see the same θ′ at
all three stages.

Reference:
    Wicker, L. J., and W. C. Skamarock, 2002: Time-splitting methods
    for elastic models using forward time schemes. Mon. Wea. Rev.,
    130, 2088-2097.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from oracle_v8.backend import xp, wrap_base

# --- Profiling / diag-sync controls (V8.6.3) ----------------------------------
# ORACLE_PROF=1        → per-section wall timers (sync'd) printed every 100 steps
# ORACLE_DIAG_EVERY=N  → compute the 6 host-sync step diagnostics only every N
#                        steps (default 1 = every step, the historical behavior;
#                        set 60 to match DIAG_EVERY in the run scripts).
import os as _os
import time as _time
_PROF       = _os.environ.get("ORACLE_PROF", "0") == "1"
_DIAG_EVERY = max(1, int(_os.environ.get("ORACLE_DIAG_EVERY", "1")))


def _dev_sync():
    """Drain the GPU pipeline so wall timers attribute correctly."""
    try:
        import cupy as _cp
        _cp.cuda.Stream.null.synchronize()
    except Exception:
        pass

from oracle_v8.solver.tendency import State, Tendency, StateVar

if TYPE_CHECKING:
    from oracle_v8.solver.operator_config import OperatorConfig


# --------------------------------------------------------------------------
# Per-step diagnostics
# --------------------------------------------------------------------------


@dataclass
class StepDiagnostics:
    """
    Scalar metrics collected during a single RK3 step.

    Designed to be lightweight: everything is a Python float or a short
    list, so a list of StepDiagnostics over N steps uses negligible
    memory compared to storing full State arrays.

    Fields
    ------
    step : int
        Step index (0-based).
    t : float
        Simulation time at the START of this step (seconds).
    dt : float
        Timestep used for this step (seconds).
    stage_compat_residuals : list of float, length 3
        Compatibility residual (see AnelasticProjection.last_solve)
        from each of the three RK3 sub-stages.  Should be at machine
        precision (< 1e-8) for a well-posed source.
    stage_disc_op_residuals : list of float, length 3
        Discrete-operator residual from each sub-stage.  Should be at
        Thomas-algorithm precision (~1e-12 or better).
    max_u, max_v, max_w : float
        Max absolute value of each velocity component at the END of
        the step.  For stability monitoring: if these exceed O(10)
        m/s on a 10-km grid, check CFL.
    max_theta_prime : float
        Max |θ′| at end of step.  Should be constant while advection
        is a stub (θ′ doesn't evolve yet).
    surface_phi_min : float
        Minimum of gauge-corrected φ on the lowest full level.
        This is the headline diagnostic: the hydrostatic-adjustment
        pressure deficit.  Expected to deepen toward the equilibrium
        value (~-123 Pa for the B&F warm bubble).
    surface_phi_max : float
        Maximum of gauge-corrected φ on the lowest full level.
        Should stay near zero away from the bubble.
    """
    step: int
    t: float
    dt: float
    stage_compat_residuals: list[float] = field(default_factory=list)
    stage_disc_op_residuals: list[float] = field(default_factory=list)
    max_u: float = 0.0
    max_v: float = 0.0
    max_w: float = 0.0
    max_theta_prime: float = 0.0
    surface_phi_min: float = 0.0
    surface_phi_max: float = 0.0


# --------------------------------------------------------------------------
# Integrator
# --------------------------------------------------------------------------


class RK3Integrator:
    """
    Wicker-Skamarock RK3 integrator for the V8 anelastic core.

    Construction
    ------------
        integrator = RK3Integrator(config=config, base=base)

    where `config` is an OperatorConfig with:
        - equation_set, staggering    (required)
        - projection                  (required: AnelasticProjection,
                                       pre-configured with grid params)
        - buoyancy, pressure_gradient (optional: if present, MUST be
                                       configured with grid params)
        - advection, coriolis, etc.   (optional stubs; no-op when absent)

    `base` is any object with `rho0` (1D array, shape (nz,)) and
    `theta0` (1D array) attributes — the hydrostatic base state.

    Usage
    -----
        state, diag = integrator.step(state, dt=5.0, step_number=n)

    The integrator is STATELESS between calls: no internal arrays are
    stored, no hidden step counter.  Pass the step index explicitly so
    the integrator can be checkpointed and resumed without
    re-instantiation.

    Timestep
    --------
    `dt` is passed to each `step()` call, not stored at construction.
    This allows CFL-adaptive timestepping: compute the CFL-stable dt
    from the current state, then call step() with that dt.  The
    integrator adjusts stage_dt = alpha * dt automatically.
    """

    # WS-RK3 stage coefficients: fraction of dt used in each stage advance.
    _STAGE_ALPHAS: tuple[float, ...] = (1.0 / 3.0, 1.0 / 2.0, 1.0)

    def __init__(self, config: "OperatorConfig", base) -> None:
        self.config = config
        # Move base state arrays (z, rho0, theta0) to compute device so
        # every tendency component receives device arrays from self.base.
        self.base = wrap_base(base)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def step(
        self,
        state: State,
        dt: float,
        step_number: int = 0,
    ) -> tuple[State, StepDiagnostics]:
        """
        Advance state by one full RK3 timestep of duration dt.

        Parameters
        ----------
        state : State
            Model state at time t.  Not modified; returns a new State
            at t + dt.
        dt : float
            Timestep in seconds.
        step_number : int
            Step index used only for diagnostics (default 0).

        Returns
        -------
        new_state : State
            Model state at t + dt with ∇·(ρ̄u) = 0 enforced.
        diagnostics : StepDiagnostics
            Scalar metrics from this step.
        """
        diag = StepDiagnostics(step=step_number, t=state.t, dt=dt)

        # ---- SLOW-mode Strang splitting (first half-step) ----
        # _apply_slow is a no-op when no SLOW components are configured.
        if _PROF and not hasattr(self, "_prof"):
            from collections import defaultdict
            self._prof, self._prof_steps = defaultdict(float), 0
        if _PROF:
            _dev_sync(); _t0 = _time.perf_counter()
        state = self._apply_slow(state, dt * 0.5)
        if _PROF:
            _dev_sync(); self._prof["slow_half_1"] += _time.perf_counter() - _t0

        # ---- Three fast-mode RK3 sub-stages ----
        state_n = state
        state_stage = state_n

        if self.config.fast_components():
            for stage_idx, alpha in enumerate(self._STAGE_ALPHAS):
                state_stage, compat, disc_op = self._rk3_stage(
                    state_n=state_n,
                    state_eval=state_stage,
                    alpha=alpha,
                    dt=dt,
                )
                diag.stage_compat_residuals.append(compat)
                diag.stage_disc_op_residuals.append(disc_op)
        else:
            # No fast physics.  V8.6.3 skipped the 3 no-op RK3+projection
            # sub-stages here as redundant — true for a normal run, but for a
            # PROJECTION-ONLY config (the pre-balance) those sub-stages ARE the
            # projection.  Skipping them left the field un-projected and blew up
            # the first real step (φ_min ~7x physical).  Apply ONE projection
            # (not 3 — no residual accumulation) so ∇·(ρ̄u)=0 is enforced.
            # (orig. note: Five, third ensemble review May 2026.)
            state_stage = state_n
            if self.config.projection is not None:
                self.config.projection.apply_projection(
                    state_stage,
                    self.config.equation_set,
                    self.config.staggering,
                    self.base,
                    dt,
                )


        # ---- SLOW-mode Strang splitting (second half-step) ----
        # Capture the projection potential BEFORE the second SLOW step.
        # When fast_components() is non-empty: phi_fast is the physically
        # meaningful buoyancy-driven signal (original design intent).
        # When fast_components() is empty (barotropic runs): phi_fast is
        # near-zero residual from 3 no-op projections; the real signal is
        # phi_slow from the post-SLOW projection above.  We capture
        # state_stage.projection_potential regardless — in the empty-fast
        # case state_stage IS state_n (post-first-SLOW), so phi reflects
        # the actual -βu driven divergence correction.
        phi_fast = state_stage.projection_potential.copy()
        if _PROF:
            _dev_sync(); _t1 = _time.perf_counter()
        state_stage = self._apply_slow(state_stage, dt * 0.5)
        if _PROF:
            _dev_sync(); self._prof["slow_half_2"] += _time.perf_counter() - _t1

        # ---- Collect end-of-step diagnostics ----
        # Pressure: use the FAST-stage φ (the hydrostatic adjustment signal).
        # Velocities and θ′: use the final post-SLOW state.
        if _PROF:
            _dev_sync(); _t2 = _time.perf_counter()
        # V8.6.3: the 6 float() reductions below force a device→host sync each;
        # ORACLE_DIAG_EVERY=N computes them only on steps the caller will read
        # ((step+1) % N == 0 covers run scripts printing after step n-1; % N == 0
        # covers direct users), filling NaN otherwise.
        if (_DIAG_EVERY == 1 or step_number % _DIAG_EVERY == 0
                or (step_number + 1) % _DIAG_EVERY == 0):
            phi_c = phi_fast - xp.mean(phi_fast)
            diag.max_u           = float(xp.max(xp.abs(state_stage.u)))
            diag.max_v           = float(xp.max(xp.abs(state_stage.v)))
            diag.max_w           = float(xp.max(xp.abs(state_stage.w)))
            diag.max_theta_prime = float(xp.max(xp.abs(state_stage.theta_prime)))
            diag.surface_phi_min = float(xp.min(phi_c[:, :, 0]))
            diag.surface_phi_max = float(xp.max(phi_c[:, :, 0]))
        else:
            _nan = float("nan")
            diag.max_u = diag.max_v = diag.max_w = _nan
            diag.max_theta_prime = _nan
            diag.surface_phi_min = diag.surface_phi_max = _nan
        if _PROF:
            _dev_sync(); self._prof["diag_sync"] += _time.perf_counter() - _t2
            self._prof_steps += 1
            if self._prof_steps % 100 == 0:
                tot = sum(self._prof.values())
                parts = "  ".join(f"{k}={v:.1f}s({100*v/tot:.0f}%)"
                                  for k, v in sorted(self._prof.items()))
                print(f"[prof] {self._prof_steps} steps, "
                      f"sum {tot:.1f}s: {parts}", flush=True)

        return state_stage, diag

    # ------------------------------------------------------------------
    # Strang slow-mode operator
    # ------------------------------------------------------------------

    def _apply_slow(self, state: State, dt_slow: float) -> State:
        """
        Apply all SLOW-stage components for dt_slow (one Strang half-step).

        Returns state unchanged when no SLOW components are configured
        (backward compatible with buoyancy+projection-only configs).

        Post-SLOW projection
        --------------------
        If any SLOW component writes momentum fields (U, V, or W), the
        updated velocity may no longer satisfy ∇·(ρ̄u) = 0.  In that
        case a projection step is applied immediately after the SLOW
        advance to restore the constraint before handing state to the
        fast RK3 loop (first half-step) or returning to the caller
        (second half-step).

        AdvectionComponent (V8.1) and CoriolisComponent both write
        momentum, so each Strang half-step includes a projection call.
        The projection with a correctly small SLOW-step dt_slow produces
        a small φ_slow that captures the pressure adjustment driven by
        the advection/Coriolis-induced divergence.

        state.t is not advanced here — simulation time advances in the
        FAST RK3 stages.
        """
        slow_comps = self.config.slow_components()
        if not slow_comps:
            return state

        total = Tendency.zeros_like(state)
        for comp in slow_comps:
            if _PROF:
                _dev_sync(); _tc = _time.perf_counter()
            t = comp.compute_tendency(
                state,
                self.config.equation_set,
                self.config.staggering,
                self.base,
                dt_slow,
            )
            total.add_(t)
            if _PROF:
                _dev_sync()
                self._prof[f"tend_{type(comp).__name__}"] += \
                    _time.perf_counter() - _tc

        state_out = State(
            u             = state.u + dt_slow * total.du_dt,
            v             = state.v + dt_slow * total.dv_dt,
            w             = state.w + dt_slow * total.dw_dt,
            theta_prime   = state.theta_prime + dt_slow * total.dtheta_prime_dt,
            projection_potential = state.projection_potential,
            t             = state.t,
        )

        # Restore ∇·(ρ̄u) = 0 if any SLOW component modified momentum.
        momentum_vars = {StateVar.U, StateVar.V, StateVar.W}
        slow_writes_momentum = any(
            bool(set(comp.writes()) & momentum_vars)
            for comp in slow_comps
        )
        if slow_writes_momentum and self.config.projection is not None:
            if _PROF:
                _dev_sync(); _tp = _time.perf_counter()
            self.config.projection.apply_projection(
                state_out,
                self.config.equation_set,
                self.config.staggering,
                self.base,
                dt_slow,
            )
            if _PROF:
                _dev_sync()
                self._prof["projection"] += _time.perf_counter() - _tp

        return state_out

    # ------------------------------------------------------------------
    # Internal sub-stage logic
    # ------------------------------------------------------------------

    def _rk3_stage(
        self,
        state_n: State,
        state_eval: State,
        alpha: float,
        dt: float,
    ) -> tuple[State, float, float]:
        """
        Execute one WS-RK3 sub-stage with lagged-pressure stability.

        Tendencies are evaluated at state_eval (most up-to-date θ′, u, v, w),
        but the projection_potential passed to PressureGradientComponent is
        taken from state_n — the φ at the START of the full timestep.

        Why lagged pressure?
        --------------------
        If we pass state_eval.projection_potential to the PG component,
        each stage uses φ from the *previous* stage's projection.  For
        dt ≳ 3 s this creates a within-step amplifying feedback:

            Stage 1 → φ₁ (deeper than φ⁰)
            Stage 2 sees PG(φ₁) → larger horizontal inflow
                   → projection deepens to φ₂ > φ₁
            Stage 3 sees PG(φ₂) → still larger inflow → catastrophic φ growth

        Using state_n.projection_potential fixes this: the PG tendency is
        the SAME across all three stages (it's not updated mid-step), so
        the within-step amplification loop is broken.  The scheme becomes
        first-order accurate in the PG but stable for any physically
        reasonable dt.

        This mirrors standard split-explicit NWP practice (WRF, CM1): the
        outer RK3 PG uses the lagged pressure from the previous outer step;
        within-step pressure updates are the acoustic substep's job (not
        applicable here in the anelastic projection formulation).

        Higher accuracy restores when AdvectionComponent is active: θ′ then
        changes between stages, making the multi-stage evaluation genuinely
        useful for buoyancy, and the PG accuracy naturally improves as the
        lagged φ becomes a good predictor of the current-stage φ.

        Parameters
        ----------
        state_n : State
            Anchor: state at the start of the full timestep.
        state_eval : State
            Evaluation point: result of the previous sub-stage.
        alpha : float
            Sub-stage coefficient (1/3, 1/2, or 1).
        dt : float
            Full timestep in seconds.

        Returns
        -------
        state_out : State
        compat_residual : float
        disc_op_residual : float
        """
        stage_dt = alpha * dt

        # 1. Tendency-evaluation state.
        #    θ′ from state_eval (current stage's best estimate).
        #    φ  from state_n    (lagged — see docstring).
        state_for_tendencies = State(
            u=state_eval.u,
            v=state_eval.v,
            w=state_eval.w,
            theta_prime=state_eval.theta_prime,
            projection_potential=state_n.projection_potential,   # ← lagged φ
            t=state_eval.t,
        )

        # 2. Accumulate PRE_PROJECTION tendencies.
        total = Tendency.zeros_like(state_n)
        for comp in self.config.fast_components():
            t = comp.compute_tendency(
                state_for_tendencies,
                self.config.equation_set,
                self.config.staggering,
                self.base,
                stage_dt,
            )
            total.add_(t)

        # 3. Advance FROM state_n.
        state_provisional = State(
            u             = state_n.u + stage_dt * total.du_dt,
            v             = state_n.v + stage_dt * total.dv_dt,
            w             = state_n.w + stage_dt * total.dw_dt,
            theta_prime   = state_n.theta_prime + stage_dt * total.dtheta_prime_dt,
            projection_potential = state_eval.projection_potential.copy(),
            t             = state_n.t + stage_dt,
        )

        # 4. Enforce ∇·(ρ̄u) = 0.
        proj = self.config.projection
        proj.apply_projection(
            state_provisional,
            self.config.equation_set,
            self.config.staggering,
            self.base,
            stage_dt,
        )

        compat   = proj.last_solve.compatibility_residual    if proj.last_solve else 0.0
        disc_op  = proj.last_solve.discrete_operator_residual if proj.last_solve else 0.0

        return state_provisional, compat, disc_op
