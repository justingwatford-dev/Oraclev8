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

from oracle_v8.solver.tendency import State, Tendency

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
        self.base = base

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
        # When SLOW components are implemented, apply here for dt/2.
        # For now: no-op.

        # ---- Three fast-mode RK3 sub-stages ----
        # state_n: anchor — all stage advances depart from here.
        # state_stage: evaluation point — tendencies evaluated here.
        state_n = state
        state_stage = state_n

        for stage_idx, alpha in enumerate(self._STAGE_ALPHAS):
            state_stage, compat, disc_op = self._rk3_stage(
                state_n=state_n,
                state_eval=state_stage,
                alpha=alpha,
                dt=dt,
            )
            diag.stage_compat_residuals.append(compat)
            diag.stage_disc_op_residuals.append(disc_op)

        # ---- SLOW-mode Strang splitting (second half-step) ----
        # When SLOW components are implemented, apply here for dt/2.

        # ---- Collect end-of-step diagnostics ----
        phi = state_stage.projection_potential
        phi_c = phi - np.mean(phi)   # gauge-subtract for physical interpretation
        diag.max_u          = float(np.max(np.abs(state_stage.u)))
        diag.max_v          = float(np.max(np.abs(state_stage.v)))
        diag.max_w          = float(np.max(np.abs(state_stage.w)))
        diag.max_theta_prime = float(np.max(np.abs(state_stage.theta_prime)))
        diag.surface_phi_min = float(np.min(phi_c[:, :, 0]))
        diag.surface_phi_max = float(np.max(phi_c[:, :, 0]))

        return state_stage, diag

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
        Execute one WS-RK3 sub-stage.

        Tendencies are evaluated at state_eval; the advance departs
        from state_n.  This is the WS-RK3 "evaluate late, advance from
        anchor" property that gives third-order accuracy.

        Parameters
        ----------
        state_n : State
            Anchor: the state at the start of the full timestep.
        state_eval : State
            Evaluation point: the result of the previous sub-stage
            (or state_n for the first sub-stage).
        alpha : float
            Sub-stage coefficient (1/3, 1/2, or 1).
        dt : float
            Full timestep in seconds (sub-stage dt = alpha * dt).

        Returns
        -------
        state_out : State
            Post-projection state for this sub-stage.
        compat_residual : float
        disc_op_residual : float
            Projection solver diagnostics for this sub-stage.
        """
        stage_dt = alpha * dt

        # 1. Accumulate PRE_PROJECTION (fast-mode) tendencies at state_eval.
        total = Tendency.zeros_like(state_n)
        for comp in self.config.fast_components():
            t = comp.compute_tendency(
                state_eval,
                self.config.equation_set,
                self.config.staggering,
                self.base,
                stage_dt,
            )
            total.add_(t)

        # 2. Advance FROM state_n (not from state_eval).
        #    θ′ is carried from state_n unchanged until AdvectionComponent
        #    is implemented.  When AdvectionComponent is active, it will
        #    return a non-zero dtheta_prime_dt, and the update below picks
        #    it up automatically.
        state_provisional = State(
            u             = state_n.u + stage_dt * total.du_dt,
            v             = state_n.v + stage_dt * total.dv_dt,
            w             = state_n.w + stage_dt * total.dw_dt,
            theta_prime   = state_n.theta_prime + stage_dt * total.dtheta_prime_dt,
            # Warm-start: give the Poisson solver the previous φ as a hint.
            # For a direct (Thomas-algorithm) solver this makes no difference;
            # it matters when we switch to an iterative solver.
            projection_potential = state_eval.projection_potential.copy(),
            t             = state_n.t + stage_dt,
        )

        # 3. Enforce the anelastic constraint ∇·(ρ̄u) = 0.
        proj = self.config.projection
        proj.apply_projection(
            state_provisional,
            self.config.equation_set,
            self.config.staggering,
            self.base,
            stage_dt,
        )

        compat   = proj.last_solve.compatibility_residual   if proj.last_solve else 0.0
        disc_op  = proj.last_solve.discrete_operator_residual if proj.last_solve else 0.0

        return state_provisional, compat, disc_op
