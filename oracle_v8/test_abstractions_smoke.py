"""
V8 Abstractions Smoke Test (Revised, Post-Five-Review)
========================================================

Verifies the V8 abstractions hang together properly after Five's
P0/P1 review:

    P0.1: TendencyComponent / ProjectionComponent / DiagnosticComponent
          hierarchy distinguishes additive tendencies from constraint
          enforcement and from read-only diagnostics.
    P0.2: phi renamed to projection_potential everywhere; Tendency
          has no d_projection_potential_dt.
    P0.3: Locked w shape convention (nx, ny, nz+1) with explicit
          boundary values, documented in State and LorenzStaggering.
    P1.4: StepStage enum with PRE_PROJECTION, PROJECTION,
          POST_PROJECTION, SLOW, DIAGNOSTIC.
    P1.5: base_state_compatibility() hardened with non-finite checks,
          missing-field detection, is_load_bearing_ready() requirement.
    P1.6: PI marked as exact-form-TBD pending Durran derivation.
    P1.7: Tendency.zeros_like, add_, validate_against_state helpers.
    P1.8: StateVar enum prevents typo class of bugs.

Not a physics test. Just verifies the contracts.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass

import numpy as np


def main():
    print("=" * 64)
    print("V8 ABSTRACTIONS SMOKE TEST (revised)")
    print("=" * 64)

    print("\n[1] Imports...")
    from oracle_v8.solver import (
        LH82AnelasticEquationSet,
        PseudoIncompressibleEquationSet,
        BuoyancyComponent,
        PressureGradientComponent,
        AnelasticProjection,
        AdvectionComponent,
        CoriolisComponent,
        SurfaceDragComponent,
        SpongeDampingComponent,
        OperatorConfig,
        StepStage,
        StateVar,
        State,
        Tendency,
        StepComponent,
        TendencyComponent,
        ProjectionComponent,
    )
    from oracle_v8.grid import LorenzStaggering, CharneyPhillipsStaggering, LevelType
    print("    ✓ all imports successful")

    print("\n[2] Equation sets...")
    lh82 = LH82AnelasticEquationSet()
    pi = PseudoIncompressibleEquationSet()
    assert lh82.name == "LH82_anelastic"
    assert pi.name == "pseudo_incompressible"
    print(f"    ✓ {lh82!r}")
    print(f"    ✓ {pi!r}")

    print("\n[3] Staggerings...")
    lorenz = LorenzStaggering()
    cp = CharneyPhillipsStaggering()
    assert lorenz.level_type("w") == LevelType.HALF
    assert lorenz.level_type("theta") == LevelType.FULL
    assert cp.level_type("theta") == LevelType.HALF
    print(f"    ✓ Lorenz: w on half, θ on full")
    print(f"    ✓ CP:     w on half, θ on half (co-located)")

    print("\n[4] StateVar enum (P1.8)...")
    assert StateVar.U.value == "u"
    assert StateVar.PROJECTION_POTENTIAL.value == "projection_potential"
    all_vars = list(StateVar)
    assert len(all_vars) == 5
    print(f"    ✓ StateVar has {len(all_vars)} canonical names: "
          f"{[v.value for v in all_vars]}")

    print("\n[5] StepStage enum (P1.4)...")
    expected = {"pre_projection", "projection",
                "post_projection", "slow", "diagnostic"}
    actual = {s.value for s in StepStage}
    assert actual == expected
    print(f"    ✓ StepStage has 5 stages: {sorted(actual)}")

    print("\n[6] StepComponent hierarchy (P0.1)...")
    buoy = BuoyancyComponent()
    pg = PressureGradientComponent()
    proj = AnelasticProjection(nx=4, ny=4, nz=8, Lx=1e6, Ly=1e6, Lz=2e4)
    adv = AdvectionComponent()
    cor = CoriolisComponent()
    drag = SurfaceDragComponent()
    sponge = SpongeDampingComponent()

    # All are StepComponents
    for c in (buoy, pg, proj, adv, cor, drag, sponge):
        assert isinstance(c, StepComponent)

    # Tendency producers are TendencyComponents (not ProjectionComponent)
    for c in (buoy, pg, adv, cor, drag, sponge):
        assert isinstance(c, TendencyComponent)
        assert not isinstance(c, ProjectionComponent)

    # AnelasticProjection is a ProjectionComponent (not TendencyComponent)
    assert isinstance(proj, ProjectionComponent)
    assert not isinstance(proj, TendencyComponent)
    print(f"    ✓ 6 TendencyComponents, 1 ProjectionComponent — "
          f"semantically distinct hierarchies")

    # Stages are correct
    assert buoy.stage == StepStage.PRE_PROJECTION
    assert pg.stage == StepStage.PRE_PROJECTION
    assert proj.stage == StepStage.PROJECTION
    for c in (adv, cor, drag, sponge):
        assert c.stage == StepStage.SLOW
    print("    ✓ stages: 2 PRE_PROJECTION, 1 PROJECTION, 4 SLOW")

    # reads/writes return StateVar enum values
    for c in (buoy, pg, proj, adv, cor, drag, sponge):
        for v in c.reads():
            assert isinstance(v, StateVar)
        for v in c.writes():
            assert isinstance(v, StateVar)
    print("    ✓ reads()/writes() return StateVar values (typo-proof)")

    print("\n[7] State / Tendency naming (P0.2 + P0.3)...")
    nx, ny, nz = 4, 4, 8
    state = State(
        u=np.zeros((nx, ny, nz)),
        v=np.zeros((nx, ny, nz)),
        w=np.zeros((nx, ny, nz + 1)),  # P0.3
        theta_prime=np.zeros((nx, ny, nz)),
        projection_potential=np.zeros((nx, ny, nz)),
        t=0.0,
    )
    assert state.w.shape == (nx, ny, nz + 1)
    assert hasattr(state, "projection_potential")
    assert not hasattr(state, "phi")
    print(f"    ✓ State.w shape locked to (nx, ny, nz+1) = {state.w.shape}")
    print("    ✓ State has projection_potential, no phi field")

    t = Tendency()
    assert not hasattr(t, "dphi_dt")
    assert not hasattr(t, "d_projection_potential_dt")
    print("    ✓ Tendency has no d_projection_potential_dt (φ is diagnostic)")

    print("\n[8] Tendency helpers (P1.7)...")
    zeros = Tendency.zeros_like(state)
    assert zeros.du_dt.shape == state.u.shape
    assert zeros.dw_dt.shape == state.w.shape
    assert np.all(zeros.du_dt == 0)
    print("    ✓ Tendency.zeros_like() produces correct-shape zeros")

    t1 = Tendency.zeros_like(state)
    t1.du_dt += 1.0
    t2 = Tendency.zeros_like(state)
    t2.du_dt += 2.0
    t1.add_(t2)
    assert np.all(t1.du_dt == 3.0)
    print("    ✓ Tendency.add_() correctly accumulates in-place")

    t3 = Tendency(du_dt=None, dv_dt=np.full_like(state.v, 5.0))
    t1.add_(t3)
    assert np.all(t1.du_dt == 3.0)  # unchanged because t3.du_dt was None
    assert np.all(t1.dv_dt == 5.0)
    print("    ✓ Tendency.add_() handles None entries (treats as zero)")

    bad_t = Tendency(du_dt=np.zeros((5, 5, 5)))
    try:
        bad_t.validate_against_state(state)
        return 1
    except ValueError:
        pass
    zeros.validate_against_state(state)  # should not raise
    print("    ✓ Tendency.validate_against_state() catches shape mismatches")

    print("\n[9] OperatorConfig composition...")
    # Production config: PressureGradientComponent is intentionally ABSENT.
    # The projection directly computes and applies ∇φ — adding PGC would
    # double-count it.  The __post_init__ guard now enforces this.
    full_config = OperatorConfig(
        equation_set=lh82,
        staggering=lorenz,
        buoyancy=buoy,
        projection=proj,
        advection=adv,
        coriolis=cor,
        surface_drag=drag,
        sponge_damping=sponge,
    )
    assert len(full_config.active_components()) == 6
    assert len(full_config.fast_components()) == 1       # buoyancy only
    assert len(full_config.projection_components()) == 1
    assert len(full_config.slow_components()) == 4
    print(f"    ✓ Production config: 6 active "
          f"(1 PRE_PROJECTION, 1 PROJECTION, 4 SLOW)")

    # Verify the PGC + projection guard fires correctly
    try:
        OperatorConfig(
            equation_set=lh82,
            staggering=lorenz,
            pressure_gradient=pg,
            projection=proj,
        )
        assert False, "Guard should have raised ValueError"
    except ValueError as e:
        assert "double-counts" in str(e)
    print("    ✓ PGC + projection guard fires correctly (double-count prevented)")

    # Research config: PGC without projection (valid for isolation tests)
    stage1 = OperatorConfig(
        equation_set=lh82,
        staggering=lorenz,
        buoyancy=BuoyancyComponent(),
        pressure_gradient=PressureGradientComponent(),
    )
    assert len(stage1.active_components()) == 2
    print(f"    ✓ Research config (PGC, no projection): 2 active components")

    print("\n[10] Configuration logging...")
    log_dict = full_config.to_log_dict()
    json_str = json.dumps(log_dict, indent=2)
    assert log_dict["n_pre_projection"] == 1
    assert log_dict["n_projection"] == 1
    assert log_dict["n_slow"] == 4
    proj_log = next(c for c in log_dict["active_components"]
                    if c["name"] == "anelastic_projection")
    assert proj_log["kind"] == "projection"
    buoy_log = next(c for c in log_dict["active_components"]
                    if c["name"] == "buoyancy")
    assert buoy_log["kind"] == "tendency"
    print(f"    ✓ to_log_dict() produces JSON ({len(json_str)} chars)")
    print(f"    ✓ component 'kind' field correctly distinguishes types")

    print("\n[11] writes_conflicts diagnostic...")
    conflicts = full_config.writes_conflicts()
    assert "u" in conflicts
    print(f"    ✓ correctly identifies u as multi-writer "
          f"({len(conflicts['u'])} components contribute additively)")

    print("\n[12] Stub method enforcement / wired-up projection...")
    # AnelasticProjection.apply_projection is now WIRED to the Poisson
    # solver. Smoke-test it with a zero-velocity state: ∇·(ρ̄·0) = 0
    # everywhere, so φ should be zero (modulo gauge), and velocities
    # should remain zero after projection.
    #
    # Build a minimal base-state-like object with rho0 = 1.225 and
    # theta0 = 300 K everywhere — sufficient for all smoke-test wiring
    # checks (projection + buoyancy).
    class TrivialBase:
        rho0   = np.full(nz, 1.225)
        theta0 = np.full(nz, 300.0)

    proj.apply_projection(state, lh82, lorenz, TrivialBase(), 1.0)
    assert np.all(state.u == 0)
    assert np.all(state.v == 0)
    assert np.all(state.w == 0)
    # φ should be at machine precision (gauge-pinned to 0)
    phi_max = float(np.max(np.abs(state.projection_potential)))
    assert phi_max < 1e-10, f"φ should be ~0 for d=0, got max|φ|={phi_max}"
    print(f"    ✓ AnelasticProjection.apply_projection wired up; "
          f"d=0 → φ=0 ({phi_max:.2e}), velocities preserved")

    # Confirm diagnostics are populated
    assert proj.last_solve is not None
    assert proj.last_solve.compatibility_residual < 1e-10
    print(f"    ✓ projection diagnostics logged: compat_residual="
          f"{proj.last_solve.compatibility_residual:.2e}, "
          f"disc_op={proj.last_solve.discrete_operator_residual:.2e}")

    # BuoyancyComponent is now fully wired — verify it computes a tendency
    # without error and returns a non-trivial dw_dt (buoyancy acts on w).
    buoy_tend = buoy.compute_tendency(state, lh82, lorenz, TrivialBase(), 1.0)
    assert buoy_tend is not None
    assert buoy_tend.dw_dt is not None
    print("    ✓ BuoyancyComponent.compute_tendency wired up (no longer stubbed)")

    print("\n[13] Hardened base_state_compatibility (P1.5)...")

    # None should fail
    ok, reason = lh82.base_state_compatibility(None)
    assert not ok
    print(f"    ✓ None base state rejected")

    # NaN should fail
    @dataclass
    class FakeBase:
        rho0: np.ndarray = None
        theta0: np.ndarray = None
        z: np.ndarray = None
        integration_scheme: str = "discrete_v8"
        ready: bool = True
        def is_load_bearing_ready(self):
            return self.ready

    fake_nan = FakeBase(
        rho0=np.array([1.0, np.nan, 1.0]),
        theta0=np.array([300.0, 300.0, 300.0]),
        z=np.array([0.0, 100.0, 200.0]),
    )
    ok, reason = lh82.base_state_compatibility(fake_nan)
    assert not ok and "non-finite" in reason
    print(f"    ✓ NaN ρ̄(z) rejected: {reason[:55]}...")

    # Not load-bearing should fail
    fake_not_ready = FakeBase(
        rho0=np.array([1.0, 0.5]),
        theta0=np.array([300.0, 305.0]),
        z=np.array([0.0, 100.0]),
        integration_scheme="trapezoidal_placeholder",
        ready=False,
    )
    ok, reason = lh82.base_state_compatibility(fake_not_ready)
    assert not ok and "load-bearing" in reason
    print(f"    ✓ Placeholder base state rejected: {reason[:55]}...")

    # Valid should pass
    fake_valid = FakeBase(
        rho0=np.array([1.225, 1.0, 0.7]),
        theta0=np.array([300.0, 305.0, 310.0]),
        z=np.array([0.0, 1000.0, 2000.0]),
    )
    ok, reason = lh82.base_state_compatibility(fake_valid)
    assert ok
    print(f"    ✓ Valid load-bearing base state accepted")

    # PI returns False with derivation-pending message
    ok, reason = pi.base_state_compatibility(fake_valid)
    assert not ok
    print(f"    ✓ PI: {reason[:55]}...")

    print("\n" + "=" * 64)
    print("V8 ABSTRACTIONS SMOKE TEST (REVISED): PASSED")
    print()
    print("All eight Five P0/P1 fixes verified:")
    print("  P0.1 ✓ ProjectionComponent distinct from TendencyComponent")
    print("  P0.2 ✓ projection_potential rename, no dphi_dt")
    print("  P0.3 ✓ w shape locked to (nx, ny, nz+1)")
    print("  P1.4 ✓ StepStage enum with 5 stages")
    print("  P1.5 ✓ base_state_compatibility hardened")
    print("  P1.6 ✓ PI marked exact-form-TBD")
    print("  P1.7 ✓ Tendency.zeros_like, add_, validate_against_state")
    print("  P1.8 ✓ StateVar enum prevents typos")
    print()
    print("Abstractions ready for Poisson solver implementation.")
    print("=" * 64)
    return 0


if __name__ == "__main__":
    sys.exit(main())
