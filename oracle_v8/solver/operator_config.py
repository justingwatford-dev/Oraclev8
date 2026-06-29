"""
Oracle V8 OperatorConfig
=========================

OperatorConfig is the structured composition that defines which
TendencyComponent objects are active in a given V8 run. It replaces
the boolean-flag soup proposed in the original V8 architecture sketch
(rejected by both Gemini and Five for good reason).

A production V8 run instantiates a config with the full physics:

    config = OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        buoyancy=BuoyancyComponent(),
        pressure_gradient=PressureGradientComponent(),
        projection=AnelasticProjection(),
        advection=AdvectionComponent(),
        coriolis=CoriolisComponent(),
        surface_drag=SurfaceDragComponent(),
        sponge_damping=SpongeDampingComponent(),
    )

A staged validation test instantiates a smaller config with only the
components it isolates:

    # Stage 1: buoyancy and pressure-gradient only
    config_stage1 = OperatorConfig(
        equation_set=LH82AnelasticEquationSet(),
        staggering=LorenzStaggering(),
        buoyancy=BuoyancyComponent(),
        pressure_gradient=PressureGradientComponent(),
        # everything else None
    )

The solver iterates over `config.active_components()` rather than
checking flags. The validation framework calls `config.to_log_dict()`
to record the active configuration in every test result.

This addresses Decision 7-8 of the V8 architecture review: Gemini's
compositional structure, Five's structured logging, and the rejection
of the original boolean-flag interface from both.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from oracle_v8.solver.equation_set import EquationSet
    from oracle_v8.grid.staggering import GridStaggering
    from oracle_v8.solver.tendency import (
        StepComponent,
        TendencyComponent,
        StepStage,
        BuoyancyComponent,
        PressureGradientComponent,
        AnelasticProjection,
        AdvectionComponent,
        CoriolisComponent,
        SurfaceDragComponent,
        SpongeDampingComponent,
    )


@dataclass
class OperatorConfig:
    """
    The structured configuration of a V8 run.

    Two REQUIRED fields:
        equation_set: which dynamical-core equations to integrate
        staggering: which vertical staggering to use

    Eight OPTIONAL component slots, each None by default. A None slot
    means that component is INACTIVE in this run. The solver iterates
    only over non-None components.

    For validation tests, instantiate with only the components needed
    for the stage being tested. For production runs, instantiate with
    all components active.

    All slots use the abstract TendencyComponent base class type. The
    concrete component types (BuoyancyComponent etc.) are subclasses;
    swapping in alternative implementations (e.g. higher-order advection)
    is just a matter of passing a different concrete class instance.
    """

    # Required: the dynamical-core context
    equation_set: "EquationSet"
    staggering: "GridStaggering"

    # Optional component slots
    buoyancy: "BuoyancyComponent | None" = None
    pressure_gradient: "PressureGradientComponent | None" = None
    projection: "AnelasticProjection | None" = None
    advection: "AdvectionComponent | None" = None
    coriolis: "CoriolisComponent | None" = None
    surface_drag: "SurfaceDragComponent | None" = None
    sponge_damping: "SpongeDampingComponent | None" = None
    horiz_diffusion: "HorizontalDiffusionComponent | None" = None
    newtonian_cooling: "NewtonianCoolingComponent | None" = None
    divergence_damping: "DivergenceDampingComponent | None" = None
    diabatic_heating: "DiabaticHeatingComponent | None" = None
    # Future: microphysics, radiation, turbulence closure

    # Research escape hatch — set True ONLY for diagnostic experiments that
    # intentionally activate both PGC and projection to demonstrate the
    # double-counting instability (e.g. test_rk3_hydrostatic_adjustment.py
    # section D.2).  Never set this in production runs.
    _unsafe_pgc_override: bool = False

    def __post_init__(self) -> None:
        """
        Validate configuration integrity on construction.

        Guards against the PGC double-counting failure mode:
        AnelasticProjection solves ∇·(ρ̄∇φ)=∇·(ρ̄u*) and applies
        u_final = u* − ∇φ.  That correction IS the pressure gradient.
        Adding PressureGradientComponent on top returns −∇φ again as a
        tendency, producing alternating-sign divergence and eventual
        instability (confirmed empirically via the discriminating
        dt-scaling test in test_rk3_hydrostatic_adjustment.py).
        """
        if (self.pressure_gradient is not None
                and self.projection is not None
                and not self._unsafe_pgc_override):
            raise ValueError(
                "Invalid V8 configuration: PressureGradientComponent and "
                "AnelasticProjection cannot both be active.  The projection "
                "solves for φ and corrects u_final = u* − ∇φ, which IS the "
                "pressure gradient force.  Enabling PressureGradientComponent "
                "simultaneously double-counts it and will cause instability.\n"
                "  → Set pressure_gradient=None for all production runs.\n"
                "  → Use pressure_gradient only in isolation (no projection) "
                "for incremental-pressure experiments."
            )

    def active_components(self) -> list["StepComponent"]:
        """
        Return the list of currently-active components.

        The solver does its own ordering based on each component's
        stage; this method just collects the non-None ones.
        """
        component_slots = [
            self.buoyancy,
            self.pressure_gradient,
            self.projection,
            self.advection,
            self.coriolis,
            self.surface_drag,
            self.sponge_damping,
            self.horiz_diffusion,
            self.newtonian_cooling,
            self.divergence_damping,
            self.diabatic_heating,
        ]
        return [c for c in component_slots if c is not None]

    def fast_components(self) -> list["StepComponent"]:
        """Return PRE_PROJECTION-stage components (gravity, pressure gradient).

        These are the FAST-mode components in RK3 split-explicit.
        """
        from oracle_v8.solver.tendency import StepStage
        return [c for c in self.active_components()
                if c.stage == StepStage.PRE_PROJECTION]

    def slow_components(self) -> list["StepComponent"]:
        """Return SLOW-stage components (advection, Coriolis, surface, sponge)."""
        from oracle_v8.solver.tendency import StepStage
        return [c for c in self.active_components()
                if c.stage == StepStage.SLOW]

    def projection_components(self) -> list["StepComponent"]:
        """Return PROJECTION-stage components (the constraint enforcement)."""
        from oracle_v8.solver.tendency import StepStage
        return [c for c in self.active_components()
                if c.stage == StepStage.PROJECTION]

    def to_log_dict(self) -> dict:
        """
        Produce a serializable description of this configuration for
        logging in test results. Captures which equation set, staggering,
        and components are active. Used by the validation framework to
        record exactly what dynamical-core configuration ran any given
        test — addressing Five's recommendation that every validation
        result include the active operator configuration.
        """
        return {
            "equation_set": {
                "name": self.equation_set.name,
                "class": type(self.equation_set).__name__,
            },
            "staggering": {
                "name": self.staggering.name,
                "class": type(self.staggering).__name__,
            },
            "active_components": [
                c.to_log_dict() for c in self.active_components()
            ],
            "n_pre_projection": len(self.fast_components()),
            "n_projection": len(self.projection_components()),
            "n_slow": len(self.slow_components()),
        }

    def writes_conflicts(self) -> dict[str, list[str]]:
        """
        Diagnostic: detect cases where multiple components write the
        same prognostic variable. This isn't necessarily an error
        (tendencies are additive — that's how physics composes), but
        it's worth knowing for debugging and configuration auditing.

        Returns a dict mapping variable name → list of components that
        write it (only entries with multiple writers are included).
        """
        writers: dict[str, list[str]] = {}
        for component in self.active_components():
            for var in component.writes():
                # var is a StateVar enum; .value gives its string name.
                key = var.value if hasattr(var, "value") else str(var)
                writers.setdefault(key, []).append(component.name)
        return {var: names for var, names in writers.items() if len(names) > 1}

    def __repr__(self) -> str:
        active = self.active_components()
        return (f"<OperatorConfig: {self.equation_set.name} on "
                f"{self.staggering.name} grid, {len(active)} components>")
