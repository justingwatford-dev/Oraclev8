"""V8 solver: dynamical core, equation sets, step components."""
from oracle_v8.solver.equation_set import (
    EquationSet,
    LH82AnelasticEquationSet,
    PseudoIncompressibleEquationSet,
    ConstraintResidual,
)
from oracle_v8.solver.tendency import (
    StateVar,
    StepStage,
    State,
    Tendency,
    StepComponent,
    TendencyComponent,
    ProjectionComponent,
    DiagnosticComponent,
    BuoyancyComponent,
    PressureGradientComponent,
    AnelasticProjection,
    AdvectionComponent,
    CoriolisComponent,
    SurfaceDragComponent,
    SpongeDampingComponent,
    HorizontalDiffusionComponent,
    HyperDiffusionComponent,
    NewtonianCoolingComponent,
    HelmholtzDivergenceDampingComponent,
)
from oracle_v8.solver.operator_config import OperatorConfig
