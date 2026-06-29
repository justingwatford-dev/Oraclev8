"""
Oracle V8 Validation Test Harness
==================================

Common infrastructure for analytical-solution validation tests.

Each test is a subclass of ValidationTest that implements:
  - initial_state(): build the initial fields
  - analytical_solution(t): the known answer at time t
  - run_step(state, dt): one solver step (will be V8 solver, currently mockable)

The harness handles:
  - Time integration loop
  - L2 / Linf error norm computation against analytical solution
  - Convergence-order analysis across grid resolutions
  - Pass/fail thresholds with reproducible logging
  - JSON-formatted result records (timestamp, git hash, test, error, pass/fail, runtime)

Design intent: every test run produces a record sufficient for the AIES paper's
validation table. "We ran the suite at every commit" should be a verifiable
fact, not a claim.

References for the analytical solutions are documented in analytical_solutions.py
and must be cross-checked against primary sources before the suite is treated
as authoritative.
"""

from __future__ import annotations

import abc
import json
import subprocess
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


# --- Result records ------------------------------------------------------------


@dataclass
class TestResult:
    """A single test's outcome, written to disk for reproducibility."""

    test_name: str
    timestamp_utc: str
    git_commit: str | None
    grid_shape: tuple[int, ...]
    n_steps: int
    dt: float
    runtime_seconds: float
    error_l2: float
    error_linf: float
    pass_threshold_l2: float
    passed: bool
    notes: str = ""
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class ConvergenceResult:
    """Order-of-accuracy result across multiple grid resolutions."""

    test_name: str
    timestamp_utc: str
    git_commit: str | None
    grid_resolutions: list[tuple[int, ...]]
    errors_l2: list[float]
    measured_order: float
    expected_order: float
    passed: bool
    notes: str = ""


# --- Utilities -----------------------------------------------------------------


def get_git_commit() -> str | None:
    """Return the current git commit hash, or None if not in a git repo."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=2.0,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return None


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def l2_error(numerical: np.ndarray, analytical: np.ndarray,
             mean_subtract: bool = True) -> float:
    """Discrete L2 norm of the difference, normalized by analytical magnitude.

    mean_subtract: if True, subtract the spatial mean from both fields before
        comparison. Required for problems whose operator has a constant
        nullspace (e.g. periodic Poisson with all-Neumann or all-periodic
        boundaries: any constant added to φ leaves both ∇φ and ∇·(ρ₀∇φ)
        unchanged, so the solver's solution is determined only up to an
        additive constant). Default True because the most common harness
        use case is elliptic-problem verification.
    """
    if mean_subtract:
        numerical = numerical - np.mean(numerical)
        analytical = analytical - np.mean(analytical)
    diff = numerical - analytical
    norm_diff = float(np.sqrt(np.mean(diff**2)))
    norm_ref = float(np.sqrt(np.mean(analytical**2)))
    if norm_ref < 1e-30:
        return norm_diff  # absolute error if reference is ~zero
    return norm_diff / norm_ref


def linf_error(numerical: np.ndarray, analytical: np.ndarray,
               mean_subtract: bool = True) -> float:
    """Discrete Linf norm of the difference, normalized by analytical max.

    mean_subtract: see l2_error. Default True for the same reason.
    """
    if mean_subtract:
        numerical = numerical - np.mean(numerical)
        analytical = analytical - np.mean(analytical)
    diff = numerical - analytical
    norm_diff = float(np.max(np.abs(diff)))
    norm_ref = float(np.max(np.abs(analytical)))
    if norm_ref < 1e-30:
        return norm_diff
    return norm_diff / norm_ref


# --- Base class ----------------------------------------------------------------


class ValidationTest(abc.ABC):
    """
    Base class for an analytical-solution validation test.

    A test defines an initial state, a known analytical solution as a function
    of time, and a comparison protocol. The harness handles running the solver,
    measuring error, and logging the result.

    Subclasses must implement:
        - initial_state(grid)
        - analytical_solution(grid, t)
        - field_to_compare(state)   # which scalar/vector field is the test on
    """

    name: str = "unnamed_test"
    expected_convergence_order: float = 2.0   # most tests are 2nd order in space
    pass_threshold_l2: float = 1e-3            # default; subclasses override

    @abc.abstractmethod
    def initial_state(self, grid: dict) -> dict:
        """Return the initial field dictionary for the solver."""

    @abc.abstractmethod
    def analytical_solution(self, grid: dict, t: float) -> dict:
        """Return the analytical solution fields at time t."""

    @abc.abstractmethod
    def field_to_compare(self, state: dict) -> np.ndarray:
        """Extract the scalar/vector field that the test compares on."""

    def make_grid(self, shape: tuple[int, ...], domain: dict) -> dict:
        """Build a uniform Cartesian grid descriptor.

        domain: dict with keys 'Lx', 'Ly', 'Lz' (and 2D drops 'Lz').
        """
        if len(shape) == 2:
            nx, ny = shape
            Lx, Ly = domain["Lx"], domain["Ly"]
            x = np.linspace(0, Lx, nx, endpoint=False)
            y = np.linspace(0, Ly, ny, endpoint=False)
            X, Y = np.meshgrid(x, y, indexing="ij")
            return {"shape": shape, "x": x, "y": y, "X": X, "Y": Y,
                    "dx": Lx / nx, "dy": Ly / ny, "domain": domain}
        elif len(shape) == 3:
            nx, ny, nz = shape
            Lx, Ly, Lz = domain["Lx"], domain["Ly"], domain["Lz"]
            x = np.linspace(0, Lx, nx, endpoint=False)
            y = np.linspace(0, Ly, ny, endpoint=False)
            z = np.linspace(0, Lz, nz, endpoint=False)
            X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
            return {"shape": shape, "x": x, "y": y, "z": z,
                    "X": X, "Y": Y, "Z": Z,
                    "dx": Lx / nx, "dy": Ly / ny, "dz": Lz / nz,
                    "domain": domain}
        else:
            raise ValueError(f"Unsupported grid shape: {shape}")

    def run(
        self,
        solver_step_fn,
        grid_shape: tuple[int, ...],
        domain: dict,
        n_steps: int,
        dt: float,
        log_dir: Path | None = None,
    ) -> TestResult:
        """
        Run the test once at fixed grid resolution.

        solver_step_fn: callable(state, dt) -> state, the V8 solver step.
                        For pre-V8 / harness-only testing, can be an identity
                        function or analytical-substitution function.
        """
        grid = self.make_grid(grid_shape, domain)
        state = self.initial_state(grid)

        t_start = time.perf_counter()
        for step in range(n_steps):
            state = solver_step_fn(state, dt)
        runtime = time.perf_counter() - t_start

        t_final = n_steps * dt
        analytical = self.analytical_solution(grid, t_final)

        num_field = self.field_to_compare(state)
        ana_field = self.field_to_compare(analytical)

        err_l2 = l2_error(num_field, ana_field)
        err_linf = linf_error(num_field, ana_field)
        passed = err_l2 < self.pass_threshold_l2

        result = TestResult(
            test_name=self.name,
            timestamp_utc=utc_timestamp(),
            git_commit=get_git_commit(),
            grid_shape=tuple(grid_shape),
            n_steps=n_steps,
            dt=dt,
            runtime_seconds=runtime,
            error_l2=err_l2,
            error_linf=err_linf,
            pass_threshold_l2=self.pass_threshold_l2,
            passed=passed,
        )

        if log_dir is not None:
            self._write_record(result, log_dir)

        return result

    def run_convergence(
        self,
        solver_step_fn,
        grid_shapes: list[tuple[int, ...]],
        domain: dict,
        n_steps: int,
        dt_at_coarsest: float,
        log_dir: Path | None = None,
    ) -> ConvergenceResult:
        """
        Run the test at multiple grid resolutions and measure the convergence
        order. dt is scaled with the grid to keep CFL constant.
        """
        errors = []
        for shape in grid_shapes:
            # CFL-preserving timestep: dt scales with dx
            ratio = grid_shapes[0][0] / shape[0]
            dt_this = dt_at_coarsest * ratio
            n_steps_this = int(round(n_steps / ratio))
            result = self.run(
                solver_step_fn, shape, domain, n_steps_this, dt_this,
                log_dir=None,  # only log the convergence summary
            )
            errors.append(result.error_l2)

        # Order = log(err_coarse / err_fine) / log(refinement_ratio)
        # Use the finest two for the measured order.
        if len(errors) >= 2 and errors[-1] > 0:
            ratio = grid_shapes[-1][0] / grid_shapes[-2][0]
            measured_order = float(np.log(errors[-2] / errors[-1]) / np.log(ratio))
        else:
            measured_order = float("nan")

        passed = (
            not np.isnan(measured_order)
            and measured_order > 0.8 * self.expected_convergence_order
        )

        result = ConvergenceResult(
            test_name=f"{self.name}_convergence",
            timestamp_utc=utc_timestamp(),
            git_commit=get_git_commit(),
            grid_resolutions=[tuple(s) for s in grid_shapes],
            errors_l2=errors,
            measured_order=measured_order,
            expected_order=self.expected_convergence_order,
            passed=passed,
        )

        if log_dir is not None:
            self._write_record(result, log_dir)

        return result

    @staticmethod
    def _write_record(record, log_dir: Path) -> None:
        log_dir.mkdir(parents=True, exist_ok=True)
        ts_safe = record.timestamp_utc.replace(":", "-")
        path = log_dir / f"{record.test_name}_{ts_safe}.json"
        with path.open("w") as f:
            json.dump(asdict(record), f, indent=2)


def compute_anelastic_residual(
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    rho0: np.ndarray,
    dx: float,
    dy: float,
    dz: float,
    periodic: tuple[bool, bool, bool] = (True, True, False),
) -> dict:
    """
    Compute the anelastic constraint residual ∇·(ρ̄ u) using the SOLVER'S
    discrete operator — not a re-derived stencil.

    Delegates to ``oracle_v8.solver.tendency.anelastic_divergence``, the exact
    operator ``AnelasticProjection`` inverts (spectral horizontal divergence +
    flux-form vertical on the staggered half-level w).  This guarantees the
    residual a test reports is the constraint the solver actually enforces.
    The previous implementation re-derived a 2Δz centered stencil that
    disagreed with the solver by ~50× on the same field and so could mask (or
    fabricate) constraint error.

    Parameters
    ----------
    u, v : ndarray, shape (nx, ny, nz)
        Full-level horizontal velocity.
    w : ndarray
        Vertical velocity.  Accepts the solver's staggered half-level form
        (nx, ny, nz+1) directly, or full-level (nx, ny, nz), in which case it
        is interpolated to half levels (interior 0.5-average; rigid w=0 at the
        surface and lid faces).
    rho0 : ndarray, shape (nz,)
        Base-state density ρ̄(z).  The solver operator assumes a horizontally
        homogeneous base state; a 3-D rho0 is not supported (raises).
    dx, dy, dz : float
        Grid spacings (uniform).
    periodic : tuple of bool
        Retained for backward compatibility.  The solver operator is periodic
        in x, y (spectral) and rigid in z (flux form); this argument does not
        change the operator.

    Returns
    -------
    dict with keys: max_abs_residual, l2_residual, residual_field.
    """
    from oracle_v8.solver.tendency import anelastic_divergence
    from oracle_v8.backend import xp, asarray, to_numpy

    nx, ny, nz = u.shape
    if rho0.ndim != 1:
        raise ValueError(
            "compute_anelastic_residual delegates to the solver operator, "
            "which assumes a horizontally homogeneous base state ρ̄(z); got "
            f"rho0 with shape {rho0.shape}. Pass the 1-D (nz,) profile."
        )

    # Wavenumbers in the solver's convention (2π · fftfreq).
    kx = 2.0 * xp.pi * xp.fft.fftfreq(nx, d=dx)
    ky = 2.0 * xp.pi * xp.fft.fftfreq(ny, d=dy)

    rho_full = asarray(rho0)
    rho_half = xp.zeros(nz + 1, dtype=rho_full.dtype)
    rho_half[1:-1] = 0.5 * (rho_full[:-1] + rho_full[1:])
    rho_half[0] = rho_full[0]
    rho_half[-1] = rho_full[-1]

    u_d, v_d = asarray(u), asarray(v)
    if w.shape[-1] == nz + 1:
        w_half = asarray(w)
    elif w.shape[-1] == nz:
        w_d = asarray(w)
        w_half = xp.zeros((nx, ny, nz + 1), dtype=w_d.dtype)
        w_half[:, :, 1:-1] = 0.5 * (w_d[:, :, :-1] + w_d[:, :, 1:])
        # surface/lid faces stay 0 (rigid)
    else:
        raise ValueError(
            f"w last dimension must be nz ({nz}) or nz+1 ({nz + 1}); got {w.shape}"
        )

    residual = to_numpy(
        anelastic_divergence(u_d, v_d, w_half, rho_full, rho_half, kx, ky, dz)
    )

    return {
        "max_abs_residual": float(np.max(np.abs(residual))),
        "l2_residual": float(np.sqrt(np.mean(residual**2))),
        "residual_field": residual,
    }
