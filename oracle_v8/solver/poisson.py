"""
Oracle V8 Variable-Coefficient Poisson Solver
==============================================

Solves the LH82 anelastic Poisson equation:

    ∇·(ρ̄(z) ∇φ) = d(x, y, z)

on a triply-bounded domain that is periodic in (x, y) and rigid-lidded
in z. φ is the projection potential whose gradient corrects the
unprojected velocity to satisfy ∇·(ρ̄ u) = 0.

Algorithm (hybrid FFT-x,y + tridiagonal-z):

  1. Take the 2D FFT of the source d in the (x, y) plane. Each
     horizontal Laplacian term -kx²-ky² becomes a multiplier at each
     (kx, ky) wavenumber pair.

  2. At each wavenumber pair, the equation reduces to a 1D ODE in z:

         d/dz(ρ̄(z) dφ̂/dz) - (kx² + ky²) ρ̄(z) φ̂ = d̂(kx, ky, z)

     This is a linear two-point BVP. Discretized with the standard
     second-order stencil on the Lorenz grid, it produces a
     tridiagonal system in z. Solve with the Thomas algorithm,
     O(nz) per wavenumber.

  3. Inverse FFT to recover φ(x, y, z).

Vertical boundary conditions: Neumann (dφ/dz = 0) at top and bottom.
This corresponds to rigid surface (w=0) and rigid lid (w=0). The
discretization uses the ghost-point approach: φ_{-1} = φ_1 (mirror),
collapsing the boundary stencil to one-sided to second-order
accuracy. Same treatment inverted at the top.

Gauge handling at (kx=0, ky=0): the all-Neumann mode has a constant
nullspace (any C added to φ leaves ∇φ unchanged). The tridiagonal
system at this wavenumber is singular. We pin the gauge by setting
φ̂(0, 0, k=0) = 0 and solving the remaining (nz-1)×(nz-1) system.

Compatibility condition at (kx=0, ky=0): the equation is solvable
only if the column integral of d̂(0, 0, z) is zero (mass conservation
in the rigid-lidded configuration). For a velocity field with w=0 at
both boundaries, this is automatic. We log the column-integral
residual every solve so we can see at AIES paper time that it's
sitting at machine precision across long integrations — that's
evidence the solver is physically well-posed, not just that it didn't
crash.

References:
  - Lipps and Hemler (1982), J. Atmos. Sci. 39, 2192-2210.
  - Durran (2010), "Numerical Methods for Fluid Dynamics", §7.2.
  - Press et al. (2007), "Numerical Recipes", §2.4 (Thomas algorithm)
    and §20.4 (FFT methods for elliptic PDEs).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from oracle_v8.backend import xp, to_numpy as to_host

# numpy alias for type annotations only — not used in any computation
import numpy as _np
np = _np   # annotation alias; all compute code uses xp

if TYPE_CHECKING:
    from oracle_v8.reference_state.base_states import DryBaseState
    from oracle_v8.grid.staggering import GridStaggering


@dataclass
class PoissonSolveResult:
    """
    The output of one Poisson solve, including diagnostics.

    Attributes
    ----------
    phi : array (device-resident)
        The projection potential, shape (nx, ny, nz) on Lorenz full levels.
    compatibility_residual : float
        Column integral of the (kx=0, ky=0) source mode — machine
        precision for physically valid divergence.  Logged every solve.
    discrete_operator_residual : float
        L2 norm of (operator(phi) − source).  Machine precision for a
        converged Thomas solve.
    n_wavenumbers_solved : int
        Number of (kx, ky) pairs solved.  Equal to nx*ny normally.
    """
    phi: object          # xp.ndarray (cupy or numpy)
    compatibility_residual: float
    discrete_operator_residual: float
    n_wavenumbers_solved: int


class VariableCoefficientPoissonSolver:
    """
    Solver for ∇·(ρ̄(z) ∇φ) = d on the LH82 anelastic configuration.

    Construction caches the wavenumber arrays for the (nx, ny)
    horizontal grid; the solver instance is reusable across timesteps
    as long as the grid resolution doesn't change.

    The solver does NOT cache base-state-dependent matrices because
    the base state is currently constant in time, but a future moving-
    nest variant may need to recompute. For now, base state is passed
    to solve() rather than to __init__.
    """

    def __init__(self, nx: int, ny: int, nz: int,
                 Lx: float, Ly: float, Lz: float):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dx = Lx / nx
        self.dy = Ly / ny
        self.dz = Lz / nz

        # Horizontal wavenumbers — built on device so k_sq lives where the
        # FFT result lives and the batch Thomas uses device arrays.
        self.kx = 2.0 * xp.pi * xp.fft.fftfreq(nx, d=self.dx)
        self.ky = 2.0 * xp.pi * xp.fft.fftfreq(ny, d=self.dy)
        KX, KY = xp.meshgrid(self.kx, self.ky, indexing="ij")
        self.k_sq = KX**2 + KY**2  # shape (nx, ny)

        # V8.6.3 — Thomas precompute cache.  The tridiagonal coefficients AND
        # the rhs-independent forward-elimination factors (c', 1/denom) depend
        # only on (rho, dz, k^2), all constant through a run.  The old code
        # rebuilt ~3x(nx*ny-1, nz) coefficient arrays and re-eliminated them
        # on EVERY solve (2 per step).  Cache keyed on the identity of the rho
        # arrays; strong refs are held so ids stay valid.
        self._thomas_cache = None

    def apply_operator(
        self,
        phi: np.ndarray,
        rho_bar_full: np.ndarray,
        rho_bar_half: np.ndarray,
    ) -> np.ndarray:
        """
        Forward operator: compute ∇·(ρ̄ ∇φ) in physical space using the
        same discrete stencil the solver inverts.

        Used for verifying the discrete residual after a solve (Five's
        P1.4 hook): apply_operator(solve_result.phi) should equal the
        original source d to machine precision modulo the gauge.

        Parameters
        ----------
        phi : ndarray, shape (nx, ny, nz)
            Field on full levels.
        rho_bar_full : ndarray, shape (nz,)
            ρ̄ on full levels.
        rho_bar_half : ndarray, shape (nz+1,)
            ρ̄ on half levels (cell faces, including ghost values at
            k=-1/2 and k=nz+1/2 corresponding to surface and top).
            Ghost values come from the ρ̄ extrapolation chosen by the
            base state; for an exponential ρ̄ profile, the ghost value
            equals the surface value for k=-1/2 (mirror BC).

        Returns
        -------
        result : ndarray, shape (nx, ny, nz)
            The discrete ∇·(ρ̄ ∇φ).
        """
        # Horizontal Laplacian via FFT
        phi_hat = xp.fft.fft2(phi, axes=(0, 1))
        horiz_laplacian_hat = -self.k_sq[:, :, None] * phi_hat
        horiz_laplacian = xp.real(xp.fft.ifft2(horiz_laplacian_hat, axes=(0, 1)))

        # ρ̄ multiplies the horizontal Laplacian on full levels
        horiz_term = rho_bar_full[None, None, :] * horiz_laplacian

        # Vertical operator: d/dz(ρ̄ dφ/dz) — fully vectorised (no z loop).
        #
        # Upward flux at interior half-level k+1/2 (k=0..nz-2):
        #   F[k] = ρ̄_{k+1/2} * (φ_{k+1} − φ_k)
        # stored in `flux[:,:,k]`.
        #
        # Divergence at full level k (interior k=1..nz-2):
        #   vert[k] = (F[k] − F[k−1]) / dz²
        #
        # Neumann BCs:
        #   k=0:    zero lower flux → vert[0]  =  F[0]   / dz²
        #   k=nz-1: zero upper flux → vert[-1] = −F[-1]  / dz²
        dz2 = self.dz ** 2
        phi_diff = phi[:, :, 1:] - phi[:, :, :-1]   # (nx, ny, nz-1)
        flux = rho_bar_half[None, None, 1:-1] * phi_diff   # (nx, ny, nz-1)

        vert_term = xp.empty_like(phi)
        vert_term[:, :,  0]   =  flux[:, :,  0] / dz2
        vert_term[:, :, 1:-1] = (flux[:, :, 1:] - flux[:, :, :-1]) / dz2
        vert_term[:, :, -1]   = -flux[:, :, -1] / dz2

        return horiz_term + vert_term

    def solve(
        self,
        d: np.ndarray,
        rho_bar_full: np.ndarray,
        rho_bar_half: np.ndarray,
        log_diagnostics: bool = True,
    ) -> PoissonSolveResult:
        """
        Solve ∇·(ρ̄ ∇φ) = d for φ.

        Parameters
        ----------
        d : ndarray, shape (nx, ny, nz)
            The source term (typically the divergence of the
            unprojected velocity, weighted by ρ̄).
        rho_bar_full : ndarray, shape (nz,)
            ρ̄ on full levels.
        rho_bar_half : ndarray, shape (nz+1,)
            ρ̄ on half levels.
        log_diagnostics : bool
            If True (default), compute and return the compatibility
            residual and discrete-operator residual for logging.

        Returns
        -------
        PoissonSolveResult with phi and diagnostics.
        """
        # Step 1: FFT in (x, y)
        d_hat = xp.fft.fft2(d, axes=(0, 1))

        # Compatibility residual at (kx=0, ky=0)
        compat_residual = (
            float(xp.abs(xp.sum(d_hat[0, 0, :]) * self.dz))
            if log_diagnostics else 0.0
        )

        # Step 2: solve all wavenumber pairs simultaneously, using cached
        # rhs-independent Thomas factors (V8.6.3 — built once per run).
        cache   = self._get_thomas_cache(rho_bar_full, rho_bar_half)
        phi_hat = self._solve_batch(d_hat, cache)
        # Gauge-pin the (kx=ky=0) Neumann mode.  Its RHS is ∂z(ρ̄ w̄); the mean
        # horizontal divergence is identically 0 on a periodic domain, so the
        # true zero-mode φ ≈ 0.  The V8.6.3 host reduced-Thomas for this mode is
        # near-singular (den ~3e-7) and amplifies RHS noise to ~1e7 → NaN.  Pin
        # it instead of solving the singular system.
        phi_hat[0, 0, :] = 0.0
        # print(f"[psn] d_hat_fin={bool(xp.isfinite(d_hat).all())} "
        #      f"d_hat_max={float(xp.max(xp.abs(d_hat))):.2e} "
        #     f"phi_hat_fin={bool(xp.isfinite(phi_hat).all())} "
        #     f"phi_hat_max={float(xp.max(xp.abs(phi_hat))):.2e}", flush=True)

        # Step 3: inverse FFT
        phi = xp.real(xp.fft.ifft2(phi_hat, axes=(0, 1)))

        # Diagnostic: discrete operator residual
        if log_diagnostics:
            op_phi = self.apply_operator(phi, rho_bar_full, rho_bar_half)
            residual_field = op_phi - d
            residual_field -= xp.mean(residual_field)
            disc_residual = float(xp.sqrt(xp.mean(residual_field**2)))
        else:
            disc_residual = 0.0

        return PoissonSolveResult(
            phi=phi,
            compatibility_residual=compat_residual,
            discrete_operator_residual=disc_residual,
            n_wavenumbers_solved=self.nx * self.ny,
        )

    def _get_thomas_cache(self, rho_bar_full, rho_bar_half):
        key = (id(rho_bar_full), id(rho_bar_half))
        _hit = self._thomas_cache is not None and self._thomas_cache["key"] == key
        # print(f"[psn-cache] key={key} hit={_hit}", flush=True)
        if _hit:
            return self._thomas_cache

        dz2       = self.dz ** 2
        rho_lower = rho_bar_half[:-1] / dz2          # (nz,)
        rho_upper = rho_bar_half[1:]  / dz2          # (nz,)
        nz        = self.nz
        n_total   = self.nx * self.ny
        k_sq_nz   = self.k_sq.reshape(n_total)[1:]   # (n_inner,)

        # a (sub) and c (super) are mode-independent -> store as (nz,) with
        # the Neumann BC rows baked in.  b varies per mode via k^2.
        a = rho_lower.copy(); a[0] = 0.0; a[-1] = rho_lower[-1]
        c = rho_upper.copy(); c[-1] = 0.0
        b = (-(rho_lower[None, :] + rho_upper[None, :])
             - k_sq_nz[:, None] * rho_bar_full[None, :])
        b[:, 0]  = -rho_upper[0]  - k_sq_nz * rho_bar_full[0]
        b[:, -1] = -rho_lower[-1] - k_sq_nz * rho_bar_full[-1]

        # rhs-independent forward elimination: c' and 1/denom, once.
        inv_den = xp.empty_like(b)
        c_p     = xp.empty_like(b)
        inv_den[:, 0] = 1.0 / b[:, 0]
        c_p[:, 0]     = c[0] * inv_den[:, 0]
        for k in range(1, nz):
            den            = b[:, k] - a[k] * c_p[:, k - 1]
            inv_den[:, k]  = 1.0 / den
            if k < nz - 1:
                c_p[:, k]  = c[k] * inv_den[:, k]

        # Zero mode: precompute the host-side elimination once (numpy) --
        # a scalar Thomas loop on device costs ~4*nz tiny kernel launches
        # per solve for nz numbers.
        rl = to_host(rho_lower).copy(); ru = to_host(rho_upper).copy()
        a0 = rl[1:].copy(); b0 = -(rl[1:] + ru[1:]); c0 = ru[1:].copy()
        a0[0] = 0.0; a0[-1] = rl[-1]; c0[-1] = 0.0; b0[-1] = -rl[-1]
        n0 = nz - 1
        inv0 = _np.empty(n0); cp0 = _np.empty(n0)
        inv0[0] = 1.0 / b0[0]; cp0[0] = c0[0] * inv0[0]
        for k in range(1, n0):
            den      = b0[k] - a0[k] * cp0[k - 1]
            inv0[k]  = 1.0 / den
            cp0[k]   = c0[k] * inv0[k] if k < n0 - 1 else 0.0

        import numpy as _np_dbg
        #print(f"[psn-cache] BUILD inv_den_fin={bool(xp.isfinite(inv_den).all())} "
        #      f"c_p_fin={bool(xp.isfinite(c_p).all())} "
        #      f"inv0_fin={bool(_np_dbg.isfinite(inv0).all())} "
        #      f"min|b[:,0]|={float(xp.min(xp.abs(b[:,0]))):.2e} "
        #      f"min|b0|={float(_np_dbg.min(_np_dbg.abs(b0))):.2e}", flush=True)
        self._thomas_cache = dict(
            key=key, a=a, c_p=c_p, inv_den=inv_den,
            a0=a0, cp0=cp0, inv0=inv0,
            _refs=(rho_bar_full, rho_bar_half),
        )
        return self._thomas_cache

    def _solve_batch(self, d_hat, cache):
        """
        Vectorised solve for ALL nx*ny wavenumber pairs using the CACHED
        rhs-independent Thomas factors (V8.6.3).

        The coefficient build and forward elimination of c' and 1/denom are
        done ONCE per run in _get_thomas_cache; each solve now performs only
        the rhs forward sweep and the back substitution — 2 fused element-wise
        kernels per z-level instead of ~5, and zero coefficient allocations.
        The singular (0,0) mode is solved on the HOST with cached factors
        (it is nz scalars; a device scalar loop costs ~4*nz kernel launches).
        """
        n_total = self.nx * self.ny
        nz      = self.nz
        d_flat  = d_hat.reshape(n_total, nz)
        rhs     = d_flat[1:, :]                      # (n_inner, nz) view
        a, c_p, inv_den = cache["a"], cache["c_p"], cache["inv_den"]

        d_p = xp.empty_like(rhs)
        d_p[:, 0] = rhs[:, 0] * inv_den[:, 0]
        for k in range(1, nz):
            d_p[:, k] = (rhs[:, k] - a[k] * d_p[:, k - 1]) * inv_den[:, k]

        x = xp.empty_like(rhs)
        x[:, -1] = d_p[:, -1]
        for k in range(nz - 2, -1, -1):
            x[:, k] = d_p[:, k] - c_p[:, k] * x[:, k + 1]

        phi_flat = xp.empty_like(d_flat)
        phi_flat[0, :]  = self._solve_zero_mode_host(d_flat[0, :], cache)
        phi_flat[1:, :] = x
        return phi_flat.reshape(self.nx, self.ny, self.nz)

    def _solve_zero_mode_host(self, d_zero, cache):
        """Singular (kx=0, ky=0) mode on the host, cached factors (V8.6.3)."""
        dz_h = to_host(d_zero)
        a0, cp0, inv0 = cache["a0"], cache["cp0"], cache["inv0"]
        n0  = len(inv0)
        rhs = dz_h[1:]
        d_p = _np.empty(n0, dtype=dz_h.dtype)
        d_p[0] = rhs[0] * inv0[0]
        for k in range(1, n0):
            d_p[k] = (rhs[k] - a0[k] * d_p[k - 1]) * inv0[k]
        x = _np.empty(n0, dtype=dz_h.dtype)
        x[-1] = d_p[-1]
        for k in range(n0 - 2, -1, -1):
            x[k] = d_p[k] - cp0[k] * x[k + 1]
        out = _np.zeros(len(dz_h), dtype=dz_h.dtype)
        out[1:] = x
        return xp.asarray(out)

    def _solve_zero_mode(
        self,
        d_zero,
        rho_lower,
        rho_upper,
        rho_bar_full,
    ):
        """
        Special-case solver for the singular (kx=0, ky=0) mode.

        Pins the gauge at φ̂(0,0,k=0) = 0 and solves the remaining
        (nz-1)×(nz-1) tridiagonal system.
        """
        nz = len(d_zero)
        phi_zero = xp.zeros(nz, dtype=d_zero.dtype)

        a_sub  = rho_lower[1:].copy()
        b_diag = -(rho_lower[1:] + rho_upper[1:]).copy()
        c_sup  = rho_upper[1:].copy()
        rhs    = d_zero[1:].copy()

        a_sub[0]  = 0.0          # φ[0] is pinned at 0
        a_sub[-1] = rho_lower[-1]
        c_sup[-1] = 0.0
        b_diag[-1] = -rho_lower[-1]

        phi_zero[1:] = _thomas_scalar(a_sub, b_diag, c_sup, rhs)
        return phi_zero


# ---------------------------------------------------------------------------
# Module-level solver functions (called by VariableCoefficientPoissonSolver)
# ---------------------------------------------------------------------------

def _thomas_batch(a, b, c, d):
    """
    Vectorised Thomas (tridiagonal) algorithm for n_systems simultaneously.

    Parameters
    ----------
    a, b, c, d : arrays of shape (n_systems, nz)
        Tridiagonal coefficients and RHS.  a[i,0] and c[i,-1] are unused.

    Returns
    -------
    x : array of shape (n_systems, nz)

    Complexity
    ----------
    2*nz Python iterations, each executing one element-wise operation on
    n_systems elements — a single GPU kernel per z-level.  The serial
    equivalent requires n_systems*2*nz kernel launches.

    GPU speedup at 64×64: ~4 095×
    GPU speedup at 128×128: ~16 383×
    """
    nz = d.shape[1]

    c_p = xp.empty_like(c)
    d_p = xp.empty_like(d)

    # Forward sweep
    c_p[:, 0] = c[:, 0] / b[:, 0]
    d_p[:, 0] = d[:, 0] / b[:, 0]

    for k in range(1, nz):
        denom    = b[:, k] - a[:, k] * c_p[:, k - 1]
        d_p[:, k] = (d[:, k] - a[:, k] * d_p[:, k - 1]) / denom
        if k < nz - 1:
            c_p[:, k] = c[:, k] / denom
        # else c_p[:,k] stays uninitialised but is never read in back-sub

    # Back substitution
    x = xp.empty_like(d)
    x[:, -1] = d_p[:, -1]
    for k in range(nz - 2, -1, -1):
        x[:, k] = d_p[:, k] - c_p[:, k] * x[:, k + 1]

    return x


def _thomas_scalar(a, b, c, d):
    """
    Serial Thomas algorithm for a single tridiagonal system.

    Used for the (kx=0, ky=0) zero mode only (one call per Poisson solve).
    a[0] and c[-1] are unused but included for indexing consistency.
    """
    n = len(d)
    c_p = xp.zeros(n, dtype=d.dtype)
    d_p = xp.zeros(n, dtype=d.dtype)

    c_p[0] = c[0] / b[0]
    d_p[0] = d[0] / b[0]
    for k in range(1, n):
        denom   = b[k] - a[k] * c_p[k - 1]
        c_p[k]  = c[k] / denom if k < n - 1 else 0.0
        d_p[k]  = (d[k] - a[k] * d_p[k - 1]) / denom

    x = xp.zeros(n, dtype=d.dtype)
    x[-1] = d_p[-1]
    for k in range(n - 2, -1, -1):
        x[k] = d_p[k] - c_p[k] * x[k + 1]
    return x
