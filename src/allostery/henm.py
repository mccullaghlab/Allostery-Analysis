"""Heterogeneous elastic network model (hENM) fitting.

Given the average positions of ``N`` coarse-grained nodes (e.g. residue
centers of mass from an all-atom MD simulation), find spring constants
``k_ij`` such that the elastic network reproduces the fluctuations of all pair
distances observed in the simulation.  The model fluctuations are always the
harmonic (normal-mode) fluctuations ``C^ENM = kT H^+`` projected on the pair
vectors ``u_ij`` of the average structure,

    ``V_ij^ENM = u_ij^T (C_ii + C_jj - C_ij - C_ji) u_ij``,

with the six rigid-body modes removed.  Two variants of the iterative
refinement are available (``method``):

``"projected"`` (default; the algorithm of the original code of Lake,
Davidson et al., JCTC 2020)
    The target is the same projection applied to the MD covariance (a
    linearization of the distance variance), the spring constants are
    updated as

        ``k_ij <- k_ij - alpha (1/V_ij^ENM - 1/V_ij^target)``

    and the fit stops when the Frobenius norm of ``1/V^ENM - 1/V^target``
    drops below ``threshold`` (in A^-2).

``"lyman"`` (Lyman, Pfaendtner & Voth, Biophys. J. 95, 4183 (2008))
    The target is the variance of the scalar distance ``|x_i - x_j|`` over the
    trajectory (eq. 1.1; see :func:`allostery.covariance.distance_variance`),
    the compliances ``1/(4 k_ij)`` are updated as (eq. 1.2)

        ``1/(4 k_ij) <- 1/(4 k_ij) - alpha (V_ij^ENM - V_ij^target)``

    and the fit stops when the relative RMS deviation of ``V^ENM`` from
    ``V^target`` drops below ``threshold`` (the paper used 1%).  Springs are
    driven towards zero asymptotically but never removed.

In both cases negative spring constants and springs beyond the distance
cutoff are set to zero.

The ``"projected"`` method is a vectorized rewrite of
``legacy/henm/henm_functions.py``, with two bugs of the original fixed:

* the distance-cutoff mask was inverted (it kept only springs *beyond* the
  cutoff, so with the default cutoff of 9999 A every spring was set to zero
  after the first iteration);
* ``nNodes = nCartCoords/3`` is a float in Python 3 and made the code crash.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .hessian import (
    kT,
    laplacian_from_springs,
    n_nodes_from_3n,
    projected_variance,
    spring_constants,
    unit_vectors,
)

METHODS = ("projected", "lyman")

#: default step size and convergence threshold of each method
DEFAULT_ALPHA = {"projected": 1e-2}  # "lyman": 1 / (4 kT), see fit_henm
DEFAULT_THRESHOLD = {"projected": 1e-4, "lyman": 1e-2}


@dataclass
class HENMResult:
    """Result of :func:`fit_henm`.

    Attributes
    ----------
    hessian
        ``N x N`` ENM Hessian (``-k_ij`` off the diagonal, row sums on the
        diagonal), in kcal/(mol A^2).
    converged
        whether ``residual < threshold`` was reached.
    n_iterations
        number of iterations performed.
    residuals
        convergence measure at every iteration: the Frobenius norm of
        ``1/V^ENM - 1/V^target`` (``"projected"``) or the relative RMS
        deviation of ``V^ENM`` from ``V^target`` (``"lyman"``), over the fitted
        pairs.
    relative_errors
        relative RMS deviation of ``V^ENM`` from ``V^target`` at every
        iteration (for both methods, to compare them).
    target_variance
        ``N x N`` target pair-distance variances.
    model_variance
        the projected ENM variances at the last evaluated iteration.
    method
        the method used.
    """

    hessian: np.ndarray
    converged: bool
    n_iterations: int
    residuals: list[float] = field(default_factory=list)
    relative_errors: list[float] = field(default_factory=list)
    target_variance: np.ndarray | None = None
    model_variance: np.ndarray | None = None
    method: str = "projected"

    @property
    def spring_constants(self) -> np.ndarray:
        return spring_constants(self.hessian)


def fit_henm(covariance: np.ndarray | None, positions: np.ndarray, temperature: float = 303.0,
             cutoff: float | None = None, alpha: float | None = None,
             max_iterations: int = 1000, threshold: float | None = None, initial_guess=10.0,
             n_zero_modes: int = 6, method: str = "projected",
             target_variance: np.ndarray | None = None, max_spring: float = 1e4,
             callback: Callable[[int, np.ndarray, float], None] | None = None,
             verbose: bool = False) -> HENMResult:
    """Fit an hENM Hessian to MD fluctuations.

    Parameters
    ----------
    covariance
        ``3N x 3N`` covariance of node positions (A^2), from structures
        aligned to a common reference.  Used to build the target unless
        ``target_variance`` is given (then it may be ``None``).
    positions
        ``(N, 3)`` average node positions (A).
    temperature
        Temperature of the simulation (K).
    cutoff
        Springs are only placed between nodes closer than ``cutoff`` (A).
        ``None`` means no cutoff (all pairs), as in Lyman et al.
    alpha
        Step size.  Default: 0.01 kcal/mol for ``"projected"``; ``1/(4 kT)``
        for ``"lyman"``, for which the update of an isolated spring is exact.
    max_iterations, threshold
        Stop after ``max_iterations`` or when the residual drops below
        ``threshold`` (default 1e-4 A^-2 for ``"projected"``, 0.01 relative
        RMS for ``"lyman"``).
    initial_guess
        Either a uniform initial spring constant (kcal/(mol A^2); default 10,
        as in the original code) or an ``N x N`` Hessian to restart from.
    n_zero_modes
        Number of zero modes removed in the pseudo-inverse.
    method
        ``"projected"`` or ``"lyman"`` (see the module documentation).
    target_variance
        ``N x N`` target pair-distance variances (A^2), e.g. from
        :func:`allostery.covariance.distance_variance`.  This is the target of
        Lyman et al.; without it the projected variance of ``covariance`` is
        used (for ``"lyman"`` that is a linearized approximation).
    max_spring
        Upper bound on spring constants for ``"lyman"`` (the compliance update
        can otherwise overshoot to negative compliances).
    callback
        ``callback(iteration, hessian, residual)`` is called at every
        iteration (e.g. to write checkpoints).
    """
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}")
    x = np.asarray(positions, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    n = x.shape[0]
    kbt = kT(temperature)
    if alpha is None:
        alpha = DEFAULT_ALPHA["projected"] if method == "projected" else 1.0 / (4.0 * kbt)
    if threshold is None:
        threshold = DEFAULT_THRESHOLD[method]

    u, d = unit_vectors(x)
    uu = u[:, :, :, None] * u[:, :, None, :]  # (N,N,3,3)
    mask = ~np.eye(n, dtype=bool)
    if cutoff is not None:
        mask &= d < cutoff

    if target_variance is not None:
        target_var = np.asarray(target_variance, dtype=np.float64)
        if target_var.shape != (n, n):
            raise ValueError(f"target_variance must have shape ({n}, {n})")
    else:
        if covariance is None:
            raise ValueError("either covariance or target_variance must be given")
        cov = np.asarray(covariance, dtype=np.float64)
        if n_nodes_from_3n(cov) != n:
            raise ValueError(f"positions must have shape ({n_nodes_from_3n(cov)}, 3) "
                             "to match the covariance")
        if not np.allclose(cov, cov.T, rtol=1e-6, atol=1e-10 * np.abs(cov).max()):
            raise ValueError("the covariance matrix is not symmetric")
        target_var = projected_variance(cov, u)
    if np.any(target_var[mask] <= 0):
        raise ValueError("non-positive pair-distance variance in the target; "
                         "is it a proper (aligned) positional covariance?")
    target_inv = np.zeros_like(target_var)
    target_inv[mask] = 1.0 / target_var[mask]

    if np.isscalar(initial_guess):
        k = np.full((n, n), float(initial_guess))
    else:
        k = spring_constants(np.asarray(initial_guess, dtype=np.float64))
        if k.shape != (n, n):
            raise ValueError("initial_guess has the wrong shape")
    k = np.where(mask, k.clip(min=0.0), 0.0)

    def update(k, model_var, step):
        if method == "projected":
            diff = np.where(mask, 1.0 / np.where(mask, model_var, 1.0) - target_inv, 0.0)
            return np.where(mask, (k - step * diff).clip(min=0.0), 0.0)
        with np.errstate(divide="ignore"):
            compliance = np.where(k > 0, 1.0 / (4.0 * np.where(k > 0, k, 1.0)), np.inf)
        compliance = compliance - step * (model_var - target_var)
        compliance = np.maximum(compliance, 1.0 / (4.0 * max_spring))
        return np.where(mask & np.isfinite(compliance), 1.0 / (4.0 * compliance), 0.0)

    idx = np.arange(n)
    residuals: list[float] = []
    rel_errors: list[float] = []
    converged = False
    model_var = None
    prev = None  # (k, model variance) of the last accepted iteration
    it = 0
    while it < max_iterations:
        it += 1
        blocks = -k[:, :, None, None] * uu  # off-diagonal 3x3 blocks, -k u u^T
        h3 = blocks.transpose(0, 2, 1, 3).copy()
        h3[idx, :, idx, :] = -blocks.sum(axis=1)
        e, v = np.linalg.eigh(h3.reshape(3 * n, 3 * n))
        if e[n_zero_modes] <= 1e-10 * e[-1]:
            # the update removed too many springs: the network became floppy
            if prev is None:
                raise ValueError("the initial elastic network has more than "
                                 f"{n_zero_modes} zero modes; increase the cutoff")
            alpha *= 0.5
            warnings.warn(f"hENM update produced a floppy network; reducing alpha to {alpha:g}",
                          RuntimeWarning, stacklevel=2)
            if alpha < 1e-12:
                raise RuntimeError("hENM fit failed: cannot find a rigid update")
            k = update(prev[0], prev[1], alpha)
            it -= 1
            continue
        cov_model = kbt * (v[:, n_zero_modes:] / e[n_zero_modes:]) @ v[:, n_zero_modes:].T
        model_var = projected_variance(cov_model, u)
        if np.any(model_var[mask] <= 0):
            raise RuntimeError("non-positive ENM pair variance (numerically singular network)")
        rel = float(np.sqrt(np.mean(((model_var - target_var)[mask] / target_var[mask]) ** 2)))
        if method == "projected":
            dev = float(np.linalg.norm(np.where(mask, 1.0 / np.where(mask, model_var, 1.0)
                                                - target_inv, 0.0)))
        else:
            dev = rel
        residuals.append(dev)
        rel_errors.append(rel)
        if verbose:
            print(f"hENM iteration {it:5d}  residual {dev:.6e}  relative error {rel:.3e}")
        if callback is not None:
            callback(it, laplacian_from_springs(k), dev)
        if dev < threshold:
            converged = True
            break
        prev = (k, model_var)
        k = update(k, model_var, alpha)

    return HENMResult(hessian=laplacian_from_springs(k), converged=converged,
                      n_iterations=it, residuals=residuals, relative_errors=rel_errors,
                      target_variance=target_var, model_variance=model_var, method=method)
