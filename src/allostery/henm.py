"""Heterogeneous elastic network model (hENM) fitting.

Given the ``3N x 3N`` positional covariance ``C`` and average positions of
``N`` coarse-grained nodes (e.g. residue centers of mass from an all-atom MD
simulation), find spring constants ``k_ij`` such that the ENM reproduces the
variances of all pair distances projected on the pair vectors,

    ``V_ij = u_ij^T (C_ii + C_jj - C_ij - C_ji) u_ij``.

The spring constants are refined iteratively (Lyman, Pfaendtner & Voth,
Biophys. J. 2008; the variant described in the SI of Lake, Davidson et al.,
JCTC 2020):

    ``k_ij <- k_ij - alpha (1/V_ij^ENM - 1/V_ij^target)``,

after which negative spring constants and springs beyond the distance cutoff
are set to zero.  ``V^ENM`` is computed from ``C^ENM = kT H^+`` with the six
rigid-body modes removed.

This is a vectorized rewrite of ``legacy/henm/henm_functions.py``, with two
bugs of the original fixed:

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
        Frobenius norm of ``1/V^ENM - 1/V^target`` (over the fitted pairs) at
        every iteration.
    target_variance
        ``N x N`` projected pair-distance variances of the input covariance.
    model_variance
        the same quantity for the ENM at the last evaluated iteration.
    """

    hessian: np.ndarray
    converged: bool
    n_iterations: int
    residuals: list[float] = field(default_factory=list)
    target_variance: np.ndarray | None = None
    model_variance: np.ndarray | None = None

    @property
    def spring_constants(self) -> np.ndarray:
        return spring_constants(self.hessian)


def fit_henm(covariance: np.ndarray, positions: np.ndarray, temperature: float = 303.0,
             cutoff: float | None = None, alpha: float = 1e-2, max_iterations: int = 1000,
             threshold: float = 1e-4, initial_guess=10.0, n_zero_modes: int = 6,
             callback: Callable[[int, np.ndarray, float], None] | None = None,
             verbose: bool = False) -> HENMResult:
    """Fit an hENM Hessian to a positional covariance matrix.

    Parameters
    ----------
    covariance
        ``3N x 3N`` covariance of node positions (A^2), from structures
        aligned to a common reference.
    positions
        ``(N, 3)`` average node positions (A).
    temperature
        Temperature of the simulation (K).
    cutoff
        Springs are only placed between nodes closer than ``cutoff`` (A).
        ``None`` means no cutoff (all pairs).
    alpha
        Step size of the iterative refinement (kcal/mol).
    max_iterations, threshold
        Stop after ``max_iterations`` or when the residual drops below
        ``threshold``.
    initial_guess
        Either a uniform initial spring constant (kcal/(mol A^2); default 10,
        as in the original code) or an ``N x N`` Hessian to restart from.
    n_zero_modes
        Number of zero modes removed in the pseudo-inverse.
    callback
        ``callback(iteration, hessian, residual)`` is called at every
        iteration (e.g. to write checkpoints).
    """
    cov = np.asarray(covariance, dtype=np.float64)
    n = n_nodes_from_3n(cov)
    x = np.asarray(positions, dtype=np.float64)
    if x.shape != (n, 3):
        raise ValueError(f"positions must have shape ({n}, 3) to match the covariance")
    if not np.allclose(cov, cov.T, rtol=1e-6, atol=1e-10 * np.abs(cov).max()):
        raise ValueError("the covariance matrix is not symmetric")
    kbt = kT(temperature)

    u, d = unit_vectors(x)
    uu = u[:, :, :, None] * u[:, :, None, :]  # (N,N,3,3)
    mask = ~np.eye(n, dtype=bool)
    if cutoff is not None:
        mask &= d < cutoff

    target_var = projected_variance(cov, u)
    if np.any(target_var[mask] <= 0):
        raise ValueError("non-positive pair-distance variance in the covariance; "
                         "is it a proper (aligned) positional covariance?")
    target = np.zeros_like(target_var)
    target[mask] = 1.0 / target_var[mask]

    if np.isscalar(initial_guess):
        k = np.full((n, n), float(initial_guess))
    else:
        k = spring_constants(np.asarray(initial_guess, dtype=np.float64))
        if k.shape != (n, n):
            raise ValueError("initial_guess has the wrong shape")
    k = np.where(mask, k.clip(min=0.0), 0.0)

    idx = np.arange(n)
    residuals: list[float] = []
    converged = False
    model_var = None
    k_prev = None
    it = 0
    while it < max_iterations:
        it += 1
        blocks = -k[:, :, None, None] * uu  # off-diagonal 3x3 blocks, -k u u^T
        h3 = blocks.transpose(0, 2, 1, 3).copy()
        h3[idx, :, idx, :] = -blocks.sum(axis=1)
        e, v = np.linalg.eigh(h3.reshape(3 * n, 3 * n))
        if e[n_zero_modes] <= 1e-10 * e[-1]:
            # the update removed too many springs: the network became floppy
            if k_prev is None:
                raise ValueError("the initial elastic network has more than "
                                 f"{n_zero_modes} zero modes; increase the cutoff")
            alpha *= 0.5
            warnings.warn(f"hENM update produced a floppy network; reducing alpha to {alpha:g}",
                          RuntimeWarning, stacklevel=2)
            k = np.where(mask, (k_prev - alpha * diff_prev).clip(min=0.0), 0.0)
            it -= 1
            if alpha < 1e-12:
                raise RuntimeError("hENM fit failed: cannot find a rigid update")
            continue
        cov_model = kbt * (v[:, n_zero_modes:] / e[n_zero_modes:]) @ v[:, n_zero_modes:].T
        model_var = projected_variance(cov_model, u)
        diff = np.zeros_like(model_var)
        ok = mask & (model_var > 0)
        diff[ok] = 1.0 / model_var[ok] - target[ok]
        dev = float(np.linalg.norm(diff))
        residuals.append(dev)
        if verbose:
            print(f"hENM iteration {it:5d}  residual {dev:.6e}")
        if callback is not None:
            callback(it, laplacian_from_springs(k), dev)
        if dev < threshold:
            converged = True
            break
        k_prev, diff_prev = k, diff
        k = np.where(mask, (k - alpha * diff).clip(min=0.0), 0.0)

    return HENMResult(hessian=laplacian_from_springs(k), converged=converged,
                      n_iterations=it, residuals=residuals, target_variance=target_var,
                      model_variance=model_var)
