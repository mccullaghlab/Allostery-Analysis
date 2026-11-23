"""Positional covariance of coarse-grained nodes from aligned structures."""

from __future__ import annotations

import numpy as np


def kabsch_rotation(mobile: np.ndarray, reference: np.ndarray,
                    weights: np.ndarray | None = None) -> np.ndarray:
    """Rotation matrix ``R`` minimizing ``|mobile @ R - reference|`` (both centered)."""
    if weights is None:
        h = mobile.T @ reference
    else:
        h = (mobile * weights[:, None]).T @ reference
    v, _, wt = np.linalg.svd(h)
    dsign = np.sign(np.linalg.det(v @ wt))
    v[:, -1] *= dsign if dsign != 0 else 1.0
    return v @ wt


def align(positions: np.ndarray, reference: np.ndarray | None = None,
          weights: np.ndarray | None = None, max_iterations: int = 50,
          tol: float = 1e-8) -> tuple[np.ndarray, np.ndarray]:
    """Superimpose every frame on a reference structure.

    If no ``reference`` is given, frames are iteratively aligned to their own
    average structure (starting from the first frame) until the average
    changes by less than ``tol`` (RMS, A).

    Parameters
    ----------
    positions
        ``(n_frames, N, 3)`` node coordinates.
    weights
        optional per-node weights for the fit (e.g. masses).

    Returns
    -------
    aligned, average
        aligned positions and their average structure.
    """
    x = np.array(positions, dtype=np.float64)
    if x.ndim != 3 or x.shape[2] != 3:
        raise ValueError("positions must have shape (n_frames, N, 3)")
    w = None if weights is None else np.asarray(weights, dtype=np.float64) / np.sum(weights)

    def center(a):
        c = (a.mean(axis=-2, keepdims=True) if w is None
             else np.einsum("...ik,i->...k", a, w)[..., None, :])
        return a - c

    x = center(x)
    fixed = reference is not None
    ref = center(np.asarray(reference, dtype=np.float64)) if fixed else x[0].copy()
    for _ in range(max_iterations if not fixed else 1):
        for f in range(x.shape[0]):
            x[f] = x[f] @ kabsch_rotation(x[f], ref, w)
        avg = x.mean(axis=0)
        if fixed:
            break
        change = np.sqrt(np.mean(np.sum((avg - ref) ** 2, axis=1)))
        ref = avg
        if change < tol:
            break
    return x, x.mean(axis=0)


def covariance(positions: np.ndarray, aligned: bool = False, **align_kwargs):
    """``3N x 3N`` positional covariance and the average structure.

    Parameters
    ----------
    positions
        ``(n_frames, N, 3)`` node coordinates.
    aligned
        set to True if the frames are already superimposed; otherwise they are
        aligned with :func:`align` first.

    Returns
    -------
    cov, average
        ``(3N, 3N)`` covariance in A^2 and ``(N, 3)`` average positions.
    """
    x = np.asarray(positions, dtype=np.float64)
    if not aligned:
        x, _ = align(x, **align_kwargs)
    avg = x.mean(axis=0)
    dx = (x - avg).reshape(x.shape[0], -1)
    cov = dx.T @ dx / x.shape[0]
    return cov, avg


def distance_variance(positions: np.ndarray) -> np.ndarray:
    """Variance of every pair distance over a trajectory.

    ``<(x_ij - <x_ij>)^2>`` with ``x_ij = |r_i - r_j|`` (eq. 1.1 of Lyman et al.,
    Biophys. J. 2008), the target of the ``"lyman"`` hENM method.  Distances
    are invariant to rigid-body motion, so no alignment is needed.

    Parameters
    ----------
    positions
        ``(n_frames, N, 3)`` node coordinates.

    Returns
    -------
    ``(N, N)`` variances in A^2 (zero diagonal).
    """
    x = np.asarray(positions, dtype=np.float64)
    if x.ndim != 3 or x.shape[2] != 3:
        raise ValueError("positions must have shape (n_frames, N, 3)")
    n = x.shape[1]
    # accumulate deviations from the first frame (better conditioned than
    # raw second moments)
    d0 = np.linalg.norm(x[0][:, None] - x[0][None], axis=-1)
    s1 = np.zeros((n, n))
    s2 = np.zeros((n, n))
    for f in x:
        dd = np.linalg.norm(f[:, None] - f[None], axis=-1) - d0
        s1 += dd
        s2 += dd * dd
    m = x.shape[0]
    var = s2 / m - (s1 / m) ** 2
    np.fill_diagonal(var, 0.0)
    return var.clip(min=0.0)


def pearson_adjacency(cov3n: np.ndarray) -> np.ndarray:
    """Absolute normalized Pearson correlation of node displacements.

    ``|tr C_ij| / sqrt(tr C_ii tr C_jj)`` (Fig. 2B of the paper).  The
    diagonal is set to zero (no self-loops).
    """
    from .hessian import block_view

    tr = np.einsum("iaja->ij", block_view(np.asarray(cov3n, dtype=np.float64)))
    dg = np.sqrt(np.diag(tr))
    a = np.abs(tr) / np.outer(dg, dg)
    np.fill_diagonal(a, 0.0)
    return a


def linear_mutual_information(cov3n: np.ndarray) -> np.ndarray:
    """Linear (Gaussian) mutual information between node displacements.

    ``MI_ij = 1/2 (ln det C_ii + ln det C_jj - ln det C_(ij))`` with
    ``C_(ij)`` the 6x6 joint covariance (Lange & Grubmueller, Proteins 2006).
    """
    from .hessian import block_view

    c = block_view(np.asarray(cov3n, dtype=np.float64))
    n = c.shape[0]
    d = np.einsum("iaib->iab", c)
    logdet_i = np.linalg.slogdet(d)[1]
    mi = np.zeros((n, n))
    for i in range(n - 1):
        j = np.arange(i + 1, n)
        joint = np.zeros((j.size, 6, 6))
        joint[:, :3, :3] = d[i]
        joint[:, 3:, 3:] = d[j]
        joint[:, :3, 3:] = c[i, :, j, :]
        joint[:, 3:, :3] = np.swapaxes(c[i, :, j, :], 1, 2)
        mi[i, j] = 0.5 * (logdet_i[i] + logdet_i[j] - np.linalg.slogdet(joint)[1])
    return mi + mi.T


def lmi_adjacency(cov3n: np.ndarray) -> np.ndarray:
    """Normalized linear mutual information ``r_MI = sqrt(1 - exp(-2 MI / 3))``
    (Fig. 2A of the paper); zero diagonal."""
    mi = linear_mutual_information(cov3n)
    a = np.sqrt(1.0 - np.exp(-2.0 * mi.clip(min=0.0) / 3.0))
    np.fill_diagonal(a, 0.0)
    return a
