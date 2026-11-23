"""Normalized adjacency matrices (graph edge weights in ``[0, 1]``).

Three choices are compared in the paper (Fig. 2): the linear mutual
information ``r_MI``, the Pearson correlation and the effective Hessian.  For
path sampling, edge costs are ``-ln A_ij``.
"""

from __future__ import annotations

import numpy as np

from .covariance import lmi_adjacency, pearson_adjacency

__all__ = ["hessian_adjacency", "pearson_adjacency", "lmi_adjacency", "edge_costs"]


def hessian_adjacency(hessian: np.ndarray) -> np.ndarray:
    """Normalized Hessian adjacency ``|H_ij| / sqrt(H_ii H_jj)`` (Fig. 2C).

    For an ENM Hessian this equals ``k_ij / sqrt(K_i K_j)`` with
    ``K_i = sum_j k_ij``, the trace-based normalization used for the Pearson
    correlation applied to the ``3N x 3N`` Hessian.  Nodes without springs
    get zero adjacency; the diagonal is zero.
    """
    h = np.asarray(hessian, dtype=np.float64)
    if h.ndim != 2 or h.shape[0] != h.shape[1]:
        raise ValueError("hessian must be a square N x N matrix")
    k = -h.copy()
    np.fill_diagonal(k, 0.0)
    k = k.clip(min=0.0)
    deg = np.diag(h).copy()
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(deg > 0, 1.0 / np.sqrt(np.where(deg > 0, deg, 1.0)), 0.0)
    a = k * np.outer(s, s)
    np.fill_diagonal(a, 0.0)
    return a


def edge_costs(adjacency: np.ndarray, threshold: float = 1e-10) -> np.ndarray:
    """Edge costs ``-ln A_ij``; ``inf`` where ``A_ij < threshold`` (no edge)."""
    a = np.asarray(adjacency, dtype=np.float64)
    with np.errstate(divide="ignore"):
        return np.where(a < threshold, np.inf, -np.log(np.where(a < threshold, 1.0, a)))
