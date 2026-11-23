"""Derivative centrality metrics of the Hessian (section 2.3 of the paper).

For a source ``m`` and a sink ``n`` the edge metric is the derivative of the
squared Frobenius norm of the ``3 x 3`` covariance block ``C_mn`` with respect
to the spring constant ``k_ij`` (eq. 5):

    ``delta_edge^(ij) = d ||C_mn||^2 / d k_ij
                      = -2 kT^2 sum_ab [H^+]_mn,ab (H^+ e_ij)_m,a (H^+ e_ij)_n,b``

where ``e_ij`` is the ``3N`` vector with ``+u_ij`` on node ``i`` and
``-u_ij`` on node ``j`` (the spring contributes ``k_ij e_ij e_ij^T`` to the
Hessian).  Contributions of several source/sink pairs are summed.  The node
metric (eq. 6) is

    ``delta_node^(i) = sum_j k_ij delta_edge^(ij) / max_i' |sum_j k_i'j delta_edge^(i'j)|``.

As a complement (also computed by the original ``edge_effect_multi.py``),
the finite change of ``sum ||C_mn||^2`` upon *removing* spring ``ij`` is
obtained exactly from the Sherman-Morrison formula.

This is a vectorized rewrite of ``legacy/paths/edge_effect_multi.py`` and
``edge_effect_node_multi.py``, with the hard-coded system size, file names,
source/sink lists and ``kT = 0.6`` turned into arguments.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

from .hessian import (
    block_view,
    hessian_3n,
    kT,
    projected_variance,
    pseudo_inverse,
    spring_constants,
    unit_vectors,
)


@dataclass
class DerivativeResult:
    """Result of :func:`derivative_analysis`.

    Attributes
    ----------
    edge
        ``(N, N)`` edge metric ``delta_edge`` summed over source/sink pairs,
        in A^4 / (kcal/mol/A^2) (i.e. A^6 mol/kcal).
    removal
        ``(N, N)`` change of ``sum ||C_mn||^2`` (A^4) if spring ``ij`` is
        removed (only meaningful where a spring exists).
    node_raw
        ``(N,)`` ``sum_j k_ij delta_edge^(ij)`` (A^4).
    node
        ``node_raw`` normalized to ``[-1, 1]`` (eq. 6).
    node_removal
        ``(N,)`` ``sum_j removal_ij``.
    spring_constants
        ``(N, N)`` spring constants used.
    per_pair
        if requested, ``{(m, n): (edge, removal)}`` for every pair.
    """

    edge: np.ndarray
    removal: np.ndarray
    node_raw: np.ndarray
    node: np.ndarray
    node_removal: np.ndarray
    spring_constants: np.ndarray
    sources: list[int] = field(default_factory=list)
    sinks: list[int] = field(default_factory=list)
    per_pair: dict | None = None

    def edge_list(self, only_springs: bool = True) -> np.ndarray:
        """``(M, 5)`` array ``i, j, delta_edge, removal, k_ij`` for ``i < j``."""
        k = self.spring_constants
        mask = np.triu(k > 0 if only_springs else np.ones_like(k, dtype=bool), 1)
        i, j = np.nonzero(mask)
        return np.column_stack([i, j, self.edge[i, j], self.removal[i, j], k[i, j]])


def derivative_analysis(hessian: np.ndarray, positions: np.ndarray, sources, sinks,
                        temperature: float = 303.0, n_zero_modes: int = 6,
                        per_pair: bool = False) -> DerivativeResult:
    """Edge and node derivative metrics for all (source, sink) combinations.

    Parameters
    ----------
    hessian
        ``N x N`` ENM Hessian (see :mod:`allostery.hessian`); clean it first
        with :func:`allostery.hessian.clean_hessian` if needed.
    positions
        ``(N, 3)`` average node positions.
    sources, sinks
        0-based node indices; every combination is used.
    temperature
        Temperature in K (``kT = k_B T``; the original script used
        ``kT = 0.6`` kcal/mol).
    """
    h = np.asarray(hessian, dtype=np.float64)
    n = h.shape[0]
    sources = [int(s) for s in np.atleast_1d(sources)]
    sinks = [int(s) for s in np.atleast_1d(sinks)]
    for v in sources + sinks:
        if not 0 <= v < n:
            raise IndexError(f"node {v} out of range")
    k = spring_constants(h)
    u, _ = unit_vectors(positions)
    kbt2 = kT(temperature) ** 2

    g = pseudo_inverse(hessian_3n(h, positions), n_zero_modes)  # H^+
    gb = block_view(g)  # (N,3,N,3)
    ucu = projected_variance(g, u)  # e_ij^T H^+ e_ij

    def response(node):
        # (H^+ e_ij)_node as a 3-vector for every pair (i, j): (N, N, 3)
        b = gb[node, :, :, :].transpose(1, 0, 2)  # b[i] = block (node, i), (N,3,3)
        bu = np.einsum("iab,ijb->ija", b, u)  # block(node,i) u_ij
        bu_j = np.einsum("jab,ijb->ija", b, u)  # block(node,j) u_ij
        return bu - bu_j

    responses = {v: response(v) for v in set(sources) | set(sinks)}
    edge = np.zeros((n, n))
    removal = np.zeros((n, n))
    pairs = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        for m, t in itertools.product(sources, sinks):
            cmn = gb[m, :, t, :]
            r1, r2 = responses[m], responses[t]
            proj = np.einsum("ab,ija,ijb->ij", cmn, r1, r2, optimize=True)
            d_edge = -2.0 * kbt2 * proj
            # Sherman-Morrison: C_mn after removing k_ij e e^T from H
            factor = k / (1.0 - k * ucu)
            c0 = cmn[None, None] + factor[:, :, None, None] * r1[:, :, :, None] * r2[:, :, None, :]
            d_rem = kbt2 * (np.sum(c0 ** 2, axis=(2, 3)) - np.sum(cmn ** 2))
            np.fill_diagonal(d_edge, 0.0)
            np.fill_diagonal(d_rem, 0.0)
            edge += d_edge
            removal += d_rem
            if per_pair:
                pairs[(m, t)] = (d_edge, d_rem)

    node_raw = np.sum(k * edge, axis=1)
    scale = np.max(np.abs(node_raw))
    node = node_raw / scale if scale > 0 else node_raw.copy()
    return DerivativeResult(edge=edge, removal=removal, node_raw=node_raw, node=node,
                            node_removal=removal.sum(axis=1), spring_constants=k,
                            sources=sources, sinks=sinks, per_pair=pairs if per_pair else None)
