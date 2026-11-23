"""Analysis of sampled path ensembles.

These functions replace the Fortran post-processing programs of the original
code:

========================================  =====================================
function                                  original program
========================================  =====================================
:func:`path_probabilities`                ``paths_analyze*_multi.f``
:func:`distill`                           ``paths_distill*_multi.f``
:func:`degeneracy`                        ``paths_deg*_multi.f``
:func:`thermodynamics`                    ``paths_partition*.f``
========================================  =====================================

All functions accept a sequence of :class:`~.ensemble.EnsembleSummary` or
:class:`~.ensemble.PathEnsemble` objects (e.g. one per source/sink pair), or a
dict whose values are such objects.

Combining ensembles
-------------------
Paths between different source/sink pairs are combined with weights equal to
the estimated partition function of each pair,
``Z_k = (N_k / N_min,k) exp(-w_min,k / T)``, where ``w_min,k`` is the weight of
the lowest-weight sampled path and ``N_min,k`` the number of times it was
observed.  Each path of ensemble ``k`` therefore carries the probability
``Z_k / (N_k sum_l Z_l)``.  This is the reweighting of eq. 4 of the paper (the
covariance between the sources and sinks is the sum over all pairs).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .ensemble import EnsembleSummary, PathEnsemble, TopPath, _truncate_top


def _summaries(ensembles, n_top: int = 1000) -> list[EnsembleSummary]:
    if isinstance(ensembles, dict):
        ensembles = list(ensembles.values())
    elif isinstance(ensembles, (EnsembleSummary, PathEnsemble)):
        ensembles = [ensembles]
    out = []
    for e in ensembles:
        out.append(e.summary(n_top) if isinstance(e, PathEnsemble) else e)
    if not out:
        raise ValueError("no ensembles given")
    first = out[0]
    for s in out[1:]:
        if s.n_nodes != first.n_nodes or s.broken != first.broken:
            raise ValueError("all ensembles must have the same number of nodes and path type")
        if s.temperature != first.temperature:
            raise ValueError("all ensembles must have been sampled at the same temperature")
    return out


def ensemble_weights(ensembles) -> np.ndarray:
    """Probability ``Z_k / sum_l Z_l`` of each ensemble (source/sink pair)."""
    s = _summaries(ensembles)
    logz = np.array([e.log_partition_function for e in s])
    w = np.exp(logz - logz.max())
    return w / w.sum()


@dataclass
class PathProbabilities:
    """Node and edge visit probabilities of a (combined) path ensemble.

    Attributes
    ----------
    node
        ``(N,)`` probability ``P_node`` that a node appears in a path.
    edge
        ``(N, N)`` symmetric probability that an edge is traversed.
    node_by_segment
        broken paths only: ``(2, N)`` contributions of the source and the sink
        segment to ``node``.
    edge_by_kind
        broken paths only: ``(3, N, N)`` contributions of edges within the
        source segment, within the sink segment, and of the jump between the
        two segments.
    ensemble_weights
        weight of each input ensemble.
    """

    node: np.ndarray
    edge: np.ndarray
    ensemble_weights: np.ndarray
    node_by_segment: np.ndarray | None = None
    edge_by_kind: np.ndarray | None = None

    def top_nodes(self, k: int = 20) -> list[tuple[int, float]]:
        order = np.argsort(-self.node)[:k]
        return [(int(i), float(self.node[i])) for i in order]

    def edge_list(self, threshold: float = 1e-16) -> np.ndarray:
        """``(M, 3)`` array of ``i, j, P_edge`` for ``i < j`` above threshold."""
        i, j = np.nonzero(np.triu(self.edge, 1) > threshold)
        return np.column_stack([i, j, self.edge[i, j]])


def path_probabilities(ensembles) -> PathProbabilities:
    """Node and edge probabilities, ``P_node`` and ``P_edge`` (Fig. 3B of the paper)."""
    s = _summaries(ensembles)
    w = ensemble_weights(s)
    n = s[0].n_nodes
    node = np.zeros(s[0].node_counts.shape)
    edge = np.zeros(s[0].edge_counts.shape)
    for wk, e in zip(w, s):
        node += wk / e.n_samples * e.node_counts
        edge += wk / e.n_samples * e.edge_counts
    edge = edge + np.swapaxes(edge, -1, -2)  # stored as upper triangle
    if s[0].broken:
        return PathProbabilities(node=node.sum(axis=0), edge=edge.sum(axis=0),
                                 ensemble_weights=w, node_by_segment=node,
                                 edge_by_kind=edge)
    return PathProbabilities(node=node.reshape(n), edge=edge.reshape(n, n), ensemble_weights=w)


@dataclass(frozen=True)
class DistilledPath(TopPath):
    """A unique path with its (re-weighted) probability."""

    source: int = -1
    sink: int = -1
    probability: float = 0.0


def distill(ensembles, n_top: int = 1000) -> list[DistilledPath]:
    """The ``n_top`` lowest-weight unique paths over all ensembles.

    ``count`` is the raw number of times the path was sampled and
    ``probability`` its re-weighted probability in the combined ensemble.
    """
    s = _summaries(ensembles, n_top)
    w = ensemble_weights(s)
    out = []
    for wk, e in zip(w, s):
        for p in e.top_paths:
            out.append(DistilledPath(p.weight, p.count, p.nodes, p.segments, source=e.source,
                                     sink=e.sink, probability=wk * p.count / e.n_samples))
    out.sort(key=lambda p: p.weight)
    return _truncate_top(out, n_top)[:n_top]


@dataclass
class Degeneracy:
    """Histogram of path weights (Fig. 3A of the paper).

    Attributes
    ----------
    weight
        bin centers.
    probability
        re-weighted probability density of the path weight.
    degeneracy
        estimated number of distinct paths in each bin,
        ``g(w) = sum Z_k exp(w/T) / N_k`` over the samples in the bin.
    """

    weight: np.ndarray
    probability: np.ndarray
    degeneracy: np.ndarray
    bin_width: float


def degeneracy(ensembles, bin_width: float = 0.01) -> Degeneracy:
    """Path weight distribution and path degeneracy (``paths_deg*``)."""
    s = _summaries(ensembles)
    w = ensemble_weights(s)
    T = s[0].temperature
    wmin = min(float(e.weights.min()) for e in s)
    wmax = max(float(e.weights.max()) for e in s)
    b0 = int(np.floor(wmin / bin_width))
    nbins = int(np.floor(wmax / bin_width)) - b0 + 1
    prob = np.zeros(nbins)
    deg = np.zeros(nbins)
    for wk, e in zip(w, s):
        b = np.floor(e.weights / bin_width).astype(np.int64) - b0
        prob += np.bincount(b, minlength=nbins) * (wk / e.n_samples / bin_width)
        # Z_k exp(w/T) / N_k = exp((w - w_min,k)/T) / N_min,k
        contrib = np.exp((e.weights - e.min_weight) / T) / e.min_count
        deg += np.bincount(b, weights=contrib, minlength=nbins)
    centers = (np.arange(nbins) + b0 + 0.5) * bin_width
    return Degeneracy(weight=centers, probability=prob, degeneracy=deg, bin_width=bin_width)


@dataclass
class Thermodynamics:
    """Per-ensemble path thermodynamics (``paths_partition*``).

    ``partition_function`` Z, Helmholtz free energy ``A = -T ln Z``,
    internal energy ``U = <w>`` and entropy term ``TS = U - A``; one entry per
    ensemble, in input order.
    """

    sources: np.ndarray
    sinks: np.ndarray
    partition_function: np.ndarray
    free_energy: np.ndarray
    internal_energy: np.ndarray
    entropy: np.ndarray

    def as_table(self) -> np.ndarray:
        return np.column_stack([self.sources, self.sinks, self.partition_function,
                                self.free_energy, self.internal_energy, self.entropy])


def thermodynamics(ensembles) -> Thermodynamics:
    s = _summaries(ensembles)
    logz = np.array([e.log_partition_function for e in s])
    T = s[0].temperature
    A = -T * logz
    U = np.array([e.mean_weight for e in s])
    return Thermodynamics(
        sources=np.array([e.source for e in s]), sinks=np.array([e.sink for e in s]),
        partition_function=np.exp(logz), free_energy=A, internal_energy=U, entropy=U - A)
