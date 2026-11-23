"""Containers for sampled path ensembles.

:class:`PathEnsemble`
    every sampled path (node sequence and weights).
:class:`EnsembleSummary`
    the sufficient statistics needed by all analyses in
    :mod:`allostery.paths.analysis` (node/edge visit counts, the minimum-weight
    path and its multiplicity, the per-sample path weights, the ``n_top``
    lowest-weight unique paths).  Summaries of chunks of one simulation can be
    merged exactly, so long simulations never need to hold all paths in
    memory.

Node indices are 0-based.  For *broken* paths every node also carries a
segment label: 0 for the segment attached to the source, 1 for the segment
attached to the sink.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# edge kinds (only used for broken paths)
EDGE_SOURCE_SEGMENT = 0
EDGE_SINK_SEGMENT = 1
EDGE_CROSS = 2


@dataclass(frozen=True)
class TopPath:
    """A unique path and how often it was sampled."""

    weight: float
    count: int
    nodes: tuple[int, ...]
    segments: tuple[int, ...] | None = None

    def signed(self) -> list[int]:
        """1-based, sign-encoded node list as used in the legacy files."""
        if self.segments is None:
            return [n + 1 for n in self.nodes]
        return [(n + 1) * (1 - 2 * s) for n, s in zip(self.nodes, self.segments)]


def _path_key(nodes: np.ndarray, segments: np.ndarray | None) -> bytes:
    if segments is None:
        return nodes.astype(np.int32).tobytes()
    return (nodes.astype(np.int64) * 2 + segments).tobytes()


@dataclass
class EnsembleSummary:
    """Sufficient statistics of a path ensemble (see module docstring)."""

    source: int
    sink: int
    n_nodes: int
    temperature: float
    angular: bool
    broken: bool
    weights: np.ndarray  # (n_samples,) total path weight of each sample
    node_counts: np.ndarray  # (n_nodes,) or (2, n_nodes) for broken paths
    edge_counts: np.ndarray  # (n_kinds, n_nodes, n_nodes) upper triangle, int64
    top_paths: list[TopPath] = field(default_factory=list)
    n_top: int = 1000

    # ------------------------------------------------------------------
    @property
    def n_samples(self) -> int:
        return int(self.weights.size)

    @property
    def min_weight(self) -> float:
        return self.top_paths[0].weight

    @property
    def min_count(self) -> int:
        return self.top_paths[0].count

    @property
    def min_path(self) -> TopPath:
        return self.top_paths[0]

    @property
    def log_partition_function(self) -> float:
        """Estimate of ``ln Z`` with ``Z = sum_paths exp(-w/T)``.

        The minimum-weight path is observed with probability
        ``exp(-w_min/T) / Z``, hence ``Z = N/N_min * exp(-w_min/T)``.
        """
        return float(np.log(self.n_samples) - np.log(self.min_count)
                     - self.min_weight / self.temperature)

    @property
    def mean_weight(self) -> float:
        return float(np.mean(self.weights))

    # ------------------------------------------------------------------
    def merge(self, other: "EnsembleSummary") -> "EnsembleSummary":
        """Exact summary of the concatenation of two samples."""
        for attr in ("source", "sink", "n_nodes", "temperature", "angular", "broken"):
            if getattr(self, attr) != getattr(other, attr):
                raise ValueError(f"cannot merge summaries with different {attr}")
        # Only the min(n_top) best paths of the union have exact counts: a
        # path missing from one summary's list may still have been sampled
        # there (below that summary's cut-off).
        n_top = min(self.n_top, other.n_top)
        return EnsembleSummary(
            source=self.source, sink=self.sink, n_nodes=self.n_nodes,
            temperature=self.temperature, angular=self.angular, broken=self.broken,
            weights=np.concatenate([self.weights, other.weights]),
            node_counts=self.node_counts + other.node_counts,
            edge_counts=self.edge_counts + other.edge_counts,
            top_paths=_merge_top(self.top_paths, other.top_paths, n_top),
            n_top=n_top,
        )

    # ------------------------------------------------------------------
    def save(self, filename) -> None:
        """Save to a compressed ``.npz`` file."""
        lengths = np.array([len(p.nodes) for p in self.top_paths], dtype=np.int64)
        flat = (np.concatenate([np.asarray(p.nodes, dtype=np.int32) for p in self.top_paths])
                if self.top_paths else np.zeros(0, np.int32))
        segs = (np.concatenate([np.asarray(p.segments, dtype=np.int8) for p in self.top_paths])
                if self.broken and self.top_paths else np.zeros(0, np.int8))
        np.savez_compressed(
            filename, kind="EnsembleSummary",
            meta=np.array([self.source, self.sink, self.n_nodes, int(self.angular),
                           int(self.broken), self.n_top], dtype=np.int64),
            temperature=self.temperature, weights=self.weights,
            node_counts=self.node_counts, edge_counts=self.edge_counts,
            top_weights=np.array([p.weight for p in self.top_paths]),
            top_counts=np.array([p.count for p in self.top_paths], dtype=np.int64),
            top_lengths=lengths, top_nodes=flat, top_segments=segs,
        )

    @classmethod
    def load(cls, filename) -> "EnsembleSummary":
        d = np.load(filename)
        source, sink, n_nodes, angular, broken, n_top = (int(v) for v in d["meta"])
        offsets = np.concatenate([[0], np.cumsum(d["top_lengths"])])
        top = []
        for k in range(d["top_weights"].size):
            sl = slice(offsets[k], offsets[k + 1])
            top.append(TopPath(float(d["top_weights"][k]), int(d["top_counts"][k]),
                               tuple(int(v) for v in d["top_nodes"][sl]),
                               tuple(int(v) for v in d["top_segments"][sl]) if broken else None))
        return cls(source=source, sink=sink, n_nodes=n_nodes,
                   temperature=float(d["temperature"]), angular=bool(angular),
                   broken=bool(broken), weights=d["weights"], node_counts=d["node_counts"],
                   edge_counts=d["edge_counts"], top_paths=top, n_top=n_top)


def _merge_top(a: list[TopPath], b: list[TopPath], n_top: int) -> list[TopPath]:
    # A path that belongs to the n_top best paths of the union belongs to the
    # n_top best of every subsample it appears in, so the merge is exact.
    merged: dict[tuple, TopPath] = {}
    for p in list(a) + list(b):
        key = (p.nodes, p.segments)
        if key in merged:
            q = merged[key]
            merged[key] = TopPath(q.weight, q.count + p.count, q.nodes, q.segments)
        else:
            merged[key] = p
    return _truncate_top(sorted(merged.values(), key=lambda p: p.weight), n_top)


def _truncate_top(paths: list[TopPath], n_top: int) -> list[TopPath]:
    """Keep the n_top lowest-weight paths (plus exact ties of the last one)."""
    if len(paths) <= n_top:
        return paths
    cut = paths[n_top - 1].weight
    k = n_top
    while k < len(paths) and paths[k].weight == cut:
        k += 1
    return paths[:k]


@dataclass
class PathEnsemble:
    """All paths sampled between one source and one sink."""

    source: int
    sink: int
    n_nodes: int
    temperature: float
    angular: bool
    broken: bool
    weights: np.ndarray  # (n_samples,) total weight
    edge_weights: np.ndarray  # (n_samples,) sum of edge costs
    angular_weights: np.ndarray  # (n_samples,) angular contribution (0 if not angular)
    offsets: np.ndarray  # (n_samples + 1,) int64
    nodes: np.ndarray  # (offsets[-1],) int32, 0-based
    segments: np.ndarray | None = None  # (offsets[-1],) int8 for broken paths

    # ------------------------------------------------------------------
    def __len__(self) -> int:
        return int(self.weights.size)

    @property
    def lengths(self) -> np.ndarray:
        """Number of nodes in each path."""
        return np.diff(self.offsets)

    def path(self, k: int) -> np.ndarray:
        """Nodes of the ``k``-th sampled path."""
        return self.nodes[self.offsets[k]:self.offsets[k + 1]]

    def path_segments(self, k: int) -> np.ndarray | None:
        if self.segments is None:
            return None
        return self.segments[self.offsets[k]:self.offsets[k + 1]]

    def __iter__(self):
        for k in range(len(self)):
            yield self.path(k)

    # ------------------------------------------------------------------
    @classmethod
    def from_signed(cls, source, sink, n_nodes, temperature, angular, broken,
                    weights3, lengths, signed_nodes) -> "PathEnsemble":
        """Build from 1-based sign-encoded node lists (kernel/legacy format)."""
        weights3 = np.asarray(weights3, dtype=np.float64).reshape(-1, 3)
        signed_nodes = np.asarray(signed_nodes)
        offsets = np.zeros(len(lengths) + 1, dtype=np.int64)
        np.cumsum(lengths, out=offsets[1:])
        nodes = (np.abs(signed_nodes) - 1).astype(np.int32)
        segments = (signed_nodes < 0).astype(np.int8) if broken else None
        return cls(source=source, sink=sink, n_nodes=n_nodes, temperature=temperature,
                   angular=angular, broken=broken, weights=weights3[:, 0].copy(),
                   edge_weights=weights3[:, 1].copy(), angular_weights=weights3[:, 2].copy(),
                   offsets=offsets, nodes=nodes, segments=segments)

    def signed_nodes(self) -> np.ndarray:
        """1-based node indices, negative for sink-segment nodes."""
        s = self.nodes.astype(np.int64) + 1
        if self.segments is not None:
            s = s * (1 - 2 * self.segments.astype(np.int64))
        return s

    @classmethod
    def concatenate(cls, ensembles: list["PathEnsemble"]) -> "PathEnsemble":
        first = ensembles[0]
        lengths = np.concatenate([e.lengths for e in ensembles])
        offsets = np.zeros(lengths.size + 1, dtype=np.int64)
        np.cumsum(lengths, out=offsets[1:])
        return cls(
            source=first.source, sink=first.sink, n_nodes=first.n_nodes,
            temperature=first.temperature, angular=first.angular, broken=first.broken,
            weights=np.concatenate([e.weights for e in ensembles]),
            edge_weights=np.concatenate([e.edge_weights for e in ensembles]),
            angular_weights=np.concatenate([e.angular_weights for e in ensembles]),
            offsets=offsets, nodes=np.concatenate([e.nodes for e in ensembles]),
            segments=(np.concatenate([e.segments for e in ensembles])
                      if first.segments is not None else None),
        )

    # ------------------------------------------------------------------
    def summary(self, n_top: int = 1000) -> EnsembleSummary:
        """Compute the :class:`EnsembleSummary` of this ensemble."""
        n = self.n_nodes
        nodes = self.nodes.astype(np.int64)
        if self.broken:
            node_counts = np.bincount(self.segments.astype(np.int64) * n + nodes,
                                      minlength=2 * n).reshape(2, n)
        else:
            node_counts = np.bincount(nodes, minlength=n)

        # consecutive pairs inside each path
        valid = np.ones(max(nodes.size - 1, 0), dtype=bool)
        starts = self.offsets[1:-1]
        valid[starts[starts > 0] - 1] = False
        a, b = nodes[:-1][valid], nodes[1:][valid]
        n_kinds = 3 if self.broken else 1
        if self.broken:
            sa = self.segments[:-1][valid].astype(np.int64)
            sb = self.segments[1:][valid].astype(np.int64)
            kind = np.where(sa == sb, sa, EDGE_CROSS)
        else:
            kind = np.zeros_like(a)
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        keys = (kind * n + lo) * n + hi
        if n_kinds * n * n <= 50_000_000:
            edge_counts = np.bincount(keys, minlength=n_kinds * n * n).reshape(n_kinds, n, n)
        else:  # pragma: no cover - very large systems
            edge_counts = np.zeros((n_kinds, n, n), dtype=np.int64)
            uk, cnt = np.unique(keys, return_counts=True)
            edge_counts.reshape(-1)[uk] = cnt

        return EnsembleSummary(
            source=self.source, sink=self.sink, n_nodes=n, temperature=self.temperature,
            angular=self.angular, broken=self.broken, weights=self.weights.copy(),
            node_counts=node_counts.astype(np.int64), edge_counts=edge_counts.astype(np.int64),
            top_paths=self.top_paths(n_top), n_top=n_top,
        )

    def top_paths(self, n_top: int = 1000) -> list[TopPath]:
        """The ``n_top`` lowest-weight unique paths with their multiplicities.

        Ties are ordered by first appearance in the simulation.
        """
        if len(self) == 0:
            return []
        w = self.weights
        order = np.argsort(w, kind="stable")
        ws = w[order]
        bounds = np.flatnonzero(np.diff(ws) != 0) + 1
        bounds = np.concatenate([[0], bounds, [ws.size]])
        found: list[TopPath] = []
        for g in range(bounds.size - 1):
            idx = order[bounds[g]:bounds[g + 1]]
            groups: dict[bytes, list] = {}
            for k in idx:  # idx is in order of appearance (stable sort)
                nodes = self.path(k)
                segs = self.path_segments(k)
                key = _path_key(nodes, segs)
                if key in groups:
                    groups[key][0] += 1
                else:
                    groups[key] = [1, nodes, segs]
            for count, nodes, segs in groups.values():
                found.append(TopPath(float(ws[bounds[g]]), count,
                                     tuple(int(v) for v in nodes),
                                     tuple(int(v) for v in segs) if segs is not None else None))
            if len(found) >= n_top:
                break
        return found

    # ------------------------------------------------------------------
    def save(self, filename) -> None:
        """Save to a compressed ``.npz`` file."""
        np.savez_compressed(
            filename, kind="PathEnsemble",
            meta=np.array([self.source, self.sink, self.n_nodes, int(self.angular),
                           int(self.broken)], dtype=np.int64),
            temperature=self.temperature, weights=self.weights,
            edge_weights=self.edge_weights, angular_weights=self.angular_weights,
            offsets=self.offsets, nodes=self.nodes,
            segments=self.segments if self.segments is not None else np.zeros(0, np.int8),
        )

    @classmethod
    def load(cls, filename) -> "PathEnsemble":
        d = np.load(filename)
        source, sink, n_nodes, angular, broken = (int(v) for v in d["meta"])
        return cls(source=source, sink=sink, n_nodes=n_nodes,
                   temperature=float(d["temperature"]), angular=bool(angular),
                   broken=bool(broken), weights=d["weights"], edge_weights=d["edge_weights"],
                   angular_weights=d["angular_weights"], offsets=d["offsets"],
                   nodes=d["nodes"], segments=d["segments"] if broken else None)


def load(filename):
    """Load a :class:`PathEnsemble` or :class:`EnsembleSummary` from ``.npz``."""
    with np.load(filename) as d:
        kind = str(d["kind"])
    if kind == "PathEnsemble":
        return PathEnsemble.load(filename)
    if kind == "EnsembleSummary":
        return EnsembleSummary.load(filename)
    raise ValueError(f"{Path(filename)} does not contain a path ensemble")
