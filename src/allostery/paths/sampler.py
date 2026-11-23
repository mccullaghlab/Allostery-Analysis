"""Monte Carlo sampling of allosteric paths (Lake, Davidson et al. JCTC 2020).

Paths between a source and a sink node are sampled from the distribution
``P(path) ~ exp(-w(path)/T)`` where the path weight ``w`` is the sum of the
edge costs ``-ln A_ij`` of a normalized adjacency matrix ``A`` (typically the
normalized effective Hessian, see :func:`allostery.adjacency.hessian_adjacency`).

Four variants of the original code are available:

=================  ============  ==========================================
``angular``        ``broken``    original program
=================  ============  ==========================================
False              False         ``paths_TWISP.f``
True               False         ``paths_TWISP_3D.f``
False              True          ``paths_TWISP_broken.f``
True               True          ``paths_TWISP_broken_3D.f``
=================  ============  ==========================================

*angular* paths add ``-sum ln|cos(theta)|`` over consecutive bond angles to the
path weight (this requires node positions).  *Broken* paths consist of two
disconnected segments, one attached to the source and one to the sink, and
describe contributions to the covariance through the null space of the
Hessian.

Example
-------
>>> sampler = PathSampler(adjacency, positions, temperature=1.0, angular=True)
>>> ens = sampler.sample(source=49, sink=303, n_samples=100_000, seed=1)
>>> summaries = sampler.sample_pairs([49, 103], [303, 336], 100_000, seed=1)
"""

from __future__ import annotations

import itertools
import warnings
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterable

import numpy as np

from . import _fortran, _reference
from ._rng import as_state, random_seed
from .ensemble import EnsembleSummary, PathEnsemble

VARIANTS = {
    "direct": dict(angular=False, broken=False),
    "direct_3d": dict(angular=True, broken=False),
    "broken": dict(angular=False, broken=True),
    "broken_3d": dict(angular=True, broken=True),
}


class PathSampler:
    """Sampler of allosteric paths on a weighted graph.

    Parameters
    ----------
    adjacency
        ``(N, N)`` normalized adjacency matrix with entries in ``[0, 1]``.
        Entries below ``1e-10`` are treated as absent edges.
    positions
        ``(N, 3)`` node coordinates; required if ``angular`` is True.
    temperature
        Monte Carlo temperature (in units of the path weight).
    angular
        Include the bond-angle term in the path weight ("3D" paths).
    broken
        Sample broken paths (two segments) instead of connected paths.
    backend
        ``"fortran"``, ``"python"`` or ``"auto"`` (Fortran if available).
    """

    def __init__(self, adjacency, positions=None, temperature: float = 1.0,
                 angular: bool = False, broken: bool = False, backend: str = "auto"):
        A = np.ascontiguousarray(adjacency, dtype=np.float64)
        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ValueError("adjacency must be a square matrix")
        if not np.all(np.isfinite(A)):
            raise ValueError("adjacency contains non-finite values")
        off = A[~np.eye(A.shape[0], dtype=bool)]
        if off.size and (off.max() > 1.0 + 1e-12 or off.min() < 0.0):
            warnings.warn("adjacency entries should lie in [0, 1] (they are interpreted "
                          "as exp(-cost)); entries > 1 give negative edge costs",
                          RuntimeWarning, stacklevel=2)
        self.adjacency = A
        self.n_nodes = A.shape[0]
        if self.n_nodes < 3:
            raise ValueError("the graph needs at least 3 nodes")
        self.angular = bool(angular)
        self.broken = bool(broken)
        if temperature <= 0:
            raise ValueError("temperature must be positive")
        self.temperature = float(temperature)
        if self.angular:
            if positions is None:
                raise ValueError("angular paths require node positions")
            X = np.ascontiguousarray(positions, dtype=np.float64)
            if X.shape != (self.n_nodes, 3):
                raise ValueError(f"positions must have shape ({self.n_nodes}, 3)")
            if np.unique(X, axis=0).shape[0] != self.n_nodes:
                raise ValueError("angular paths require distinct node positions")
            self.positions = X
        else:
            self.positions = (np.ascontiguousarray(positions, dtype=np.float64)
                              if positions is not None else np.zeros((1, 3)))
        if backend not in ("auto", "fortran", "python"):
            raise ValueError("backend must be 'auto', 'fortran' or 'python'")
        if backend == "auto":
            backend = "fortran" if _fortran.load() is not None else "python"
        elif backend == "fortran" and _fortran.load() is None:
            raise RuntimeError("the Fortran kernel is not available")
        self.backend = backend
        self._kernel = _fortran if backend == "fortran" else _reference

    @classmethod
    def from_variant(cls, variant: str, adjacency, positions=None, **kwargs) -> "PathSampler":
        """Construct from a variant name (``direct``, ``direct_3d``, ``broken``,
        ``broken_3d``)."""
        try:
            flags = VARIANTS[variant]
        except KeyError:
            raise ValueError(f"unknown variant {variant!r}; choose from {list(VARIANTS)}")
        return cls(adjacency, positions, **flags, **kwargs)

    # ------------------------------------------------------------------
    def shortest_path(self, source: int, sink: int) -> tuple[np.ndarray, float]:
        """Dijkstra shortest path (0-based nodes) and its weight."""
        self._check_nodes(source, sink)
        path, dist = self._kernel.dijkstra(self.adjacency, source, sink)
        if dist >= 9e99 or path[0] != source + 1:
            raise ValueError(f"node {sink} cannot be reached from node {source}")
        return path - 1, dist

    def initial_path(self, source: int, sink: int) -> np.ndarray:
        """1-based, sign-encoded initial path used by the original programs."""
        if self.broken:
            return np.array([source + 1, -(sink + 1)], dtype=np.int32)
        return (self.shortest_path(source, sink)[0] + 1).astype(np.int32)

    def _check_nodes(self, source, sink):
        for v in (source, sink):
            if not 0 <= v < self.n_nodes:
                raise IndexError(f"node {v} out of range [0, {self.n_nodes})")
        if source == sink:
            raise ValueError("source and sink must differ")

    # ------------------------------------------------------------------
    def sample(self, source: int, sink: int, n_samples: int, seed=None,
               keep_paths: bool = True, n_top: int = 1000, n_equilibration: int = 0,
               initial_path=None, buffer_size: int = 4_000_000,
               callback: Callable[[int], None] | None = None):
        """Sample ``n_samples`` paths between ``source`` and ``sink``.

        One sample is recorded after every Monte Carlo sweep (``N`` attempted
        moves).

        Parameters
        ----------
        source, sink
            0-based node indices.
        n_samples
            Number of recorded paths.
        seed
            Integer seed, ``None`` (random) or a ``ran3`` state array (which
            is advanced in place, allowing to continue a random stream).
        keep_paths
            If True return a :class:`PathEnsemble` with every path, otherwise
            an :class:`EnsembleSummary` (computed chunk by chunk with bounded
            memory).
        n_top
            Number of lowest-weight unique paths kept in summaries.
        n_equilibration
            Number of initial sweeps that are discarded (0 in the original
            code, which starts from the shortest path).
        initial_path
            Optional starting path (0-based nodes, from source to sink; for
            broken paths a pair ``(source_segment, sink_segment)``).
        buffer_size
            Number of node entries per chunk handed to the kernel.
        callback
            Called with the number of samples generated so far.
        """
        self._check_nodes(source, sink)
        n = self.n_nodes
        if isinstance(seed, np.ndarray):
            state = seed
            if state.dtype != np.int32 or state.shape != (57,):
                raise ValueError("a ran3 state must be an int32 array of length 57")
        else:
            state = as_state(seed)

        path = np.zeros(n, dtype=np.int32)
        p0 = (self._encode_initial(initial_path, source, sink) if initial_path is not None
              else self.initial_path(source, sink))
        path[:p0.size] = p0
        npath = int(p0.size)

        buffer_size = max(int(buffer_size), 2 * n)
        if self.backend == "python":
            buffer_size = min(buffer_size, 50 * n)
        nodes_buf = np.zeros(buffer_size, dtype=np.int32)
        max_chunk = max(1, buffer_size // 2)
        lengths_buf = np.zeros(max_chunk, dtype=np.int32)
        weights_buf = np.zeros((max_chunk, 3), dtype=np.float64)

        def run(nsweeps):
            nonlocal npath
            npath, ndone, nnodes, status = self._kernel.sample(
                self.adjacency, self.positions, self.angular, source, sink,
                self.temperature, nsweeps, path, npath, state,
                nodes_buf, lengths_buf, weights_buf)
            if status != 0:  # pragma: no cover - never observed
                warnings.warn("inconsistent insertion index in the path kernel",
                              RuntimeWarning, stacklevel=3)
            return ndone, nnodes

        remaining = int(n_equilibration)
        while remaining > 0:
            ndone, _ = run(min(remaining, max_chunk))
            remaining -= ndone

        chunks: list[PathEnsemble] = []
        summary: EnsembleSummary | None = None
        done = 0
        while done < n_samples:
            ndone, nnodes = run(min(n_samples - done, max_chunk))
            chunk = PathEnsemble.from_signed(
                source, sink, n, self.temperature, self.angular, self.broken,
                weights_buf[:ndone], lengths_buf[:ndone], nodes_buf[:nnodes])
            if keep_paths:
                chunks.append(chunk)
            else:
                s = chunk.summary(n_top)
                summary = s if summary is None else summary.merge(s)
            done += ndone
            if callback is not None:
                callback(done)

        if keep_paths:
            if not chunks:
                return PathEnsemble.from_signed(source, sink, n, self.temperature,
                                                self.angular, self.broken,
                                                np.zeros((0, 3)), [], [])
            return chunks[0] if len(chunks) == 1 else PathEnsemble.concatenate(chunks)
        return summary

    def _encode_initial(self, initial_path, source, sink) -> np.ndarray:
        if self.broken:
            seg_src, seg_snk = initial_path
            seg_src = [int(v) for v in seg_src]
            seg_snk = [int(v) for v in seg_snk]
            if seg_src[0] != source or seg_snk[-1] != sink:
                raise ValueError("broken initial path must start at the source and end at the sink")
            nodes = seg_src + seg_snk
            p = np.array([v + 1 for v in seg_src] + [-(v + 1) for v in seg_snk], dtype=np.int32)
        else:
            nodes = [int(v) for v in initial_path]
            if nodes[0] != source or nodes[-1] != sink:
                raise ValueError("initial path must start at the source and end at the sink")
            p = np.array([v + 1 for v in nodes], dtype=np.int32)
        if len(set(nodes)) != len(nodes):
            raise ValueError("initial path visits a node more than once")
        return p

    # ------------------------------------------------------------------
    def sample_pairs(self, sources: Iterable[int], sinks: Iterable[int], n_samples: int,
                     seed: int | None = None, keep_paths: bool = False, n_top: int = 1000,
                     n_jobs: int = 1, **kwargs) -> dict[tuple[int, int], object]:
        """Sample every (source, sink) combination.

        Each pair gets its own random stream seeded with ``seed + k``.  With
        the Fortran backend, ``n_jobs > 1`` runs pairs in parallel threads.

        Returns a dict mapping ``(source, sink)`` to an
        :class:`EnsembleSummary` (or :class:`PathEnsemble` if
        ``keep_paths``), in the order of ``itertools.product(sources, sinks)``.
        """
        pairs = list(itertools.product(sources, sinks))
        base = random_seed() if seed is None else int(seed)

        def job(k):
            s, t = pairs[k]
            return self.sample(s, t, n_samples, seed=base + k, keep_paths=keep_paths,
                               n_top=n_top, **kwargs)

        if n_jobs != 1 and self.backend == "fortran":
            workers = None if n_jobs in (None, -1, 0) else int(n_jobs)
            with ThreadPoolExecutor(max_workers=workers) as pool:
                results = list(pool.map(job, range(len(pairs))))
        else:
            results = [job(k) for k in range(len(pairs))]
        return dict(zip(pairs, results))
