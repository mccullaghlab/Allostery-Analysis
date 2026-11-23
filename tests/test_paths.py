import itertools
import math
from collections import Counter

import numpy as np
import pytest

from allostery import hessian_adjacency
from allostery.paths import (
    EnsembleSummary,
    PathEnsemble,
    PathSampler,
    degeneracy,
    distill,
    ensemble_weights,
    fortran_available,
    legacy,
    load,
    path_probabilities,
    thermodynamics,
)
from allostery.paths._rng import ran3, seed_state

from .conftest import chain_system

needs_fortran = pytest.mark.skipif(not fortran_available(), reason="Fortran kernel unavailable")
VARIANTS = ["direct", "direct_3d", "broken", "broken_3d"]


@pytest.fixture(scope="module")
def graph():
    x, h = chain_system(14, seed=11)
    return x, hessian_adjacency(h)


def tiny_graph(n=7, seed=5):
    rng = np.random.default_rng(seed)
    a = rng.uniform(0.2, 0.9, size=(n, n))
    a = (a + a.T) / 2
    a[rng.uniform(size=(n, n)) < 0.3] = 0
    a = np.minimum(a, a.T)
    np.fill_diagonal(a, 0)
    x = rng.normal(size=(n, 3)) * 4
    return x, a


# ---------------------------------------------------------------------------
def test_ran3_reproducible():
    s1, s2 = seed_state(42), seed_state(42)
    v1 = [ran3(s1) for _ in range(100)]
    v2 = [ran3(s2) for _ in range(100)]
    assert v1 == v2
    assert all(0 < v < 1 for v in v1)
    assert abs(np.mean([ran3(s1) for _ in range(20000)]) - 0.5) < 0.01


@needs_fortran
@pytest.mark.parametrize("variant", VARIANTS)
@pytest.mark.parametrize("temperature", [1.0, 0.7])
def test_fortran_matches_python_reference(graph, variant, temperature):
    x, a = graph
    n = a.shape[0]
    f = PathSampler.from_variant(variant, a, x, temperature=temperature, backend="fortran")
    p = PathSampler.from_variant(variant, a, x, temperature=temperature, backend="python")
    sf, sp = seed_state(7), seed_state(7)
    ef = f.sample(0, n - 1, 150, seed=sf)
    ep = p.sample(0, n - 1, 150, seed=sp)
    np.testing.assert_array_equal(ef.nodes, ep.nodes)
    np.testing.assert_array_equal(ef.offsets, ep.offsets)
    np.testing.assert_array_equal(ef.weights, ep.weights)
    np.testing.assert_array_equal(sf, sp)
    if variant.startswith("broken"):
        np.testing.assert_array_equal(ef.segments, ep.segments)


def _simple_paths(a, start, end, avoid=(), forbid_end=None):
    """All simple paths from start (to end, or any end if end is None)."""
    n = a.shape[0]
    out = []

    def rec(path):
        cur = path[-1]
        if end is not None and cur == end:
            out.append(list(path))
            return
        if end is None:
            out.append(list(path))
        for nb in range(n):
            if a[cur, nb] > 0 and nb not in path and nb not in avoid and nb != forbid_end:
                path.append(nb)
                rec(path)
                path.pop()

    rec([start])
    return out


def _weights(a, x, p, angular):
    w = sum(-math.log(a[i, j]) for i, j in zip(p[:-1], p[1:]))
    if angular:
        for i, j, k in zip(p[:-2], p[1:-1], p[2:]):
            u = (x[j] - x[i]) / np.linalg.norm(x[j] - x[i])
            v = (x[k] - x[j]) / np.linalg.norm(x[k] - x[j])
            w -= math.log(abs(u @ v))
    return w


@needs_fortran
@pytest.mark.parametrize("variant", VARIANTS)
def test_sampler_is_boltzmann(variant):
    """The sampled path distribution is exp(-w/T)/Z (exact enumeration)."""
    x, a = tiny_graph()
    n, src, snk, T = a.shape[0], 0, 6, 1.0
    sampler = PathSampler.from_variant(variant, a, x, temperature=T)
    angular = sampler.angular
    exact = {}
    if sampler.broken:
        for s1 in _simple_paths(a, src, None, forbid_end=snk):
            for s2 in _simple_paths(a, snk, None, avoid=set(s1)):
                key = (tuple(s1), tuple(s2[::-1]))
                exact[key] = _weights(a, x, s1, angular) + _weights(a, x, s2, angular)
    else:
        for p in _simple_paths(a, src, snk):
            exact[(tuple(p),)] = _weights(a, x, p, angular)
    z = sum(math.exp(-w / T) for w in exact.values())
    ens = sampler.sample(src, snk, 60000, seed=3)

    def key(k):
        p = ens.path(k)
        if not sampler.broken:
            return (tuple(p),)
        s = ens.path_segments(k)
        return (tuple(p[s == 0]), tuple(p[s == 1]))

    cnt = Counter(key(k) for k in range(len(ens)))
    assert set(cnt) <= set(exact)
    tv = 0.5 * sum(abs(math.exp(-w / T) / z - cnt.get(k, 0) / len(ens)) for k, w in exact.items())
    assert tv < 0.03, tv
    # sampled path weights agree with the enumeration
    for k in range(0, len(ens), 997):
        assert ens.weights[k] == pytest.approx(exact[key(k)], rel=1e-10, abs=1e-12)
    # partition-function estimate
    s = ens.summary()
    assert math.exp(s.log_partition_function) == pytest.approx(z, rel=0.05)


def test_python_backend_runs(graph):
    x, a = graph
    e = PathSampler(a, x, angular=True, backend="python").sample(0, 13, 20, seed=1)
    assert len(e) == 20 and e.nodes[0] == 0


# ---------------------------------------------------------------------------
@needs_fortran
@pytest.mark.parametrize("variant", ["direct", "broken_3d"])
def test_streaming_summary_is_exact(graph, variant):
    x, a = graph
    s = PathSampler.from_variant(variant, a, x)
    full = s.sample(1, 12, 3000, seed=5, buffer_size=200)  # many small chunks
    streamed = s.sample(1, 12, 3000, seed=5, keep_paths=False, n_top=50, buffer_size=200)
    ref = full.summary(n_top=50)
    np.testing.assert_array_equal(streamed.weights, ref.weights)
    np.testing.assert_array_equal(streamed.node_counts, ref.node_counts)
    np.testing.assert_array_equal(streamed.edge_counts, ref.edge_counts)
    assert [(p.weight, p.count, p.nodes) for p in streamed.top_paths] == \
           [(p.weight, p.count, p.nodes) for p in ref.top_paths]
    assert streamed.min_count == ref.min_count


@needs_fortran
def test_top_paths_counts(graph):
    x, a = graph
    e = PathSampler(a).sample(0, 13, 2000, seed=2)
    top = e.top_paths(10)
    cnt = Counter(tuple(int(v) for v in e.path(k)) for k in range(len(e)))
    for p in top:
        assert cnt[p.nodes] == p.count
    assert [p.weight for p in top] == sorted(p.weight for p in top)
    assert top[0].weight == e.weights.min()


@needs_fortran
def test_save_load_and_legacy_roundtrip(graph, tmp_path):
    x, a = graph
    e = PathSampler(a, x, angular=True, broken=True).sample(0, 13, 300, seed=4)
    e.save(tmp_path / "e.npz")
    e2 = load(tmp_path / "e.npz")
    assert isinstance(e2, PathEnsemble)
    np.testing.assert_array_equal(e.nodes, e2.nodes)
    np.testing.assert_array_equal(e.segments, e2.segments)
    s = e.summary(20)
    s.save(tmp_path / "s.npz")
    s2 = load(tmp_path / "s.npz")
    assert isinstance(s2, EnsembleSummary)
    assert s2.top_paths == s.top_paths
    np.testing.assert_array_equal(s2.edge_counts, s.edge_counts)

    legacy.write_paths(tmp_path / "p.dat", e)
    e3 = legacy.read_paths(tmp_path / "p.dat", n_nodes=14)
    assert e3.broken and e3.angular and (e3.source, e3.sink) == (0, 13)
    np.testing.assert_array_equal(e3.nodes, e.nodes)
    np.testing.assert_array_equal(e3.segments, e.segments)
    np.testing.assert_allclose(e3.weights, e.weights, rtol=1e-11)


@needs_fortran
def test_sample_pairs_parallel_equals_serial(graph):
    x, a = graph
    s = PathSampler(a)
    r1 = s.sample_pairs([0, 1], [12, 13], 500, seed=9, n_jobs=1)
    r2 = s.sample_pairs([0, 1], [12, 13], 500, seed=9, n_jobs=4)
    assert list(r1) == [(0, 12), (0, 13), (1, 12), (1, 13)]
    for k in r1:
        np.testing.assert_array_equal(r1[k].weights, r2[k].weights)


def test_input_validation(graph):
    x, a = graph
    with pytest.raises(ValueError):
        PathSampler(a, angular=True)  # positions missing
    with pytest.raises(ValueError):
        PathSampler(a).sample(3, 3, 10)
    with pytest.raises(IndexError):
        PathSampler(a).sample(0, 99, 10)
    b = a.copy()
    b[13, :] = b[:, 13] = 0  # isolate the sink
    with pytest.raises(ValueError, match="cannot be reached"):
        PathSampler(b).sample(0, 13, 10)
    with pytest.warns(RuntimeWarning):
        PathSampler(a * 3)


# ---------------------------------------------------------------------------
@needs_fortran
def test_analyses(graph):
    x, a = graph
    s = PathSampler(a)
    res = s.sample_pairs([0, 1], [12, 13], 2000, seed=1)
    w = ensemble_weights(res)
    assert w.sum() == pytest.approx(1.0)
    probs = path_probabilities(res)
    # every path contains its own source and sink (and may pass through others)
    assert probs.node[0] >= w[0] + w[1] - 1e-12 and probs.node[1] >= w[2] + w[3] - 1e-12
    assert probs.node[12] >= w[0] + w[2] - 1e-12 and probs.node[13] >= w[1] + w[3] - 1e-12
    single = path_probabilities(res[(0, 12)])
    assert single.node[0] == pytest.approx(1.0) and single.node[12] == pytest.approx(1.0)
    assert np.allclose(probs.edge, probs.edge.T)
    assert np.all(probs.node <= 1 + 1e-12)
    # edge probabilities are consistent with the node probabilities:
    # sum_j P_edge(i,j) = 2 P_node(i) - P(i is an end point)
    ends = np.zeros(14)
    ends[[0, 1]] += [w[0] + w[1], w[2] + w[3]]
    ends[[12, 13]] += [w[0] + w[2], w[1] + w[3]]
    np.testing.assert_allclose(probs.edge.sum(1), 2 * probs.node - ends, atol=1e-12)

    top = distill(res, n_top=25)
    assert len(top) == 25
    assert sum(p.probability for p in distill(res, n_top=10**6)) == pytest.approx(1.0)

    deg = degeneracy(res, bin_width=0.1)
    assert np.sum(deg.probability) * 0.1 == pytest.approx(1.0)
    assert np.all(deg.degeneracy >= 0)

    td = thermodynamics(res)
    assert np.allclose(td.free_energy, -np.log(td.partition_function))
    assert np.allclose(td.entropy, td.internal_energy - td.free_energy)


@needs_fortran
def test_broken_analyses(graph):
    x, a = graph
    res = PathSampler(a, broken=True).sample_pairs([0], [13], 2000, seed=1)
    probs = path_probabilities(res)
    assert probs.node_by_segment.shape == (2, 14)
    assert probs.node_by_segment[0, 0] == pytest.approx(1.0)  # source in source segment
    assert probs.node_by_segment[1, 13] == pytest.approx(1.0)
    # exactly one jump between the segments per path
    assert np.triu(probs.edge_by_kind[2], 1).sum() == pytest.approx(1.0)
