import numpy as np
import pytest

from allostery import clean_hessian, derivative_analysis
from allostery.hessian import (
    block_view,
    hessian_3n,
    kT,
    laplacian_from_springs,
    pseudo_inverse,
    spring_constants,
)

from .conftest import chain_system

T = 300.0


def frob2(h, x, pairs):
    c = block_view(kT(T) * pseudo_inverse(hessian_3n(h, x)))
    return sum(np.sum(c[m, :, t, :] ** 2) for m, t in pairs)


@pytest.fixture(scope="module")
def system():
    return chain_system(14, seed=9)


def test_edge_metric_is_the_derivative(system):
    x, h = system
    src, snk = [0, 1], [12, 13]
    pairs = [(m, t) for m in src for t in snk]
    res = derivative_analysis(h, x, src, snk, temperature=T)
    k = spring_constants(h)
    scale = np.abs(res.edge).max()
    rng = np.random.default_rng(0)
    edges = [(i, j) for i in range(14) for j in range(i + 1, 14)]
    for idx in rng.choice(len(edges), 12, replace=False):
        i, j = edges[idx]
        eps = 1e-4
        kp, km = k.copy(), k.copy()
        kp[i, j] = kp[j, i] = k[i, j] + eps
        km[i, j] = km[j, i] = k[i, j] - eps
        fd = (frob2(laplacian_from_springs(kp), x, pairs)
              - frob2(laplacian_from_springs(km), x, pairs)) / (2 * eps)
        assert res.edge[i, j] == pytest.approx(fd, rel=1e-5, abs=1e-8 * scale)
        assert res.edge[j, i] == res.edge[i, j]


def test_removal_metric(system):
    x, h = system
    src, snk = [0], [13]
    res = derivative_analysis(h, x, src, snk, temperature=T)
    k = spring_constants(h)
    base = frob2(h, x, [(0, 13)])
    # a spring whose removal keeps the network rigid
    i, j = 3, 4
    k2 = k.copy()
    k2[i, j] = k2[j, i] = 0.0
    exact = frob2(laplacian_from_springs(k2), x, [(0, 13)]) - base
    assert res.removal[i, j] == pytest.approx(exact, rel=1e-6)


def test_node_metric(system):
    x, h = system
    res = derivative_analysis(h, x, [0, 1], [12, 13], temperature=T, per_pair=True)
    assert np.max(np.abs(res.node)) == pytest.approx(1.0)
    np.testing.assert_allclose(res.node_raw, np.sum(res.spring_constants * res.edge, axis=1))
    total = sum(e for e, _ in res.per_pair.values())
    np.testing.assert_allclose(total, res.edge)
    assert res.edge_list().shape[1] == 5


def test_clean_hessian(system):
    x, h = system
    bad = h.copy()
    bad[0, 9] = bad[9, 0] = 2.0  # negative spring
    c = clean_hessian(bad, x, cutoff=6.0)
    k = spring_constants(c)
    assert k[0, 9] == 0 and np.all(k >= 0)
    d = np.linalg.norm(x[:, None] - x[None], axis=-1)
    assert np.all(k[d > 6.0] == 0)
    np.testing.assert_allclose(c.sum(axis=1), 0, atol=1e-12)
