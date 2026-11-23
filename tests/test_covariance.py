import numpy as np
import pytest

from allostery import (
    align,
    compute_covariance,
    covariance_from_hessian,
    hessian_3n,
    hessian_adjacency,
    lmi_adjacency,
    pearson_adjacency,
)
from allostery.adjacency import edge_costs
from allostery.covariance import linear_mutual_information
from allostery.hessian import pseudo_inverse, unit_vectors

from .conftest import chain_system


def random_rotation(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def test_hessian_3n_properties():
    x, h = chain_system(12, seed=3)
    h3 = hessian_3n(h, x)
    assert np.allclose(h3, h3.T)
    e = np.linalg.eigvalsh(h3)
    assert np.all(np.abs(e[:6]) < 1e-9 * e[-1]) and e[6] > 1e-6
    # same construction as the original code (einsum + block_diag)
    u, _ = unit_vectors(x)
    uu = np.einsum("ijk,ijl->ijkl", u, u)
    legacy = np.einsum("ij,ijkl->ikjl", h, uu)
    from scipy.linalg import block_diag

    legacy = legacy.reshape(36, 36) - block_diag(*np.sum(np.einsum("ij,ijkl->ikjl", h, uu), axis=2))
    np.testing.assert_allclose(h3, legacy, atol=1e-12)


def test_pseudo_inverse_tolerance_mode():
    x, h = chain_system(10, seed=3)
    h3 = hessian_3n(h, x)
    np.testing.assert_allclose(pseudo_inverse(h3), pseudo_inverse(h3, n_zero_modes=None),
                               rtol=1e-6, atol=1e-10)


def test_align_and_covariance():
    rng = np.random.default_rng(1)
    x, h = chain_system(12, seed=2)
    cov = covariance_from_hessian(h, x, 300.0)
    frames = rng.multivariate_normal(x.ravel(), cov, size=3000, method="eigh").reshape(-1, 12, 3)
    moved = np.array([f @ random_rotation(rng) + rng.normal(size=3) * 5 for f in frames])
    aligned, avg = align(moved)
    # the average structure is recovered up to a global rotation/translation
    r_avg, _ = align(np.array([avg]), reference=x)
    assert np.sqrt(np.mean(np.sum((r_avg[0] - (x - x.mean(0))) ** 2, axis=1))) < 0.05
    c, a = compute_covariance(moved)
    # rotation-invariant comparison: projected pair-distance variances
    from allostery.hessian import projected_variance

    u, d = unit_vectors(a)
    v_est = projected_variance(c, u)
    v_ref = projected_variance(cov, unit_vectors(x)[0])
    near = (d > 0) & (d < 8)
    assert np.median(np.abs(v_est[near] / v_ref[near] - 1)) < 0.1


def test_adjacencies():
    x, h = chain_system(12, seed=2)
    a = hessian_adjacency(h)
    assert np.allclose(a, a.T) and np.all(a >= 0) and np.all(a <= 1)
    assert np.all(np.diag(a) == 0)
    k = -h.copy()
    np.fill_diagonal(k, 0)
    assert a[0, 1] == pytest.approx(k[0, 1] / np.sqrt(h[0, 0] * h[1, 1]))
    c = covariance_from_hessian(h, x, 300.0)
    p = pearson_adjacency(c)
    assert np.allclose(p, p.T) and np.all(p <= 1 + 1e-12)
    m = lmi_adjacency(c)
    assert np.allclose(m, m.T) and np.all((m >= 0) & (m < 1))
    mi = linear_mutual_information(c)
    # MI of a 1x1-block Gaussian pair equals -1/2 log(1 - rho^2) when isotropic;
    # here only check it is positive off the diagonal and zero on it
    assert np.all(mi[~np.eye(12, dtype=bool)] > 0)
    costs = edge_costs(a)
    assert np.isinf(costs[0, 0]) and np.all(costs[a > 0] >= 0)


def test_trajectory_node_positions(tmp_path):
    mda = pytest.importorskip("MDAnalysis")
    from MDAnalysis.coordinates.memory import MemoryReader

    from allostery.trajectory import node_positions, trajectory_covariance

    n_res, per_res = 6, 3
    u = mda.Universe.empty(n_res * per_res, n_residues=n_res,
                           atom_resindex=np.repeat(np.arange(n_res), per_res),
                           trajectory=True)
    u.add_TopologyAttr("names", ["N", "CA", "C"] * n_res)
    u.add_TopologyAttr("resnames", ["ALA"] * n_res)
    u.add_TopologyAttr("resids", np.arange(1, n_res + 1))
    u.add_TopologyAttr("masses", [14.0, 12.0, 12.0] * n_res)
    u.add_TopologyAttr("segids", ["A"])
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(20, n_res * per_res, 3)) + np.arange(n_res * per_res)[None, :, None]
    u.load_new(coords.astype(np.float32), format=MemoryReader)
    pos, labels = node_positions(u, selection="all", node="com")
    assert pos.shape == (20, n_res, 3) and labels[0] == "A:ALA1"
    w = np.array([14.0, 12.0, 12.0])
    expected = np.einsum("fia,i->fa", coords[:, :3].astype(np.float32), w) / w.sum()
    np.testing.assert_allclose(pos[:, 0], expected, rtol=1e-5)
    ca, _ = node_positions(u, selection="all", node="ca")
    np.testing.assert_allclose(ca[:, 1], coords[:, 4], rtol=1e-5)
    cov, avg, _ = trajectory_covariance(u, selection="all")
    assert cov.shape == (18, 18) and avg.shape == (6, 3)
