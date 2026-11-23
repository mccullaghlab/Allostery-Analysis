import numpy as np
import pytest

from allostery.hessian import hessian_3n, laplacian_from_springs, unit_vectors


def chain_system(n=24, seed=0, cutoff=7.5, bond=3.8):
    """A random, compact, protein-like chain of ``n`` nodes with an ENM
    Hessian (spring constants decaying with distance inside a cutoff)."""
    rng = np.random.default_rng(seed)
    radius = 3.0 * n ** (1 / 3) + 3.0
    i, x, tries = 1, np.zeros((n, 3)), 0
    while i < n:
        step = rng.normal(size=3)
        cand = x[i - 1] + step * bond / np.linalg.norm(step)
        # keep the chain compact but avoid overlaps
        if np.all(np.linalg.norm(x[:i] - cand, axis=1) > 3.0) and np.linalg.norm(cand) < radius:
            x[i], i, tries = cand, i + 1, 0
        else:
            tries += 1
            if tries > 500:  # dead end: restart the chain
                i, tries = 1, 0
    _, d = unit_vectors(x)
    k = np.where((d < cutoff) & (d > 0), 5.0 * np.exp(-(d - 3.8) / 2.0), 0.0)
    # stronger backbone springs; i+2 and i+3 springs make the network rigid
    for off, kk in ((1, 20.0), (2, 3.0), (3, 1.5)):
        idx = np.arange(n - off)
        k[idx, idx + off] = k[idx + off, idx] = np.maximum(k[idx, idx + off], kk)
    h = laplacian_from_springs(k)
    e = np.linalg.eigvalsh(hessian_3n(h, x))
    assert e[6] > 1e-8 * e[-1], "synthetic network is not rigid"
    return x, h


@pytest.fixture
def small_system():
    return chain_system(16, seed=1)


@pytest.fixture
def medium_system():
    return chain_system(24, seed=2)
