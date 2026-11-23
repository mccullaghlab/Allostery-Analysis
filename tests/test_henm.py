import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from allostery import covariance_from_hessian, fit_henm
from allostery.hessian import spring_constants

from .conftest import chain_system

LEGACY = Path(__file__).resolve().parents[1] / "legacy" / "henm"


def test_recovers_known_springs():
    x, h = chain_system(20, seed=4)
    cov = covariance_from_hessian(h, x, 300.0)
    res = fit_henm(cov, x, 300.0, cutoff=12.0, alpha=0.25, max_iterations=2000, threshold=1e-7)
    assert res.converged
    np.testing.assert_allclose(res.spring_constants, spring_constants(h), atol=1e-5)
    assert res.residuals[-1] < 1e-7
    np.testing.assert_allclose(res.model_variance, res.target_variance, rtol=1e-5)


def test_recovers_without_cutoff():
    x, h = chain_system(16, seed=5)
    cov = covariance_from_hessian(h, x, 310.0)
    res = fit_henm(cov, x, 310.0, alpha=0.1, max_iterations=3000, threshold=1e-6)
    assert res.converged
    np.testing.assert_allclose(res.spring_constants, spring_constants(h), atol=1e-4)


def test_noisy_covariance_gives_physical_hessian():
    rng = np.random.default_rng(0)
    x, h = chain_system(16, seed=6)
    cov = covariance_from_hessian(h, x, 300.0)
    frames = rng.multivariate_normal(x.ravel(), cov, size=4000, method="eigh")
    dx = frames - frames.mean(0)
    cov_s = dx.T @ dx / len(frames)
    res = fit_henm(cov_s, frames.mean(0).reshape(-1, 3), 300.0, cutoff=8.0, alpha=0.2,
                   max_iterations=400)
    k = res.spring_constants
    assert np.all(k >= 0)
    np.testing.assert_allclose(res.hessian.sum(axis=1), 0, atol=1e-8)
    # the strong backbone springs are recovered to within sampling noise
    i = np.arange(15)
    assert np.median(np.abs(k[i, i + 1] - 20) / 20) < 0.15


def test_floppy_initial_network_raises():
    x, h = chain_system(12, seed=1)
    cov = covariance_from_hessian(h, x, 300.0)
    with pytest.raises(ValueError, match="zero modes"):
        fit_henm(cov, x, 300.0, cutoff=1.0)


def test_shape_validation():
    x, h = chain_system(10, seed=1)
    cov = covariance_from_hessian(h, x, 300.0)
    with pytest.raises(ValueError):
        fit_henm(cov, x[:-1], 300.0)
    with pytest.raises(ValueError):
        fit_henm(cov[:-1, :-1], x, 300.0)


def test_matches_patched_legacy(tmp_path):
    """Same iterates as legacy/henm/henm_functions.py with its two bugs fixed."""
    pytest.importorskip("scipy")
    src = (LEGACY / "henm_functions.py").read_text()
    patched = (src.replace("nNodes = nCartCoords/3", "nNodes = nCartCoords//3")
               .replace("dist_cutoff_mask = (x0 < inverse_distance_cutoff)",
                        "dist_cutoff_mask = (x0 > inverse_distance_cutoff)"))
    assert patched.count("//3") == 1 and patched.count("x0 > inverse") == 1
    (tmp_path / "henm_legacy.py").write_text(patched)

    x, h = chain_system(12, seed=8)
    cov = covariance_from_hessian(h, x, 300.0)
    np.save(tmp_path / "cov.npy", cov)
    np.save(tmp_path / "x.npy", x)
    script = (
        "import numpy as np, henm_legacy\n"
        "h = henm_legacy.perform_henm(np.load('cov.npy'), np.load('x.npy'), './',"
        " max_iterations=25, alpha=0.05, kBT=1.9872e-3*300, threshold=1e-12,"
        " output_step=1000)\n"
        "np.save('h_legacy.npy', h)\n"
    )
    subprocess.run([sys.executable, "-c", script], cwd=tmp_path, check=True,
                   capture_output=True)
    h_legacy = np.load(tmp_path / "h_legacy.npy")
    res = fit_henm(cov, x, 300.0, alpha=0.05, max_iterations=25, threshold=1e-12)
    np.testing.assert_allclose(res.hessian, h_legacy, rtol=1e-8, atol=1e-10)
