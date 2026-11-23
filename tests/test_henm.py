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


def test_distance_variance():
    from allostery.covariance import distance_variance

    rng = np.random.default_rng(3)
    frames = rng.normal(size=(500, 6, 3)) + np.arange(6)[None, :, None] * 4
    dv = distance_variance(frames)
    brute = np.var(np.linalg.norm(frames[:, :, None] - frames[:, None], axis=-1), axis=0)
    np.testing.assert_allclose(dv, brute, rtol=1e-9, atol=1e-12)
    # invariant to rigid-body motion of each frame
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    np.testing.assert_allclose(distance_variance(frames @ q + 5.0), dv, rtol=1e-9, atol=1e-12)


def test_lyman_recovers_springs_with_known_topology():
    x, h = chain_system(20, seed=4)
    k_true = spring_constants(h)
    cov = covariance_from_hessian(h, x, 300.0)
    guess = np.where(k_true > 0, 10.0, 0.0)
    from allostery.hessian import laplacian_from_springs

    res = fit_henm(cov, x, 300.0, method="lyman", initial_guess=laplacian_from_springs(guess),
                   threshold=1e-7, max_iterations=3000)
    assert res.converged and res.method == "lyman"
    np.testing.assert_allclose(res.spring_constants, k_true, atol=1e-4)


def test_lyman_all_pairs_reaches_one_percent():
    """Without a known topology, absent springs only decay asymptotically,
    but the 1% criterion of the paper is reached and true springs are fitted."""
    x, h = chain_system(12, seed=4)
    k_true = spring_constants(h)
    cov = covariance_from_hessian(h, x, 300.0)
    res = fit_henm(cov, x, 300.0, method="lyman", max_iterations=20000)
    assert res.converged and res.relative_errors[-1] < 0.01
    k = res.spring_constants
    assert np.all(k >= 0)
    assert np.median(np.abs(k - k_true)[k_true > 1] / k_true[k_true > 1]) < 0.05


def test_lyman_with_distance_variance_target():
    from allostery.covariance import align, distance_variance

    rng = np.random.default_rng(1)
    x, h = chain_system(12, seed=5)
    cov = covariance_from_hessian(50 * h, x, 300.0)  # stiff: ~0.3 A fluctuations
    frames = rng.multivariate_normal(x.ravel(), cov, size=5000, method="eigh").reshape(-1, 12, 3)
    _, avg = align(frames)
    dv = distance_variance(frames)
    res = fit_henm(None, avg, 300.0, method="lyman", target_variance=dv, cutoff=12.0,
                   max_iterations=3000, threshold=0.02)
    assert res.converged
    np.testing.assert_allclose(res.target_variance, dv)
    # the backbone springs (50 * 20) are recovered to within sampling noise
    i = np.arange(11)
    assert np.median(np.abs(res.spring_constants[i, i + 1] / 1000.0 - 1)) < 0.2


def test_henm_method_validation():
    x, h = chain_system(10, seed=1)
    cov = covariance_from_hessian(h, x, 300.0)
    with pytest.raises(ValueError, match="method"):
        fit_henm(cov, x, 300.0, method="voth")
    with pytest.raises(ValueError, match="covariance or target_variance"):
        fit_henm(None, x, 300.0, method="lyman")
