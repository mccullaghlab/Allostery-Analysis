import numpy as np
import pytest

from allostery import covariance_from_hessian
from allostery.cli import main
from allostery.paths import fortran_available
from allostery.pipeline import Config, run

from .conftest import chain_system

pytestmark = pytest.mark.skipif(not fortran_available(), reason="Fortran kernel unavailable")


@pytest.fixture
def inputs(tmp_path):
    x, h = chain_system(18, seed=12)
    cov = covariance_from_hessian(h, x, 300.0)
    np.savetxt(tmp_path / "cov.dat", cov)
    np.savetxt(tmp_path / "avg.dat", x)
    return tmp_path, x, h


CONFIG = """
output_directory = "out"
temperature = 300.0
index_base = 1

[input]
covariance = "cov.dat"
average_structure = "avg.dat"

[henm]
cutoff = 12.0
alpha = 0.25
max_iterations = 2000
threshold = 1e-6

[paths]
sources = [1, 2]
sinks = [17, 18]
variants = ["direct", "broken_3d"]
n_samples = 2000
seed = 7
n_top = 20
bin_width = 0.1
legacy_output = true
n_jobs = 2

[derivative]
cutoff = 12.0
"""


def test_full_pipeline(inputs):
    d, x, h = inputs
    (d / "config.toml").write_text(CONFIG)
    res = run(d / "config.toml")
    out = d / "out"
    assert res["henm"].converged
    hess = np.loadtxt(out / "hessian.dat")
    np.testing.assert_allclose(hess, h, atol=1e-4)
    for v in ("direct", "broken_3d"):
        vd = out / "paths" / v
        node = np.loadtxt(vd / "node_probability.dat")
        assert node[0, 0] == 1  # 1-based output
        assert node[0, 1] >= 0.0 and np.all(node[:, 1] <= 1.0 + 1e-12)
        assert (vd / "edge_probability.dat").exists()
        assert len((vd / "top_paths.dat").read_text().splitlines()) == 21
        td = np.loadtxt(vd / "thermodynamics.dat")
        assert td.shape == (4, 6) and set(td[:, 0]) == {1, 2}
        assert len(list(vd.glob("paths_*.dat"))) == 4
        assert len(list(vd.glob("summary_*.npz"))) == 4
    nm = np.loadtxt(out / "derivative" / "node_metric.dat")
    assert nm.shape == (18, 4) and np.abs(nm[:, 1]).max() == pytest.approx(1.0)
    assert (out / "run_config.json").exists()

    # re-running only the derivative step reuses hessian.dat
    (out / "derivative" / "node_metric.dat").unlink()
    run(d / "config.toml", steps=["derivative"])
    np.testing.assert_allclose(np.loadtxt(out / "derivative" / "node_metric.dat"), nm)


def test_config_validation(tmp_path):
    with pytest.raises(ValueError, match="unknown option"):
        Config.from_dict({"pathz": {}})
    with pytest.raises(ValueError, match="unknown option"):
        Config.from_dict({"paths": {"n_sample": 3}})
    with pytest.raises(ValueError, match="variant"):
        Config.from_dict({"paths": {"variants": ["3d"]}})


def test_cli(inputs):
    d, x, h = inputs
    o = str(d / "cli")
    assert main(["henm", "--covariance", str(d / "cov.dat"), "--structure", str(d / "avg.dat"),
                 "-T", "300", "--cutoff", "12", "--alpha", "0.25", "--threshold", "1e-6",
                 "-o", o]) == 0
    # Lyman method with an explicit distance-variance target
    from allostery.hessian import projected_variance, unit_vectors

    np.savetxt(d / "dvar.dat", projected_variance(np.loadtxt(d / "cov.dat"), unit_vectors(x)[0]))
    assert main(["henm", "--method", "lyman", "--distance-variance", str(d / "dvar.dat"),
                 "--structure", str(d / "avg.dat"), "-T", "300", "--cutoff", "12",
                 "--max-iterations", "30", "-o", str(d / "cli_lyman")]) == 0
    assert np.loadtxt(d / "cli_lyman" / "hessian.dat").shape == (18, 18)
    with pytest.raises(SystemExit):
        main(["henm", "--structure", str(d / "avg.dat"), "-o", o])
    assert main(["adjacency", "--hessian", o + "/hessian.dat", "-o", o]) == 0
    assert main(["paths", "--adjacency", o + "/adjacency.dat", "--structure", str(d / "avg.dat"),
                 "--sources", "0", "--sinks", "16", "17", "--variant", "direct", "direct_3d",
                 "-n", "500", "--save-ensembles", "-o", o]) == 0
    assert main(["derivative", "--hessian", o + "/hessian.dat", "--structure",
                 str(d / "avg.dat"), "--sources", "0", "--sinks", "17", "-o", o]) == 0
    ens = sorted((d / "cli" / "paths" / "direct").glob("ensemble_*.npz"))
    assert len(ens) == 2
    assert main(["analyze", *map(str, ens), "-o", o + "/reanalysis"]) == 0
    a = np.loadtxt(o + "/reanalysis/node_probability.dat")
    b = np.loadtxt(o + "/paths/direct/node_probability.dat")
    np.testing.assert_allclose(a, b)
    assert main(["info", "--no-compile"]) == 0


def test_pipeline_from_trajectory(tmp_path):
    """covariance -> hENM -> adjacency starting from PDB + DCD files."""
    mda = pytest.importorskip("MDAnalysis")
    from MDAnalysis.coordinates.memory import MemoryReader

    x, h = chain_system(12, seed=3)
    # stiff network: small fluctuations keep the average structure undistorted
    cov = covariance_from_hessian(50 * h, x, 300.0)
    rng = np.random.default_rng(0)
    frames = rng.multivariate_normal(x.ravel(), cov, size=300, method="eigh").reshape(-1, 12, 3)
    u = mda.Universe.empty(12, n_residues=12, atom_resindex=np.arange(12), trajectory=True)
    u.add_TopologyAttr("names", ["CA"] * 12)
    u.add_TopologyAttr("resnames", ["GLY"] * 12)
    u.add_TopologyAttr("resids", np.arange(1, 13))
    u.add_TopologyAttr("masses", [12.0] * 12)
    u.add_TopologyAttr("segids", ["A"])
    u.load_new(frames.astype(np.float32), format=MemoryReader)
    u.atoms.write(tmp_path / "top.pdb")
    with mda.Writer(str(tmp_path / "traj.dcd"), 12) as w:
        for _ in u.trajectory:
            w.write(u.atoms)
    (tmp_path / "c.toml").write_text(
        'output_directory = "out"\ntemperature = 300.0\n'
        '[input]\ntopology = "top.pdb"\ntrajectory = "traj.dcd"\nselection = "all"\n'
        '[henm]\ncutoff = 12.0\nalpha = 0.1\nmax_iterations = 20\n')
    run(tmp_path / "c.toml", steps=["covariance", "henm", "adjacency"])
    out = tmp_path / "out"
    assert np.loadtxt(out / "covariance.dat").shape == (36, 36)
    assert (out / "labels.txt").read_text().split()[0] == "A:GLY1"
    hess = np.loadtxt(out / "hessian.dat")
    assert hess.shape == (12, 12) and np.allclose(hess.sum(1), 0, atol=1e-6)
    avg = np.loadtxt(out / "average_structure.dat")
    d_avg = np.linalg.norm(avg[:-1] - avg[1:], axis=1)
    np.testing.assert_allclose(d_avg, 3.8, atol=0.2)  # bonds of the chain preserved
    adj = np.loadtxt(out / "adjacency.dat")
    assert adj.shape == (12, 12) and adj.max() <= 1

    # pair-distance variances are written and used by the Lyman method
    dvar = np.loadtxt(out / "distance_variance.dat")
    assert dvar.shape == (12, 12) and np.all(dvar[~np.eye(12, dtype=bool)] > 0)
    (tmp_path / "l.toml").write_text(
        'output_directory = "out"\ntemperature = 300.0\n'
        '[henm]\nmethod = "lyman"\ncutoff = 12.0\nmax_iterations = 50\n')
    res = run(tmp_path / "l.toml", steps=["henm"])
    assert res["henm"].method == "lyman"
    np.testing.assert_allclose(res["henm"].target_variance, dvar, rtol=1e-8)
    hist = np.loadtxt(out / "henm_residuals.dat")
    assert hist.shape[1] == 3 and hist[-1, 2] < hist[0, 2]  # relative error decreases
