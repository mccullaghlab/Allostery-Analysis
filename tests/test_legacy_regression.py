"""Regression tests against the original Fortran programs in ``legacy/paths``.

The original programs are compiled with the hard-coded system size
(``npar=454``) replaced by the size of a small synthetic system, run with
the same ``ran3`` state as the new code, and every output is compared.
"""

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from allostery import hessian_adjacency
from allostery.paths import (
    PathSampler,
    fortran_available,
    degeneracy,
    distill,
    legacy,
    path_probabilities,
    thermodynamics,
)
from allostery.paths._rng import read_ran3_file, seed_state, write_ran3_file

from .conftest import chain_system

LEGACY = Path(__file__).resolve().parents[1] / "legacy" / "paths"
GFORTRAN = shutil.which("gfortran")

pytestmark = [
    pytest.mark.legacy,
    pytest.mark.skipif(GFORTRAN is None, reason="gfortran is needed to build the legacy code"),
    pytest.mark.skipif(not fortran_available(), reason="Fortran kernel unavailable"),
]

N = 24
SOURCES = [0, 1, 2, 3]
SINKS = [20, 21, 22, 23]
NMAX = 400

VARIANTS = {
    # name: (TWISP program, suffix of the analysis programs, angular, broken)
    "direct": ("paths_TWISP", "", False, False),
    "direct_3d": ("paths_TWISP_3D", "_3D", True, False),
    "broken": ("paths_TWISP_broken", "_broken", False, True),
    "broken_3d": ("paths_TWISP_broken_3D", "_broken_3D", True, True),
}


@pytest.fixture(scope="module")
def legacy_bin(tmp_path_factory):
    """Compile every original program with npar = N."""
    out = tmp_path_factory.mktemp("legacy_bin")
    for src in LEGACY.glob("*.f"):
        code = src.read_text().replace("npar=454", f"npar={N}")
        f = out / src.name
        f.write_text(code)
        exe = out / (src.stem + ".e")
        subprocess.run([GFORTRAN, "-O2", "-ffp-contract=off", "-ffixed-line-length-none",
                        str(f), "-o", str(exe)], check=True, capture_output=True)
    return out


@pytest.fixture(scope="module")
def system():
    x, h = chain_system(N, seed=7)
    return x, hessian_adjacency(h)


def run(exe, workdir, input_lines):
    (workdir / "input.dat").write_text("\n".join(str(v) for v in input_lines) + "\n")
    subprocess.run([str(exe)], cwd=workdir, check=True, capture_output=True)


def legacy_paths(legacy_bin, workdir, system, variant, temperature, seed):
    """Run the original sampler for all 16 pairs; return the path files."""
    x, a = system
    prog = VARIANTS[variant][0]
    legacy.write_adjacency(workdir / "adj.dat", a)
    np.savetxt(workdir / "xyz.dat", x)
    files = {}
    for k, (s, t) in enumerate((s, t) for s in SOURCES for t in SINKS):
        write_ran3_file(workdir / "ran3.dat", seed_state(seed + k))
        fout = f"paths_{variant}.{s + 1:03d}.{t + 1:03d}.dat"
        run(legacy_bin / f"{prog}.e", workdir,
            ["adj.dat", "xyz.dat", fout, s + 1, t + 1, NMAX, f"{temperature:.2f}"])
        files[(s, t)] = (workdir / fout, read_ran3_file(workdir / "ran3.dat"))
    return files


@pytest.mark.parametrize("variant", list(VARIANTS))
@pytest.mark.parametrize("temperature", [1.0, 0.8])
def test_sampler_and_analyses_match_legacy(legacy_bin, tmp_path, system, variant, temperature):
    x, a = system
    _, suffix, angular, broken = VARIANTS[variant]
    seed = 1234
    files = legacy_paths(legacy_bin, tmp_path, system, variant, temperature, seed)

    # ---------------------------------------------------------------- sampling
    # the adjacency is passed through the legacy text format (as the original
    # programs read it) so that both codes see identical numbers
    a_txt = np.loadtxt(tmp_path / "adj.dat")
    sampler = PathSampler(a_txt, x, temperature=temperature, angular=angular, broken=broken,
                          backend="fortran")
    new = {}
    for k, ((s, t), (fname, final_state)) in enumerate(files.items()):
        state = seed_state(seed + k)
        ens = sampler.sample(s, t, NMAX, seed=state)
        old = legacy.read_paths(fname, s, t, N, temperature, angular, broken)
        assert len(old) == NMAX
        np.testing.assert_array_equal(ens.offsets, old.offsets)
        np.testing.assert_array_equal(ens.nodes, old.nodes)
        if broken:
            np.testing.assert_array_equal(ens.segments, old.segments)
        np.testing.assert_allclose(ens.weights, old.weights, rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(ens.edge_weights, old.edge_weights, rtol=1e-11, atol=1e-11)
        np.testing.assert_array_equal(state, final_state)  # same random stream consumed
        new[(s, t)] = old  # analyses below use the legacy (rounded) weights

    path_list = [str(files[p][0].name) for p in files]

    # ------------------------------------------------- node/edge probabilities
    run(legacy_bin / f"paths_analyze{suffix}_multi.e", tmp_path,
        ["node.dat", "edge.dat", NMAX, f"{temperature:.2f}", len(path_list), *path_list])
    probs = path_probabilities(new)
    node_old = np.loadtxt(tmp_path / "node.dat")
    edge_old = np.atleast_2d(np.loadtxt(tmp_path / "edge.dat"))
    if broken:
        np.testing.assert_allclose(probs.node_by_segment.T, node_old[:, 1:], rtol=1e-12, atol=1e-15)
        i, j = edge_old[:, 0].astype(int) - 1, edge_old[:, 1].astype(int) - 1
        for kind in range(3):
            np.testing.assert_allclose(probs.edge_by_kind[kind][i, j], edge_old[:, 2 + kind],
                                       rtol=1e-12, atol=1e-15)
        n_edges = np.count_nonzero(np.triu(probs.edge_by_kind.sum(0), 1) > 1e-16)
    else:
        np.testing.assert_allclose(probs.node, node_old[:, 1], rtol=1e-12, atol=1e-15)
        i, j = edge_old[:, 0].astype(int) - 1, edge_old[:, 1].astype(int) - 1
        np.testing.assert_allclose(probs.edge[i, j], edge_old[:, 2], rtol=1e-12)
        n_edges = np.count_nonzero(np.triu(probs.edge, 1) > 1e-16)
    assert n_edges == len(edge_old)

    # -------------------------------------------------------------- distilled
    run(legacy_bin / f"paths_distill{suffix}_multi.e", tmp_path,
        ["distill.dat", NMAX, len(path_list), *path_list])
    old_rows = []
    for line in (tmp_path / "distill.dat").read_text().splitlines():
        tok = line.split()
        if float(tok[0].replace("+100", "E+100") if "E" not in tok[0] else tok[0]) > 1e98:
            break
        old_rows.append((float(tok[0]), int(tok[1]), [int(v) for v in tok[3:3 + int(tok[2])]]))
    top = distill(new, n_top=1000)
    assert len(top) == len(old_rows)
    # equal weights may be ordered differently; compare as multisets per weight
    as_set = lambda rows: sorted((w, c, tuple(p)) for w, c, p in rows)  # noqa: E731
    assert as_set(old_rows) == as_set([(p.weight, p.count, p.signed()) for p in top])

    # -------------------------------------------------------------- degeneracy
    dl = 0.05
    run(legacy_bin / f"paths_deg{suffix}_multi.e", tmp_path,
        ["deg.dat", NMAX, f"{temperature:.2f}", dl, len(path_list), *path_list])
    old = np.loadtxt(tmp_path / "deg.dat")
    deg = degeneracy(new, bin_width=dl)
    # the original prints bin b = int(w/dl) at (b - 0.5) dl instead of (b + 0.5) dl
    old_centers = old[:, 0] + dl
    for c, p, g in zip(deg.weight, deg.probability, deg.degeneracy):
        r = np.argmin(np.abs(old_centers - c))
        assert abs(old_centers[r] - c) < 1e-6
        np.testing.assert_allclose([p, g], old[r, 1:], rtol=1e-10, atol=1e-14)
    assert np.isclose(old[:, 1].sum(), deg.probability.sum())

    # ---------------------------------------------------------- thermodynamics
    run(legacy_bin / f"paths_partition{suffix}.e", tmp_path,
        ["part.dat", NMAX, f"{temperature:.2f}", dl, len(path_list), *path_list])
    blocks, cur = [], []
    for line in (tmp_path / "part.dat").read_text().splitlines():
        tok = line.split()
        if tok and tok[0][0].isalpha():
            if cur:
                blocks.append(cur)
            cur = []
        elif tok:
            cur.extend(float(v) for v in tok)
    blocks.append(cur)
    td = thermodynamics(new)
    np.testing.assert_allclose(td.partition_function, blocks[0], rtol=1e-12)
    np.testing.assert_allclose(td.free_energy, blocks[1], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(td.internal_energy, blocks[2], rtol=1e-12)
    np.testing.assert_allclose(td.entropy, blocks[3], rtol=1e-10, atol=1e-10)


def test_legacy_edge_effect_script(tmp_path):
    """Compare :func:`derivative_analysis` with ``edge_effect_multi.py``."""
    from allostery import clean_hessian, derivative_analysis

    n = 24
    x, h = chain_system(n, seed=3)
    h_noisy = h.copy()
    h_noisy[0, 5] = h_noisy[5, 0] = 0.3  # a "negative spring" that must be removed
    src, snk = [0, 2], [20, 23]
    script = (LEGACY / "edge_effect_multi.py").read_text()
    script = (script.replace("N=454", f"N={n}")
              .replace("src=[49,103,129,224]", f"src={src}")
              .replace("snk=[303,336,430,432]", f"snk={snk}"))
    (tmp_path / "edge_effect.py").write_text(script)
    d = tmp_path / "sys_v02" / "paths"
    d.mkdir(parents=True)
    np.savetxt(tmp_path / "sys_v02" / "IGPS_sys_v02_1us_com_average_structure.dat", x)
    np.savetxt(tmp_path / "sys_v02" / "IGPS_sys_v02_1us_com_hess.dat", h_noisy)
    edges = [(i, j) for i in range(n) for j in range(i)]
    for s in src:
        for t in snk:
            np.savetxt(d / f"edge_sys_v02_1us_com_hess_adj.{s + 1:03d}.{t + 1:03d}.dat",
                       [[i + 1, j + 1, 0.0] for i, j in edges], fmt="%d %d %f")
    np.savetxt(d / "edge_sys_v02_1us_com_hess_adj.dat",
               [[i + 1, j + 1, 0.0] for i, j in edges], fmt="%d %d %f")
    pytest.importorskip("scipy")
    subprocess.run(["python3", "edge_effect.py", "sys", "hess"], cwd=tmp_path, check=True,
                   capture_output=True)
    old = np.loadtxt(d / "dC_sys_v02_1us_com_hess_adj.dat")

    hc = clean_hessian(h_noisy, x, cutoff=12.0)
    kT_legacy = 0.6  # hard-coded in the script
    res = derivative_analysis(hc, x, src, snk, temperature=kT_legacy / 1.9872e-3)
    i, j = old[:, 0].astype(int) - 1, old[:, 1].astype(int) - 1
    scale = np.abs(old[:, 2]).max()
    np.testing.assert_allclose(res.edge[i, j], old[:, 2], rtol=1e-6, atol=1e-9 * scale)
    finite = np.isfinite(old[:, 3])
    np.testing.assert_allclose(res.removal[i, j][finite], old[finite, 3], rtol=1e-6,
                               atol=1e-9 * np.abs(old[finite, 3]).max())
    np.testing.assert_allclose(res.spring_constants[i, j], old[:, 4], atol=1e-12)

    # node metric (edge_effect_node_multi.py) uses the uncleaned Hessian's
    # negative off-diagonals as weights, which equal the cleaned ones here
    node_old = -np.sum(res.edge * hc.clip(max=0.0), axis=0)
    np.testing.assert_allclose(res.node_raw, node_old, rtol=1e-12)
