"""Configuration-driven workflow (``allostery run config.toml``).

The workflow consists of up to five steps, each reading the results of the
previous ones from the output directory if they were not computed in the same
run:

``covariance``  covariance and average structure from a trajectory
``henm``        hENM Hessian from the covariance
``adjacency``   normalized adjacency for path sampling
``paths``       path sampling and analysis for every variant
``derivative``  derivative centrality metrics

See ``examples/config.toml`` for a documented configuration file.  The step
functions (:func:`step_covariance`, :func:`step_henm`, ...) can also be called
directly from Python.
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from . import __version__
from .adjacency import hessian_adjacency, lmi_adjacency, pearson_adjacency
from .derivative import derivative_analysis
from .henm import fit_henm
from .hessian import clean_hessian
from .io import load_matrix, save_matrix
from .paths import (
    VARIANTS,
    EnsembleSummary,
    PathSampler,
    degeneracy,
    distill,
    legacy,
    path_probabilities,
    thermodynamics,
)

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    import tomli as tomllib

STEPS = ("covariance", "henm", "adjacency", "paths", "derivative")

# file names inside the output directory
F_COV = "covariance.dat"
F_AVG = "average_structure.dat"
F_LABELS = "labels.txt"
F_HESS = "hessian.dat"
F_ADJ = "adjacency.dat"


def log(msg: str) -> None:
    print(f"[allostery {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------
# configuration
# ----------------------------------------------------------------------------
DEFAULTS: dict[str, Any] = {
    "output_directory": "allostery_results",
    "temperature": 303.0,
    "index_base": 0,
    "input": {
        "covariance": None, "average_structure": None, "hessian": None,
        "adjacency": None, "labels": None,
        "topology": None, "trajectory": None, "selection": "protein", "node": "com",
        "start": None, "stop": None, "step": None,
    },
    "henm": {
        "cutoff": None, "alpha": 1e-2, "max_iterations": 1000, "threshold": 1e-4,
        "initial_guess": 10.0, "n_zero_modes": 6,
    },
    "adjacency": {"kind": "hessian"},
    "paths": {
        "enabled": True, "sources": [], "sinks": [], "variants": ["direct"],
        "n_samples": 100_000, "temperature": 1.0, "seed": 2020, "n_top": 1000,
        "bin_width": 0.01, "n_equilibration": 0, "n_jobs": 1, "backend": "auto",
        "save_ensembles": False, "legacy_output": False,
    },
    "derivative": {
        "enabled": True, "sources": None, "sinks": None, "cutoff": None,
        "per_pair": False, "n_zero_modes": 6,
    },
}


@dataclass
class Config:
    data: dict
    base_dir: Path = field(default_factory=Path.cwd)

    @classmethod
    def from_file(cls, filename) -> "Config":
        filename = Path(filename)
        with open(filename, "rb") as f:
            raw = tomllib.load(f)
        return cls.from_dict(raw, base_dir=filename.resolve().parent)

    @classmethod
    def from_dict(cls, raw: dict, base_dir=None) -> "Config":
        data = _merge(DEFAULTS, raw)
        cfg = cls(data, Path(base_dir) if base_dir else Path.cwd())
        cfg.validate()
        return cfg

    def __getitem__(self, key):
        return self.data[key]

    def path(self, value) -> Path | None:
        if value is None:
            return None
        p = Path(value)
        return p if p.is_absolute() else self.base_dir / p

    @property
    def output(self) -> Path:
        return self.path(self.data["output_directory"])

    def validate(self) -> None:
        d = self.data
        if d["index_base"] not in (0, 1):
            raise ValueError("index_base must be 0 or 1")
        for v in d["paths"]["variants"]:
            if v not in VARIANTS:
                raise ValueError(f"unknown path variant {v!r}; choose from {list(VARIANTS)}")
        if d["adjacency"]["kind"] not in ("hessian", "pearson", "lmi"):
            raise ValueError("adjacency.kind must be 'hessian', 'pearson' or 'lmi'")


def _merge(defaults: dict, user: dict) -> dict:
    out = {}
    for k, v in defaults.items():
        if isinstance(v, dict):
            sub = user.get(k, {})
            if not isinstance(sub, dict):
                raise ValueError(f"[{k}] must be a table")
            unknown = set(sub) - set(v)
            if unknown:
                raise ValueError(f"unknown option(s) in [{k}]: {sorted(unknown)}")
            out[k] = {**v, **sub}
        else:
            out[k] = user.get(k, v)
    unknown = set(user) - set(defaults)
    if unknown:
        raise ValueError(f"unknown option(s): {sorted(unknown)}")
    return out


def _to0(nodes, base):
    return [int(v) - base for v in (nodes or [])]


# ----------------------------------------------------------------------------
# steps
# ----------------------------------------------------------------------------
def step_covariance(topology, trajectory, output: Path, selection="protein", node="com",
                    start=None, stop=None, step=None):
    """Covariance, average structure and residue labels from a trajectory."""
    from .trajectory import trajectory_covariance

    log(f"reading trajectory ({selection!r}, node={node})")
    cov, avg, labels = trajectory_covariance(topology, trajectory, selection, node,
                                             start, stop, step)
    output.mkdir(parents=True, exist_ok=True)
    save_matrix(output / F_COV, cov)
    save_matrix(output / F_AVG, avg, fmt="%.6f")
    (output / F_LABELS).write_text("\n".join(labels) + "\n")
    log(f"{len(labels)} nodes; wrote {F_COV}, {F_AVG}, {F_LABELS}")
    return cov, avg, labels


def step_henm(covariance, positions, temperature, output: Path, **options):
    """Fit the hENM Hessian and write it (plus the convergence history)."""
    log(f"fitting hENM Hessian for {positions.shape[0]} nodes at T = {temperature} K")
    res = fit_henm(covariance, positions, temperature=temperature, **options)
    output.mkdir(parents=True, exist_ok=True)
    save_matrix(output / F_HESS, res.hessian,
                header="hENM Hessian (kcal/mol/A^2): -k_ij off-diagonal")
    np.savetxt(output / "henm_residuals.dat",
               np.column_stack([np.arange(1, len(res.residuals) + 1), res.residuals]),
               fmt=["%d", "%.8e"], header="iteration residual")
    status = "converged" if res.converged else "NOT converged (max_iterations reached)"
    log(f"hENM {status} after {res.n_iterations} iterations "
        f"(residual {res.residuals[-1]:.3e}); wrote {F_HESS}")
    return res


def step_adjacency(kind: str, output: Path, hessian=None, covariance=None):
    if kind == "hessian":
        a = hessian_adjacency(hessian)
    elif kind == "pearson":
        a = pearson_adjacency(covariance)
    else:
        a = lmi_adjacency(covariance)
    output.mkdir(parents=True, exist_ok=True)
    save_matrix(output / F_ADJ, a, fmt="%.12e", header=f"normalized {kind} adjacency")
    log(f"wrote {kind} adjacency to {F_ADJ}")
    return a


def step_paths(adjacency, positions, sources, sinks, output: Path, variants=("direct",),
               n_samples=100_000, temperature=1.0, seed=2020, n_top=1000, bin_width=0.01,
               n_equilibration=0, n_jobs=1, backend="auto", save_ensembles=False,
               legacy_output=False, index_base=0, labels=None):
    """Sample paths for every source/sink pair and variant; write the analyses."""
    if not sources or not sinks:
        raise ValueError("paths: sources and sinks must be given")
    results = {}
    for variant in variants:
        vdir = output / "paths" / variant
        vdir.mkdir(parents=True, exist_ok=True)
        sampler = PathSampler.from_variant(variant, adjacency, positions,
                                           temperature=temperature, backend=backend)
        log(f"paths [{variant}]: {len(sources)} x {len(sinks)} pairs, {n_samples} samples "
            f"each ({sampler.backend} backend)")
        keep = bool(save_ensembles or legacy_output)
        t0 = time.time()
        ens = sampler.sample_pairs(sources, sinks, n_samples, seed=seed, keep_paths=keep,
                                   n_top=n_top, n_jobs=n_jobs, n_equilibration=n_equilibration)
        summaries = {}
        for (s, t), e in ens.items():
            tag = f"{s + index_base:04d}_{t + index_base:04d}"
            if keep:
                if save_ensembles:
                    e.save(vdir / f"ensemble_{tag}.npz")
                if legacy_output:
                    legacy.write_paths(vdir / f"paths_{tag}.dat", e)
                e = e.summary(n_top)
            e.save(vdir / f"summary_{tag}.npz")
            summaries[(s, t)] = e
        log(f"paths [{variant}]: sampling took {time.time() - t0:.1f} s")
        results[variant] = write_path_analyses(summaries, vdir, n_top, bin_width,
                                               index_base, labels)
    return results


def write_path_analyses(summaries, directory: Path, n_top=1000, bin_width=0.01,
                        index_base=0, labels=None) -> dict:
    """Write node/edge probabilities, distilled paths, degeneracy and thermodynamics."""
    directory = Path(directory)
    b = index_base
    s0 = next(iter(summaries.values()))
    n = s0.n_nodes
    probs = path_probabilities(summaries)

    node_cols = [np.arange(n) + b, probs.node]
    head = f"node(index_base={b}) P_node"
    if probs.node_by_segment is not None:
        node_cols += [probs.node_by_segment[0], probs.node_by_segment[1]]
        head += " P_source_segment P_sink_segment"
    _write_table(directory / "node_probability.dat", np.column_stack(node_cols), head,
                 ["%d"] + ["%.10e"] * (len(node_cols) - 1), labels)

    if probs.edge_by_kind is not None:
        e = probs.edge_by_kind
        i, j = np.nonzero(np.triu(e.sum(0), 1) > 1e-16)
        tab = np.column_stack([i + b, j + b, e.sum(0)[i, j], e[0][i, j], e[1][i, j], e[2][i, j]])
        head = "i j P_edge P_source_segment P_sink_segment P_jump"
    else:
        el = probs.edge_list()
        tab = np.column_stack([el[:, 0] + b, el[:, 1] + b, el[:, 2]])
        head = "i j P_edge"
    np.savetxt(directory / "edge_probability.dat", tab,
               fmt=["%d", "%d"] + ["%.10e"] * (tab.shape[1] - 2), header=head)

    top = distill(summaries, n_top)
    with open(directory / "top_paths.dat", "w") as f:
        f.write(f"# weight count probability source sink n_nodes nodes(index_base={b})"
                + (" -- source segment | sink segment\n" if s0.broken else "\n"))
        for p in top:
            nodes = [str(v + b) for v in p.nodes]
            if p.segments is not None:
                k = sum(1 for sg in p.segments if sg == 0)
                nodes.insert(k, "|")
            f.write(f"{p.weight:.10e} {p.count} {p.probability:.6e} {p.source + b} "
                    f"{p.sink + b} {len(p.nodes)} " + " ".join(nodes) + "\n")

    deg = degeneracy(summaries, bin_width)
    np.savetxt(directory / "degeneracy.dat",
               np.column_stack([deg.weight, deg.probability, deg.degeneracy]),
               fmt="%.8e", header="path_weight probability_density degeneracy")

    td = thermodynamics(summaries)
    tab = td.as_table()
    tab[:, :2] += b
    np.savetxt(directory / "thermodynamics.dat", tab,
               fmt=["%d", "%d", "%.10e", "%.10e", "%.10e", "%.10e"],
               header="source sink Z A=-T*lnZ U=<w> TS=U-A")
    log(f"wrote path analyses to {directory}")
    return {"probabilities": probs, "top_paths": top, "degeneracy": deg,
            "thermodynamics": td, "summaries": summaries}


def step_derivative(hessian, positions, sources, sinks, temperature, output: Path,
                    cutoff=None, per_pair=False, n_zero_modes=6, index_base=0, labels=None):
    if not sources or not sinks:
        raise ValueError("derivative: sources and sinks must be given")
    h = clean_hessian(hessian, positions, cutoff)
    log(f"derivative analysis: {len(sources)} x {len(sinks)} pairs")
    res = derivative_analysis(h, positions, sources, sinks, temperature=temperature,
                              n_zero_modes=n_zero_modes, per_pair=per_pair)
    d = output / "derivative"
    d.mkdir(parents=True, exist_ok=True)
    b = index_base
    n = h.shape[0]
    _write_table(d / "node_metric.dat",
                 np.column_stack([np.arange(n) + b, res.node, res.node_raw, res.node_removal]),
                 f"node(index_base={b}) delta_node sum_j(k_ij*delta_edge) sum_j(removal_ij)",
                 ["%d", "%.8e", "%.8e", "%.8e"], labels)
    el = res.edge_list(only_springs=False)
    el[:, :2] += b
    np.savetxt(d / "edge_metric.dat", el, fmt=["%d", "%d", "%.8e", "%.8e", "%.8e"],
               header="i j delta_edge removal k_ij")
    if per_pair:
        for (m, t), (de, dr) in res.per_pair.items():
            np.savez_compressed(d / f"pair_{m + b:04d}_{t + b:04d}.npz", edge=de, removal=dr)
    log(f"wrote derivative metrics to {d}")
    return res


def _write_table(filename, table, header, fmt, labels=None):
    if labels is None:
        np.savetxt(filename, table, fmt=fmt, header=header)
        return
    with open(filename, "w") as f:
        f.write(f"# {header} label\n")
        for row, lab in zip(table, labels):
            f.write(" ".join(fm % v for fm, v in zip(fmt, row)) + f" {lab}\n")


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------
def run(config: Config | dict | str | Path, steps=None) -> dict:
    """Run the workflow described by ``config`` (file name, dict or Config)."""
    if not isinstance(config, Config):
        config = (Config.from_dict(config) if isinstance(config, dict)
                  else Config.from_file(config))
    explicit_steps = steps is not None
    steps = set(STEPS if steps is None else steps)
    unknown = steps - set(STEPS)
    if unknown:
        raise ValueError(f"unknown step(s) {sorted(unknown)}; choose from {STEPS}")
    c, inp, out, b = config, config["input"], config.output, config["index_base"]
    out.mkdir(parents=True, exist_ok=True)
    results: dict[str, Any] = {}
    log(f"allostery {__version__}; output directory {out}")

    def first_existing(*candidates):
        for p in candidates:
            if p is not None and Path(p).exists():
                return Path(p)
        return None

    # --- covariance -------------------------------------------------------
    cov = avg = labels = None
    if inp["topology"] and "covariance" in steps:
        cov, avg, labels = step_covariance(
            c.path(inp["topology"]),
            [str(c.path(t)) for t in np.atleast_1d(inp["trajectory"])] if inp["trajectory"] else None,
            out, inp["selection"], inp["node"], inp["start"], inp["stop"], inp["step"])

    def need_avg():
        nonlocal avg
        if avg is None:
            f = first_existing(c.path(inp["average_structure"]), out / F_AVG)
            if f is None:
                raise FileNotFoundError("no average structure: set input.average_structure "
                                        "or input.topology/trajectory")
            avg = load_matrix(f)
        return avg

    def need_cov():
        nonlocal cov
        if cov is None:
            f = first_existing(c.path(inp["covariance"]), out / F_COV)
            if f is None:
                raise FileNotFoundError("no covariance: set input.covariance or "
                                        "input.topology/trajectory")
            cov = load_matrix(f)
        return cov

    if labels is None:
        f = first_existing(c.path(inp["labels"]), out / F_LABELS)
        if f is not None:
            labels = f.read_text().split()

    # --- hENM -------------------------------------------------------------
    hess = None
    if inp["hessian"]:
        hess = load_matrix(c.path(inp["hessian"]))
        log(f"using Hessian from {inp['hessian']}")
    elif "henm" in steps:
        opts = dict(c["henm"])
        if isinstance(opts["initial_guess"], str):
            opts["initial_guess"] = load_matrix(c.path(opts["initial_guess"]))
        results["henm"] = step_henm(need_cov(), need_avg(), c["temperature"], out, **opts)
        hess = results["henm"].hessian

    def need_hess():
        nonlocal hess
        if hess is None:
            f = first_existing(out / F_HESS)
            if f is None:
                raise FileNotFoundError("no Hessian: run the henm step or set input.hessian")
            hess = load_matrix(f)
        return hess

    # --- adjacency --------------------------------------------------------
    adj = None
    if inp["adjacency"]:
        adj = load_matrix(c.path(inp["adjacency"]))
    elif "adjacency" in steps and (explicit_steps or (c["paths"]["enabled"] and c["paths"]["sources"])):
        kind = c["adjacency"]["kind"]
        adj = step_adjacency(kind, out, hessian=need_hess() if kind == "hessian" else None,
                             covariance=need_cov() if kind != "hessian" else None)

    # --- paths ------------------------------------------------------------
    p = c["paths"]
    if "paths" in steps and p["enabled"] and p["sources"]:
        if adj is None:
            f = first_existing(out / F_ADJ)
            if f is None:
                raise FileNotFoundError("no adjacency matrix available")
            adj = load_matrix(f)
        needs_xyz = any(VARIANTS[v]["angular"] for v in p["variants"])
        results["paths"] = step_paths(
            adj, need_avg() if needs_xyz else None, _to0(p["sources"], b), _to0(p["sinks"], b),
            out, variants=p["variants"], n_samples=int(p["n_samples"]),
            temperature=p["temperature"], seed=p["seed"], n_top=p["n_top"],
            bin_width=p["bin_width"], n_equilibration=p["n_equilibration"], n_jobs=p["n_jobs"],
            backend=p["backend"], save_ensembles=p["save_ensembles"],
            legacy_output=p["legacy_output"], index_base=b, labels=labels)

    # --- derivative -------------------------------------------------------
    d = c["derivative"]
    src = d["sources"] if d["sources"] is not None else p["sources"]
    snk = d["sinks"] if d["sinks"] is not None else p["sinks"]
    if "derivative" in steps and d["enabled"] and src:
        results["derivative"] = step_derivative(
            need_hess(), need_avg(), _to0(src, b), _to0(snk, b), c["temperature"], out,
            cutoff=d["cutoff"], per_pair=d["per_pair"], n_zero_modes=d["n_zero_modes"],
            index_base=b, labels=labels)

    with open(out / "run_config.json", "w") as f:
        json.dump({"version": __version__, "steps": sorted(steps), "config": c.data}, f,
                  indent=2, default=str)
    log("done")
    return results
