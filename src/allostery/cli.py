"""Command-line interface: ``allostery <command> ...``.

Commands
--------
run         run a full workflow from a TOML configuration file
covariance  residue-node covariance from an MD trajectory (needs MDAnalysis)
henm        fit an hENM Hessian to a covariance matrix
adjacency   normalized adjacency (hessian, pearson or lmi)
paths       sample and analyse allosteric paths
derivative  derivative centrality metrics
analyze     re-analyse saved path summaries / legacy path files
compile     compile the Fortran path-sampling kernel
info        show version and kernel status

Run ``allostery <command> -h`` for the options of each command.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


def _nodes(values, base):
    return [int(v) - base for v in values]


def _add_index_base(p):
    p.add_argument("--index-base", type=int, choices=(0, 1), default=0,
                   help="numbering of nodes on the command line and in outputs (default 0)")


def cmd_run(a):
    from .pipeline import run

    run(a.config, steps=a.steps.split(",") if a.steps else None)


def cmd_covariance(a):
    from .pipeline import step_covariance

    step_covariance(a.topology, a.trajectory, Path(a.output), a.selection, a.node,
                    a.start, a.stop, a.step)


def cmd_henm(a):
    from .io import load_matrix
    from .pipeline import step_henm

    guess = a.initial_guess
    try:
        guess = float(guess)
    except ValueError:
        guess = load_matrix(guess)
    res = step_henm(load_matrix(a.covariance), load_matrix(a.structure), a.temperature,
                    Path(a.output), cutoff=a.cutoff, alpha=a.alpha,
                    max_iterations=a.max_iterations, threshold=a.threshold,
                    initial_guess=guess, verbose=a.verbose)
    return 0 if res.converged or not a.strict else 2


def cmd_adjacency(a):
    from .io import load_matrix
    from .pipeline import step_adjacency

    if a.kind == "hessian" and not a.hessian:
        sys.exit("--hessian is required for --kind hessian")
    if a.kind != "hessian" and not a.covariance:
        sys.exit("--covariance is required for --kind pearson/lmi")
    step_adjacency(a.kind, Path(a.output),
                   hessian=load_matrix(a.hessian) if a.hessian else None,
                   covariance=load_matrix(a.covariance) if a.covariance else None)


def _labels(path):
    return Path(path).read_text().split() if path else None


def cmd_paths(a):
    from .io import load_matrix
    from .pipeline import step_paths

    step_paths(load_matrix(a.adjacency), load_matrix(a.structure) if a.structure else None,
               _nodes(a.sources, a.index_base), _nodes(a.sinks, a.index_base), Path(a.output),
               variants=a.variant, n_samples=a.n_samples, temperature=a.temperature,
               seed=a.seed, n_top=a.n_top, bin_width=a.bin_width,
               n_equilibration=a.n_equilibration, n_jobs=a.n_jobs, backend=a.backend,
               save_ensembles=a.save_ensembles, legacy_output=a.legacy_output,
               index_base=a.index_base, labels=_labels(a.labels))


def cmd_derivative(a):
    from .io import load_matrix
    from .pipeline import step_derivative

    step_derivative(load_matrix(a.hessian), load_matrix(a.structure),
                    _nodes(a.sources, a.index_base), _nodes(a.sinks, a.index_base),
                    a.temperature, Path(a.output), cutoff=a.cutoff, per_pair=a.per_pair,
                    index_base=a.index_base, labels=_labels(a.labels))


def cmd_analyze(a):
    from .paths import legacy
    from .paths.ensemble import load
    from .pipeline import write_path_analyses

    summaries = {}
    for f in a.files:
        if f.endswith(".npz"):
            e = load(f)
        else:
            e = legacy.read_paths(f, n_nodes=a.n_nodes, temperature=a.temperature)
        if hasattr(e, "summary"):
            e = e.summary(a.n_top)
        summaries[(e.source, e.sink)] = e
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    write_path_analyses(summaries, out, a.n_top, a.bin_width, a.index_base, _labels(a.labels))


def cmd_compile(a):
    from .paths._fortran import compile_library

    from .paths._fortran import package_library_path

    out = a.output or package_library_path()
    try:
        path = compile_library(out, compiler=a.fc, verbose=True)
    except PermissionError:
        path = compile_library(None, compiler=a.fc, verbose=True)
    print(f"compiled {path}")


def cmd_info(a):
    import warnings

    from .paths import _fortran

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lib = _fortran.load(compile_if_missing=not a.no_compile)
    print(f"allostery {__version__}")
    print(f"Fortran kernel: {'available' if lib is not None else 'NOT available'}"
          + (f" ({lib._name})" if lib is not None
             else f" ({_fortran._load_error or 'not compiled; run `allostery compile`'})"))
    print(f"Fortran compiler: {_fortran.find_compiler()}")
    try:
        import MDAnalysis

        print(f"MDAnalysis: {MDAnalysis.__version__}")
    except ImportError:
        print("MDAnalysis: not installed (needed only for trajectories)")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="allostery", description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("run", help="run a workflow from a TOML configuration file")
    s.add_argument("config")
    s.add_argument("--steps", help="comma-separated subset of: covariance,henm,adjacency,"
                                   "paths,derivative")
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("covariance", help="node covariance from an MD trajectory")
    s.add_argument("topology")
    s.add_argument("trajectory", nargs="*")
    s.add_argument("--selection", default="protein")
    s.add_argument("--node", choices=("com", "cog", "ca"), default="com")
    s.add_argument("--start", type=int)
    s.add_argument("--stop", type=int)
    s.add_argument("--step", type=int)
    s.add_argument("-o", "--output", default=".", help="output directory")
    s.set_defaults(func=cmd_covariance)

    s = sub.add_parser("henm", help="fit an hENM Hessian")
    s.add_argument("--covariance", required=True, help="3N x 3N covariance (A^2)")
    s.add_argument("--structure", required=True, help="N x 3 average structure (A)")
    s.add_argument("-T", "--temperature", type=float, default=303.0, help="K")
    s.add_argument("--cutoff", type=float, help="spring distance cutoff (A); default none")
    s.add_argument("--alpha", type=float, default=1e-2)
    s.add_argument("--max-iterations", type=int, default=1000)
    s.add_argument("--threshold", type=float, default=1e-4)
    s.add_argument("--initial-guess", default="10.0",
                   help="uniform spring constant or a Hessian file to restart from")
    s.add_argument("--strict", action="store_true", help="exit with status 2 if not converged")
    s.add_argument("-v", "--verbose", action="store_true")
    s.add_argument("-o", "--output", default=".", help="output directory")
    s.set_defaults(func=cmd_henm)

    s = sub.add_parser("adjacency", help="normalized adjacency matrix")
    s.add_argument("--kind", choices=("hessian", "pearson", "lmi"), default="hessian")
    s.add_argument("--hessian")
    s.add_argument("--covariance")
    s.add_argument("-o", "--output", default=".", help="output directory")
    s.set_defaults(func=cmd_adjacency)

    s = sub.add_parser("paths", help="sample and analyse allosteric paths")
    s.add_argument("--adjacency", required=True)
    s.add_argument("--structure", help="N x 3 positions (needed for *_3d variants)")
    s.add_argument("--sources", nargs="+", required=True)
    s.add_argument("--sinks", nargs="+", required=True)
    s.add_argument("--variant", nargs="+", default=["direct"],
                   choices=("direct", "direct_3d", "broken", "broken_3d"))
    s.add_argument("-n", "--n-samples", type=int, default=100_000)
    s.add_argument("--temperature", type=float, default=1.0, help="path temperature")
    s.add_argument("--seed", type=int, default=2020)
    s.add_argument("--n-top", type=int, default=1000)
    s.add_argument("--bin-width", type=float, default=0.01)
    s.add_argument("--n-equilibration", type=int, default=0)
    s.add_argument("-j", "--n-jobs", type=int, default=1, help="parallel threads (-1: all)")
    s.add_argument("--backend", choices=("auto", "fortran", "python"), default="auto")
    s.add_argument("--save-ensembles", action="store_true", help="store every sampled path")
    s.add_argument("--legacy-output", action="store_true",
                   help="also write path files in the format of the original programs")
    s.add_argument("--labels", help="file with one node label per line")
    _add_index_base(s)
    s.add_argument("-o", "--output", default=".", help="output directory")
    s.set_defaults(func=cmd_paths)

    s = sub.add_parser("derivative", help="derivative centrality metrics")
    s.add_argument("--hessian", required=True)
    s.add_argument("--structure", required=True)
    s.add_argument("--sources", nargs="+", required=True)
    s.add_argument("--sinks", nargs="+", required=True)
    s.add_argument("-T", "--temperature", type=float, default=303.0, help="K")
    s.add_argument("--cutoff", type=float, help="drop springs longer than this (A)")
    s.add_argument("--per-pair", action="store_true")
    s.add_argument("--labels")
    _add_index_base(s)
    s.add_argument("-o", "--output", default=".", help="output directory")
    s.set_defaults(func=cmd_derivative)

    s = sub.add_parser("analyze", help="analyse saved summaries (.npz) or legacy path files")
    s.add_argument("files", nargs="+")
    s.add_argument("--n-nodes", type=int, help="number of nodes (legacy files)")
    s.add_argument("--temperature", type=float, default=1.0, help="path temperature (legacy files)")
    s.add_argument("--n-top", type=int, default=1000)
    s.add_argument("--bin-width", type=float, default=0.01)
    s.add_argument("--labels")
    _add_index_base(s)
    s.add_argument("-o", "--output", default=".")
    s.set_defaults(func=cmd_analyze)

    s = sub.add_parser("compile", help="compile the Fortran kernel")
    s.add_argument("--fc", help="Fortran compiler (default $FC or gfortran)")
    s.add_argument("-o", "--output", help="library file (default: inside the package)")
    s.set_defaults(func=cmd_compile)

    s = sub.add_parser("info", help="version and kernel status")
    s.add_argument("--no-compile", action="store_true")
    s.set_defaults(func=cmd_info)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    ret = args.func(args)
    return int(ret or 0)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
