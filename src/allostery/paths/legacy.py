"""Readers/writers for the text formats of the original Fortran programs.

These allow analysing path files produced by the original code and writing
files that the original analysis programs can read.
"""

from __future__ import annotations

import numpy as np

from .ensemble import PathEnsemble


def read_paths(filename, source: int | None = None, sink: int | None = None,
               n_nodes: int | None = None, temperature: float = 1.0,
               angular: bool | None = None, broken: bool | None = None) -> PathEnsemble:
    """Read a path file written by ``paths_TWISP*.e``.

    Each line holds the path weight (1D) or the total, edge and angular
    weights (3D), the number of nodes and the 1-based nodes (negative for the
    sink segment of broken paths).  ``angular`` and ``broken`` are detected
    from the file if not given; ``source``, ``sink`` and ``n_nodes`` default
    to the first/last node and the largest node index.
    """
    weights, lengths, nodes = [], [], []
    with open(filename) as f:
        for line in f:
            tok = line.split()
            if not tok:
                continue
            if angular is None:
                angular = "." in tok[1] or "E" in tok[1].upper()
            nw = 3 if angular else 1
            w = [float(v) for v in tok[:nw]]
            npath = int(tok[nw])
            p = [int(v) for v in tok[nw + 1:nw + 1 + npath]]
            if len(p) != npath:
                raise ValueError(f"malformed line in {filename}: {line!r}")
            weights.append(w if angular else [w[0], w[0], 0.0])
            lengths.append(npath)
            nodes.extend(p)
    nodes = np.asarray(nodes, dtype=np.int64)
    if broken is None:
        broken = bool((nodes < 0).any())
    if n_nodes is None:
        n_nodes = int(np.abs(nodes).max()) if nodes.size else 0
    if source is None:
        source = abs(int(nodes[0])) - 1
    if sink is None:
        sink = abs(int(nodes[lengths[0] - 1])) - 1
    return PathEnsemble.from_signed(source, sink, n_nodes, temperature, bool(angular),
                                    bool(broken), np.asarray(weights).reshape(-1, 3),
                                    lengths, nodes)


def write_paths(filename, ensemble: PathEnsemble) -> None:
    """Write paths in the format of ``paths_TWISP*.e``."""
    width = 4 if ensemble.broken else 3
    signed = ensemble.signed_nodes()
    with open(filename, "w") as f:
        for k in range(len(ensemble)):
            if ensemble.angular:
                head = (f"{_e18(ensemble.weights[k])} {_e18(ensemble.edge_weights[k])} "
                        f"{_e18(ensemble.angular_weights[k])} ")
            else:
                head = f"{_e18(ensemble.weights[k])} "
            p = signed[ensemble.offsets[k]:ensemble.offsets[k + 1]]
            f.write(head + f"{p.size:{width}d}" + "".join(f" {v:{width}d}" for v in p) + " \n")


def _e18(x: float) -> str:
    """Fortran ``E18.12`` formatting (``0.ddddddddddddE+xx``, 12 significant digits)."""
    if x == 0:
        return "0.000000000000E+00"
    mant, e = f"{abs(x):.11e}".split("e")
    out = f"0.{mant.replace('.', '')}E{int(e) + 1:+03d}"
    return out if x > 0 else "-" + out


def write_adjacency(filename, adjacency) -> None:
    """Write an adjacency matrix in the fixed-width format read by the
    original programs (25 characters per entry)."""
    np.savetxt(filename, np.asarray(adjacency), fmt="%24.18e")
