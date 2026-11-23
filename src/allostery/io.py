"""Reading and writing matrices (text or ``.npy``)."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def load_matrix(filename) -> np.ndarray:
    """Load a matrix from ``.npy``, ``.npz`` (first array) or a text file."""
    p = Path(filename)
    if p.suffix == ".npy":
        return np.load(p)
    if p.suffix == ".npz":
        with np.load(p) as d:
            return d[d.files[0]]
    return np.loadtxt(p)


def save_matrix(filename, matrix, fmt: str = "%.10e", header: str = "") -> None:
    """Save a matrix; ``.npy`` is binary, anything else is text."""
    p = Path(filename)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.suffix == ".npy":
        np.save(p, matrix)
    else:
        np.savetxt(p, np.asarray(matrix), fmt=fmt, header=header)
