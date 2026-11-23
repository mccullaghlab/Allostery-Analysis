"""Coarse-grained node positions from all-atom trajectories (needs MDAnalysis).

Install the optional dependency with ``pip install "allostery[md]"``.
"""

from __future__ import annotations

import numpy as np

from .covariance import covariance


def _mda():
    try:
        import MDAnalysis
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError("reading trajectories requires MDAnalysis: "
                          "pip install MDAnalysis") from exc
    return MDAnalysis


NODE_TYPES = ("com", "cog", "ca")


def node_positions(topology, trajectory=None, selection: str = "protein",
                   node: str = "com", start: int | None = None, stop: int | None = None,
                   step: int | None = None) -> tuple[np.ndarray, list[str]]:
    """Per-residue node coordinates for every frame.

    Parameters
    ----------
    topology, trajectory
        anything accepted by :class:`MDAnalysis.Universe` (a topology file and
        one or several trajectory files), or an existing ``Universe`` as
        ``topology``.
    selection
        atoms forming the coarse-grained system.
    node
        ``"com"`` (residue center of mass, used in the paper), ``"cog"``
        (center of geometry) or ``"ca"`` (alpha carbon).
    start, stop, step
        frame slice.

    Returns
    -------
    positions, labels
        ``(n_frames, N, 3)`` float array and residue labels such as
        ``"A:LEU50"``.
    """
    mda = _mda()
    if isinstance(topology, mda.Universe):
        u = topology
    else:
        u = mda.Universe(topology, *([] if trajectory is None else
                                     ([trajectory] if isinstance(trajectory, str) else trajectory)))
    if node not in NODE_TYPES:
        raise ValueError(f"node must be one of {NODE_TYPES}")
    atoms = u.select_atoms(selection)
    if node == "ca":
        atoms = atoms.select_atoms("name CA")
    residues = atoms.residues
    labels = [f"{r.segid or (r.atoms[0].chainID if hasattr(r.atoms[0], 'chainID') else '')}:"
              f"{r.resname}{r.resid}" for r in residues]
    frames = []
    for _ in u.trajectory[start:stop:step]:
        if node == "ca":
            frames.append(atoms.positions.astype(np.float64))
        elif node == "com":
            frames.append(atoms.center_of_mass(compound="residues"))
        else:
            frames.append(atoms.center_of_geometry(compound="residues"))
    if not frames:
        raise ValueError("no frames selected")
    return np.asarray(frames, dtype=np.float64), labels


def trajectory_covariance(topology, trajectory=None, selection: str = "protein",
                          node: str = "com", start=None, stop=None, step=None, **align_kwargs):
    """Covariance (``3N x 3N``), average structure (``N x 3``) and residue labels
    of residue nodes from a trajectory.  Frames are aligned on the nodes
    themselves (iteratively, to the average structure)."""
    pos, labels = node_positions(topology, trajectory, selection, node, start, stop, step)
    cov, avg = covariance(pos, **align_kwargs)
    return cov, avg, labels
