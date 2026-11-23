"""Basic operations on coarse-grained (residue-level) Hessians.

Conventions
-----------
* ``N`` nodes (e.g. residues), Cartesian coordinates ordered
  ``(x1, y1, z1, x2, ...)``, i.e. a ``3N x 3N`` matrix reshapes to
  ``(N, 3, N, 3)``.
* An *ENM Hessian* is stored as an ``N x N`` matrix ``H`` with off-diagonal
  elements ``H_ij = -k_ij`` (``k_ij`` the spring constant between nodes
  ``i`` and ``j``) and diagonal ``H_ii = sum_j k_ij``.  Its ``3N x 3N``
  form is ``H_ij = -k_ij  r_ij r_ij^T`` (eq. 2 of the paper) with
  ``r_ij`` the unit vector between the average positions.
* Units: kcal/mol, Angstrom, Kelvin.  Covariance = ``kT * pinv(H)``.
"""

from __future__ import annotations

import numpy as np

#: Boltzmann constant in kcal/(mol K) (value used in the original code)
KB_KCAL = 1.9872e-3


def kT(temperature: float) -> float:
    """Thermal energy in kcal/mol."""
    return KB_KCAL * float(temperature)


def n_nodes_from_3n(matrix: np.ndarray) -> int:
    n3 = matrix.shape[0]
    if matrix.ndim != 2 or matrix.shape[1] != n3 or n3 % 3:
        raise ValueError(f"expected a square 3N x 3N matrix, got shape {matrix.shape}")
    return n3 // 3


def unit_vectors(positions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pairwise unit vectors and distances.

    Returns ``(u, d)`` with ``u[i, j] = (x_i - x_j) / |x_i - x_j|`` (zero on
    the diagonal) and ``d[i, j] = |x_i - x_j|``.
    """
    x = np.asarray(positions, dtype=np.float64)
    if x.ndim != 2 or x.shape[1] != 3:
        raise ValueError("positions must have shape (N, 3)")
    diff = x[:, None, :] - x[None, :, :]
    d = np.sqrt(np.einsum("ijk,ijk->ij", diff, diff))
    inv = np.divide(1.0, d, out=np.zeros_like(d), where=d != 0)
    return diff * inv[:, :, None], d


def block_view(matrix3n: np.ndarray) -> np.ndarray:
    """View a ``3N x 3N`` matrix as ``(N, 3, N, 3)``."""
    n = n_nodes_from_3n(matrix3n)
    return matrix3n.reshape(n, 3, n, 3)


def diagonal_blocks(matrix3n: np.ndarray) -> np.ndarray:
    """``(N, 3, 3)`` diagonal blocks of a ``3N x 3N`` matrix."""
    return np.einsum("iaib->iab", block_view(matrix3n))


def laplacian_from_springs(k: np.ndarray) -> np.ndarray:
    """``N x N`` ENM Hessian from a symmetric spring-constant matrix."""
    k = np.array(k, dtype=np.float64)
    np.fill_diagonal(k, 0.0)
    h = -k
    np.fill_diagonal(h, k.sum(axis=1))
    return h


def spring_constants(hessian: np.ndarray) -> np.ndarray:
    """Spring constants ``k_ij = -H_ij`` (zero diagonal)."""
    k = -np.array(hessian, dtype=np.float64)
    np.fill_diagonal(k, 0.0)
    return k


def hessian_3n(hessian: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """Project an ``N x N`` ENM Hessian onto ``3N x 3N`` Cartesian form.

    Only the off-diagonal elements of ``hessian`` are used; the diagonal
    blocks are set such that every row of blocks sums to zero (translational
    invariance).
    """
    h = np.asarray(hessian, dtype=np.float64)
    n = h.shape[0]
    u, _ = unit_vectors(positions)
    if u.shape[0] != n:
        raise ValueError("hessian and positions have different numbers of nodes")
    off = np.array(h)
    np.fill_diagonal(off, 0.0)
    blocks = off[:, :, None, None] * u[:, :, :, None] * u[:, :, None, :]  # (N,N,3,3)
    diag = blocks.sum(axis=1)  # (N,3,3)
    h3 = blocks.transpose(0, 2, 1, 3).copy()  # (N,3,N,3)
    idx = np.arange(n)
    h3[idx, :, idx, :] = -diag
    return h3.reshape(3 * n, 3 * n)


def pseudo_inverse(matrix: np.ndarray, n_zero_modes: int | None = 6,
                   rcond: float | None = None) -> np.ndarray:
    """Moore-Penrose pseudo-inverse of a symmetric positive semi-definite matrix.

    By default the ``n_zero_modes`` (6: rigid translations and rotations)
    smallest eigenvalues are discarded, as in the original code.  If
    ``n_zero_modes`` is None, eigenvalues below ``rcond * max(eigenvalue)``
    (default ``1e-10``) are discarded instead.
    """
    e, v = np.linalg.eigh(matrix)
    if n_zero_modes is None:
        tol = (1e-10 if rcond is None else rcond) * max(abs(e[-1]), np.finfo(float).tiny)
        keep = np.abs(e) > tol
        e, v = e[keep], v[:, keep]
    else:
        e, v = e[n_zero_modes:], v[:, n_zero_modes:]
    return (v / e) @ v.T


def covariance_from_hessian(hessian: np.ndarray, positions: np.ndarray, temperature: float,
                            n_zero_modes: int | None = 6) -> np.ndarray:
    """``3N x 3N`` covariance ``C = kT H^+`` of an ENM Hessian (in A^2)."""
    return kT(temperature) * pseudo_inverse(hessian_3n(hessian, positions), n_zero_modes)


def projected_variance(cov3n: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Variance of each pair distance projected on its unit vector.

    ``V_ij = u_ij^T (C_ii + C_jj - C_ij - C_ji) u_ij``, i.e. the variance of
    ``(x_i - x_j) . u_ij``.  The diagonal is zero.
    """
    c = block_view(cov3n)
    d = np.einsum("iaib->iab", c)
    vi = np.einsum("ija,iab,ijb->ij", u, d, u, optimize=True)
    vj = np.einsum("ija,jab,ijb->ij", u, d, u, optimize=True)
    vij = np.einsum("ija,iajb,ijb->ij", u, c, u, optimize=True)
    return vi + vj - 2.0 * vij


def clean_hessian(hessian: np.ndarray, positions: np.ndarray | None = None,
                  cutoff: float | None = None) -> np.ndarray:
    """Remove non-physical springs from an ENM Hessian.

    Positive off-diagonal elements (negative spring constants) are removed,
    as are springs between nodes further apart than ``cutoff`` (if given),
    and the diagonal is reset to the row sums.  This is the preprocessing of
    ``edge_effect_multi.py`` (which used a 12 A cutoff).
    """
    k = spring_constants(hessian).clip(min=0.0)
    if cutoff is not None:
        if positions is None:
            raise ValueError("a distance cutoff requires positions")
        _, d = unit_vectors(positions)
        k[d > cutoff] = 0.0
    return laplacian_from_springs(k)
