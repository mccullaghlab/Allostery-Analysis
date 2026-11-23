"""allostery -- residue-level allostery from the effective coarse-grained Hessian.

Implementation of the methods of

    P. T. Lake, R. B. Davidson, H. Klem, G. M. Hocky, M. McCullagh,
    "Residue-Level Allostery Propagates through the Effective Coarse-Grained
    Hessian", J. Chem. Theory Comput. 16, 3385-3395 (2020).

Workflow
--------
1. covariance of residue nodes from an all-atom trajectory
   (:mod:`allostery.trajectory`, :mod:`allostery.covariance`);
2. effective harmonic Hessian by hENM fitting (:func:`fit_henm`);
3. allosteric path sampling on the normalized Hessian adjacency
   (:class:`PathSampler`, :mod:`allostery.paths`);
4. derivative centrality metrics (:func:`derivative_analysis`).

The whole workflow can be driven from a TOML file with ``allostery run``
(see :mod:`allostery.pipeline`).
"""

from . import adjacency, covariance, derivative, hessian, henm, io, paths
from .adjacency import edge_costs, hessian_adjacency, lmi_adjacency, pearson_adjacency
from .covariance import align
from .covariance import covariance as compute_covariance
from .derivative import DerivativeResult, derivative_analysis
from .henm import HENMResult, fit_henm
from .hessian import (
    KB_KCAL,
    clean_hessian,
    covariance_from_hessian,
    hessian_3n,
    kT,
    pseudo_inverse,
)
from .paths import PathSampler

__version__ = "1.0.0"

__all__ = [
    "fit_henm", "HENMResult", "derivative_analysis", "DerivativeResult", "PathSampler",
    "hessian_adjacency", "pearson_adjacency", "lmi_adjacency", "edge_costs",
    "compute_covariance", "align", "clean_hessian", "covariance_from_hessian", "hessian_3n",
    "pseudo_inverse", "kT", "KB_KCAL",
    "adjacency", "covariance", "derivative", "hessian", "henm", "io", "paths",
]
