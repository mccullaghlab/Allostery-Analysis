# allostery — residue-level allostery from the effective coarse-grained Hessian

Python package (with a small Fortran kernel) implementing the analyses of

> P. T. Lake, R. B. Davidson, H. Klem, G. M. Hocky, M. McCullagh,
> *Residue-Level Allostery Propagates through the Effective Coarse-Grained Hessian*,
> J. Chem. Theory Comput. **16**, 3385–3395 (2020). [doi:10.1021/acs.jctc.9b01149](https://doi.org/10.1021/acs.jctc.9b01149)

1. **Covariance** of residue nodes (e.g. centers of mass) from an all-atom MD trajectory.
2. **hENM**: fit an effective harmonic (heterogeneous elastic network) Hessian to the covariance.
3. **Allosteric paths**: Monte Carlo sampling of paths between source and sink residues on
   the normalized Hessian adjacency, and the analyses built on them (node/edge
   probabilities `P_node`, lowest-weight paths, path degeneracy, path thermodynamics).
4. **Derivative metrics**: sensitivity of the source–sink covariance to every spring constant
   (`δ_edge`, eq. 5) and the derived node centrality (`δ_node`, eq. 6).

## Installation

```bash
pip install .            # or: pip install -e .  for development
pip install ".[md]"      # + MDAnalysis, to read trajectories
```

Python ≥ 3.9 and NumPy are the only requirements. The path sampler is written in Fortran:
it is compiled automatically with `gfortran` (or `$FC`) the first time it is needed and
cached in `~/.cache/allostery`. To compile it explicitly, or to check the installation:

```bash
allostery compile
allostery info
```

Without a Fortran compiler a pure-Python implementation of the same algorithm is used
(identical results, but ~50× slower).

## Quick start

### Command line: whole workflow from a configuration file

```bash
allostery run examples/config.toml
```

`examples/config.toml` documents every option. A minimal configuration:

```toml
output_directory = "results"
temperature = 303.0      # K
index_base = 1           # node numbers below and in the outputs are 1-based

[input]
topology = "protein.prmtop"
trajectory = ["prod.dcd"]
# or: covariance = "cov.dat" and average_structure = "avg.dat"

[henm]
cutoff = 12.0

[paths]
sources = [50, 104, 130, 225]
sinks = [304, 337, 431, 433]
variants = ["direct", "broken"]
n_samples = 1000000
n_jobs = -1
```

Output (in `results/`):

| file | content |
|---|---|
| `covariance.dat`, `average_structure.dat`, `labels.txt` | from the trajectory |
| `hessian.dat`, `henm_residuals.dat` | N×N hENM Hessian (−k<sub>ij</sub> off-diagonal), convergence |
| `adjacency.dat` | normalized adjacency used for the paths |
| `paths/<variant>/node_probability.dat` | `P_node` of every node |
| `paths/<variant>/edge_probability.dat` | probability of every traversed edge |
| `paths/<variant>/top_paths.dat` | lowest-weight unique paths with counts and probabilities |
| `paths/<variant>/degeneracy.dat` | path-weight distribution and degeneracy (Fig. 3A) |
| `paths/<variant>/thermodynamics.dat` | Z, A, U, TS for every source/sink pair |
| `paths/<variant>/summary_*.npz` | per-pair summaries (re-analysable with `allostery analyze`) |
| `derivative/node_metric.dat`, `derivative/edge_metric.dat` | `δ_node`, `δ_edge`, spring-removal effects |

Each step can also be run on its own:

```bash
allostery covariance protein.prmtop prod.dcd --node com -o results
allostery henm --covariance results/covariance.dat --structure results/average_structure.dat \
               -T 303 --cutoff 12 -o results
allostery adjacency --hessian results/hessian.dat -o results
allostery paths --adjacency results/adjacency.dat --structure results/average_structure.dat \
                --sources 50 104 --sinks 304 337 --index-base 1 \
                --variant direct direct_3d -n 1000000 -j -1 -o results
allostery derivative --hessian results/hessian.dat --structure results/average_structure.dat \
                     --sources 50 104 --sinks 304 337 --index-base 1 -o results
allostery analyze results/paths/direct/summary_*.npz -o reanalysis
```

### Python

```python
import allostery
from allostery.paths import PathSampler, path_probabilities, distill

cov, avg, labels = allostery.trajectory.trajectory_covariance("top.prmtop", "traj.dcd")
henm = allostery.fit_henm(cov, avg, temperature=303, cutoff=12.0)
adj = allostery.hessian_adjacency(henm.hessian)

sampler = PathSampler(adj, avg, temperature=1.0, angular=False, broken=False)
ens = sampler.sample_pairs(sources=[49, 103], sinks=[303, 336], n_samples=10**6, seed=1, n_jobs=-1)
probs = path_probabilities(ens)          # probs.node, probs.edge
top = distill(ens, n_top=100)            # lowest-weight paths

der = allostery.derivative_analysis(henm.hessian, avg, [49, 103], [303, 336], temperature=303)
der.node, der.edge                       # δ_node, δ_edge
```

Python indices are always 0-based. See `examples/quickstart.py` for a complete,
self-contained example on a synthetic system.

## Path variants

| variant | angular term | broken | original program |
|---|---|---|---|
| `direct` | – | – | `paths_TWISP.f` |
| `direct_3d` | −Σ ln\|cos θ\| | – | `paths_TWISP_3D.f` |
| `broken` | – | ✓ | `paths_TWISP_broken.f` |
| `broken_3d` | −Σ ln\|cos θ\| | ✓ | `paths_TWISP_broken_3D.f` |

Paths are sampled from P(path) ∝ exp(−w/T) with w = Σ −ln A<sub>ij</sub> over the edges
(plus the angular term for `*_3d`). *Broken* paths consist of a segment attached to the
source and one attached to the sink (contributions through the null space of the Hessian).
Ensembles of different source/sink pairs are combined with weights given by their estimated
partition functions Z = (N / N<sub>min</sub>) exp(−w<sub>min</sub>/T).

## Package layout

```
src/allostery/
    covariance.py   alignment, covariance, Pearson and linear-MI adjacencies
    trajectory.py   residue nodes from MD trajectories (MDAnalysis)
    hessian.py      N×N ↔ 3N×3N Hessians, pseudo-inverse, cleaning
    henm.py         hENM fit
    adjacency.py    normalized adjacency matrices
    derivative.py   derivative centrality metrics
    paths/          path sampler (Fortran kernel _twisp.f90 + Python reference),
                    ensembles, analyses, legacy file formats
    pipeline.py     configuration-driven workflow
    cli.py          `allostery` command
tests/              test suite (pytest)
legacy/             the original code, kept for reference and regression tests
examples/           configuration template and Python example
```

## Testing

```bash
pip install -e ".[test]"
pytest
```

The tests include regression tests against the original code. They compile the Fortran
programs in `legacy/paths` for a small system, and require identical sampled paths and
identical results for every analysis program. Other tests check against exact
enumeration that the sampler reproduces the Boltzmann distribution of paths. The hENM fit
must recover known spring constants. The derivative metrics are checked against finite
differences.

## Changes with respect to the original code

* All fixed sizes (`npar=454`), file-name conventions and `input.dat`/`ran3.dat` files
  are gone. The 20 Fortran programs are replaced by one sampling kernel and NumPy analyses.
* **hENM bug fix**: the distance-cutoff mask in `henm_functions.py` was inverted. It kept
  only springs *beyond* the cutoff, so with the default cutoff of 9999 Å every spring was
  zeroed after the first iteration. `create_hessian.py` also did not run under Python 3
  (mixed tabs/spaces, float `nNodes`).
* The hENM fit detects updates that make the network floppy (more than 6 zero modes)
  and halves the step size instead of silently producing NaNs.
* `paths_deg*` reported each histogram bin one bin width too low. Bins are now
  reported at their centers.
* The derivative analysis uses k<sub>B</sub>T from the given temperature (the original
  script hard-coded 0.6 kcal/mol) and is vectorized.

## Citation

If you use this code, please cite the paper above.
