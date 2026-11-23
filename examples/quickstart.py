"""End-to-end example of the Python API on a synthetic system.

A random protein-like chain with a known elastic network plays the role of
the all-atom simulation: "trajectory" frames are drawn from the network's
Gaussian ensemble and randomly rotated/translated.  The example then

1. computes the node covariance (aligning the frames),
2. fits the hENM Hessian and compares it with the true one,
3. samples allosteric paths between two "sites" and analyses them,
4. computes the derivative centrality metrics.

Run with ``python examples/quickstart.py`` (takes ~10 s).
"""

import numpy as np

import allostery
from allostery.paths import PathSampler, distill, path_probabilities, thermodynamics

rng = np.random.default_rng(2020)
T = 303.0

# --------------------------------------------------------------------------
# a synthetic "protein": 60 nodes on a compact random chain
n = 60
x = [np.zeros(3)]
while len(x) < n:
    step = rng.normal(size=3)
    cand = x[-1] + 3.8 * step / np.linalg.norm(step)
    if np.linalg.norm(cand) < 12 and all(np.linalg.norm(cand - y) > 3.5 for y in x):
        x.append(cand)
x = np.array(x)
d = np.linalg.norm(x[:, None] - x[None], axis=-1)
k_true = np.where((d > 0) & (d < 9.0), 4.0 * np.exp(-(d - 3.8) / 2), 0.0)
for off, kk in ((1, 20.0), (2, 3.0), (3, 1.5)):
    i = np.arange(n - off)
    k_true[i, i + off] = k_true[i + off, i] = np.maximum(k_true[i, i + off], kk)
h_true = allostery.hessian.laplacian_from_springs(k_true)

cov_true = allostery.covariance_from_hessian(h_true, x, T)
frames = rng.multivariate_normal(x.ravel(), cov_true, size=20000, method="eigh").reshape(-1, n, 3)
for f in range(len(frames)):  # random rigid-body motion
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    frames[f] = frames[f] @ (q * np.sign(np.linalg.det(q))) + rng.normal(size=3) * 10

# --------------------------------------------------------------------------
# 1. covariance (for real data: allostery.trajectory.trajectory_covariance)
cov, avg = allostery.compute_covariance(frames)
print(f"covariance: {cov.shape}, average structure: {avg.shape}")

# 2. hENM
# (with a sampled covariance the residual reaches a noise floor, so the fit
# stops at max_iterations rather than at the threshold)
res = allostery.fit_henm(cov, avg, temperature=T, cutoff=12.0, alpha=0.05,
                         max_iterations=300)
k_fit = res.spring_constants
# the fit is done in the (arbitrary) frame of `avg`; spring constants are
# rotation invariant and can be compared directly
i = np.arange(n - 1)
print(f"hENM: {res.n_iterations} iterations, residual {res.residuals[-1]:.3g}; "
      f"backbone springs {k_fit[i, i + 1].mean():.1f} (true 20.0), "
      f"mean |k_fit - k_true| = {np.abs(k_fit - k_true)[k_true > 0].mean():.2f}")

# 3. paths on the normalized Hessian adjacency
adj = allostery.hessian_adjacency(res.hessian)
sampler = PathSampler(adj, avg, temperature=1.0)
# two "sites": node 0 and its neighbor, and the two nodes that are farthest
# from node 0 on the graph (largest shortest-path weight)
dist = [sampler.shortest_path(0, j)[1] if j else 0.0 for j in range(n)]
far = np.argsort(dist)[::-1]
sources, sinks = [0, 1], [int(far[0]), int(far[1])]
print(f"sources {sources}, sinks {sinks}")
ensembles = sampler.sample_pairs(sources, sinks, n_samples=20000, seed=1)
probs = path_probabilities(ensembles)
print("most visited nodes (P_node):",
      ", ".join(f"{i}:{p:.2f}" for i, p in probs.top_nodes(8)))
best = distill(ensembles, n_top=3)
for p in best:
    print(f"  path {p.nodes}  weight {p.weight:.3f}  probability {p.probability:.3f}")
td = thermodynamics(ensembles)
print("free energy per pair:", np.round(td.free_energy, 3))

# 4. derivative metrics
der = allostery.derivative_analysis(res.hessian, avg, sources, sinks, temperature=T)
order = np.argsort(der.node)[:5]
print("largest |delta_node|:", ", ".join(f"{i}:{der.node[i]:.2f}" for i in order))
