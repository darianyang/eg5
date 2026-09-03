"""
Shared pathway-space clustering for the two Eg5 ADP-unbinding WE runs, using the
SHARED FIXED BIN GRID discretization (build_shared_grid.py + reassign_binned.py)
instead of the 6 agglomerative shared clusters.

This is the binned-grid counterpart of cluster_paths_shared.py: it pools the
successful pathways from both runs into ONE LPATH distance matrix, cuts ONE
merged dendrogram, and compares how each condition's pathway population is
distributed across the shared path clusters.

    nomon = WT Eg5 (no monastrol)      wmon = Eg5 + monastrol

The state alphabet (dictionary) is loaded from shared_grid.pkl, so a nomon path
string and a wmon path string mean the same thing under lpath.match.calc_dist
(LCS-based) -- the grid cells and the bound/unbound anchors are identical between
runs by construction.

Outputs (in shared_paths_binned/):
    distmat_shared.npy       pooled NxN pathway distance matrix
    dendrogram_shared.pdf    one merged, cluster-colored dendrogram
    path_distributions.pdf   grouped bars: weighted & count fractions per cluster
    path_cluster_summary.csv per-cluster counts/weights for both conditions
    path_assignments.csv     per-pathway origin + shared path-cluster id

Usage:
    python cluster_paths_shared_binned.py [n_clusters]
        n_clusters : optional int. If omitted, k is chosen from the largest gap
                     between the top merges in the merged dendrogram.
"""
import os
import sys
import pickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform

# LPATH lives outside the local namespace-package `lpath/` dir; import its metric
sys.path.insert(0, "/ihome/lchong/dty7/Apps/LPATH")
from lpath.match import calc_dist                      # noqa: E402
from sklearn.metrics import pairwise_distances         # noqa: E402

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = ["nomon", "wmon"]
CONDENSE = 1                       # must match run_lpath_binned.sh (--condense 1)
GRID_FILE = os.path.join(HERE, "shared_grid.pkl")
REASSIGNED = "succ_traj/reassigned_binned.pickle"
OUTDIR = "shared_paths_binned"
COLOR_LIST = ["tomato", "dodgerblue", "orchid", "orange", "mediumseagreen",
              "gold", "slateblue", "sienna", "turquoise", "hotpink"]
K_SEARCH_MAX = 12


class _Bar:
    """No-op stand-in for the tqdm bar calc_dist expects."""
    def update(self, n):
        pass


def load_dictionary():
    """The state->char dictionary is shared and fixed; load it from the grid."""
    with open(GRID_FILE, "rb") as f:
        grid = pickle.load(f)
    return dict(grid["dictionary"])


def load_run(run, dictionary):
    """Return (path_strings, weights) for one run.

    path_strings : list of the col-2 state-id sequences (LPATH metric input).
    weights      : WE weight of each pathway, taken (as gen_dist_matrix does)
                   from the last non-unknown frame's weight column.
    """
    with open(os.path.join(HERE, run, REASSIGNED), "rb") as f:
        re = np.asarray(pickle.load(f))
    unknown = len(dictionary) - 1
    path_strings, weights = [], []
    for pathway in re:
        pathway = np.asarray(pathway)
        nonzero = pathway[pathway[:, 2] < unknown]
        weights.append(float(nonzero[-1][-1]))
        path_strings.append(pathway[:, 2])
    return path_strings, np.asarray(weights)


def choose_k_from_gap(z, kmax):
    """Pick k from the largest gap between the top merges of the dendrogram."""
    d = np.sort(z[:, 2])
    best_k, best_gap, best_thr = 2, -np.inf, d[-1]
    kmax = min(kmax, len(d))
    for k in range(2, kmax + 1):
        upper = d[-(k - 1)]
        lower = d[-k]
        gap = upper - lower
        if gap > best_gap:
            best_k, best_gap, best_thr = k, gap, 0.5 * (upper + lower)
    return best_k, best_thr


def colored_dendrogram(z, labels, threshold, out_path):
    cluster_colors = [COLOR_LIST[(l - 1) % len(COLOR_LIST)] for l in labels]
    link_cols = {}
    n = len(z)
    for i, (a, b) in enumerate(z[:, :2].astype(int)):
        ca = link_cols[a] if a > n else cluster_colors[a]
        cb = link_cols[b] if b > n else cluster_colors[b]
        link_cols[i + 1 + n] = ca if ca == cb else "grey"

    fig, ax = plt.subplots(figsize=(9, 4.5))
    with plt.rc_context({"lines.linewidth": 1.5}):
        sch.dendrogram(z, no_labels=True, color_threshold=threshold,
                       link_color_func=lambda x: link_cols[x],
                       above_threshold_color="grey", ax=ax)
    ax.axhline(y=threshold, c="k", linestyle="--", linewidth=1.5)
    ax.set_ylabel("distance")
    ax.set_xlabel("pooled pathways (nomon + wmon)")
    ax.set_title("Shared pathway dendrogram (binned grid)")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def distribution_plot(summary, ks, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharex=True)
    x = np.arange(len(ks))
    w = 0.38
    for ax, key, title in [(axes[0], "wfrac", "Weighted (WE flux) fraction"),
                           (axes[1], "cfrac", "Count fraction")]:
        ax.bar(x - w / 2, [summary[k]["nomon"][key] for k in ks], w,
               label="nomon", color="dodgerblue")
        ax.bar(x + w / 2, [summary[k]["wmon"][key] for k in ks], w,
               label="wmon", color="tomato")
        ax.set_xticks(x)
        ax.set_xticklabels([f"P{k}" for k in ks])
        ax.set_xlabel("shared path cluster")
        ax.set_ylabel("fraction within condition")
        ax.set_title(title)
        ax.legend()
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main():
    outdir = os.path.join(HERE, OUTDIR)
    os.makedirs(outdir, exist_ok=True)

    # 1. shared dictionary (fixed grid alphabet, identical for both runs)
    dictionary = load_dictionary()
    print(f"Shared grid alphabet: {len(dictionary)} states (incl. unknown)")

    # 2. pool path strings + weights, tag origin
    path_strings, weights, origin = [], [], []
    for run in RUNS:
        ps, w = load_run(run, dictionary)
        path_strings.extend(ps)
        weights.extend(w)
        origin.extend([run] * len(ps))
        print(f"  {run}: {len(ps)} successful pathways, "
              f"weight sum = {w.sum():.3e}")
    weights = np.asarray(weights)
    origin = np.asarray(origin)
    n = len(path_strings)
    print(f"Pooled: {n} pathways")

    # 3. pooled distance matrix with LPATH's own metric (condense matches run).
    print("Computing pooled distance matrix ...")
    bar = _Bar()

    def _metric(i, j):
        return calc_dist(path_strings[int(i[0])], path_strings[int(j[0])],
                         dictionary, bar, CONDENSE)

    distmat = pairwise_distances(
        X=np.arange(n).reshape(-1, 1), metric=_metric, n_jobs=1)
    np.save(os.path.join(outdir, "distmat_shared.npy"), distmat)

    # 4. ward linkage on the condensed pooled matrix
    z = sch.linkage(squareform(distmat, checks=False), method="ward")

    # 5. choose k
    if len(sys.argv) > 1:
        k = int(sys.argv[1])
        d = np.sort(z[:, 2])
        threshold = 0.5 * (d[-(k - 1)] + d[-k]) if k >= 2 else d[-1]
        print(f"Using user-specified k = {k}")
    else:
        k, threshold = choose_k_from_gap(z, K_SEARCH_MAX)
        print(f"Auto-selected k = {k} from largest top-merge gap "
              f"(cut threshold = {threshold:.3f})")

    labels = sch.fcluster(z, t=k, criterion="maxclust")
    ks = sorted(np.unique(labels))

    # 6. per-cluster distributions, normalized WITHIN each condition
    tot_w = {run: weights[origin == run].sum() for run in RUNS}
    tot_n = {run: int((origin == run).sum()) for run in RUNS}
    summary = {}
    for cl in ks:
        summary[cl] = {}
        for run in RUNS:
            m = (labels == cl) & (origin == run)
            wsum = weights[m].sum()
            cnt = int(m.sum())
            summary[cl][run] = {
                "n": cnt, "wsum": wsum,
                "cfrac": cnt / tot_n[run] if tot_n[run] else 0.0,
                "wfrac": wsum / tot_w[run] if tot_w[run] else 0.0,
            }

    # 7. report
    print("\n=== Shared path-cluster populations "
          "(fractions normalized within each condition) ===")
    hdr = (f"{'clust':>5} | {'n_nomon':>7} {'n_wmon':>7} | "
           f"{'cnt%_nomon':>10} {'cnt%_wmon':>10} | "
           f"{'wt%_nomon':>10} {'wt%_wmon':>10}")
    print(hdr)
    print("-" * len(hdr))
    for cl in ks:
        s = summary[cl]
        print(f"{cl:>5} | {s['nomon']['n']:>7} {s['wmon']['n']:>7} | "
              f"{100*s['nomon']['cfrac']:>10.2f} {100*s['wmon']['cfrac']:>10.2f} | "
              f"{100*s['nomon']['wfrac']:>10.2f} {100*s['wmon']['wfrac']:>10.2f}")
    print(f"{'tot':>5} | {tot_n['nomon']:>7} {tot_n['wmon']:>7} | "
          f"{100:>10.2f} {100:>10.2f} | {100:>10.2f} {100:>10.2f}")

    # 8. CSVs
    with open(os.path.join(outdir, "path_cluster_summary.csv"), "w") as f:
        f.write("cluster,n_nomon,n_wmon,cntfrac_nomon,cntfrac_wmon,"
                "wsum_nomon,wsum_wmon,wtfrac_nomon,wtfrac_wmon\n")
        for cl in ks:
            s = summary[cl]
            f.write(f"{cl},{s['nomon']['n']},{s['wmon']['n']},"
                    f"{s['nomon']['cfrac']:.6f},{s['wmon']['cfrac']:.6f},"
                    f"{s['nomon']['wsum']:.6e},{s['wmon']['wsum']:.6e},"
                    f"{s['nomon']['wfrac']:.6f},{s['wmon']['wfrac']:.6f}\n")
    with open(os.path.join(outdir, "path_assignments.csv"), "w") as f:
        f.write("index,origin,path_cluster,weight\n")
        for i in range(n):
            f.write(f"{i},{origin[i]},{labels[i]},{weights[i]:.6e}\n")

    # 9. plots
    colored_dendrogram(z, labels, threshold,
                       os.path.join(outdir, "dendrogram_shared.pdf"))
    distribution_plot(summary, ks, os.path.join(outdir, "path_distributions.pdf"))

    print(f"\nDone. Wrote distmat_shared.npy, dendrogram_shared.pdf, "
          f"path_distributions.pdf, and CSVs to {OUTDIR}/")


if __name__ == "__main__":
    main()
