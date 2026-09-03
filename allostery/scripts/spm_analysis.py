"""Shortest Path Map (SPM) allosteric analysis per LPATH cluster.

Reimplements the SPM method of Rodriguez-Santos et al. (JCIM 2026) for the two
Eg5 ADP-unbinding WE simulations, resolved by the 6 shared LPATH clusters, and
compares the SPM between nomon (WT) and wmon (+monastrol) at each cluster.

Pipeline (per condition x cluster), on the uniform Ca ensemble from
extract_ca_coords.py (unweighted, per the current choice):

  1. iteratively Ca-superpose all frames onto their mean structure (Kabsch),
  2. displacement correlation (paper eq 2):
         C_ij = <dr_i . dr_j> / sqrt(<|dr_i|^2><|dr_j|^2>)
     where dr_i is residue i's Ca displacement from the mean structure,
  3. build a residue graph: an edge i-j exists iff the mean Ca-Ca distance is
     < ``contact`` A; its length is l_ij = -log(|C_ij|) (paper eq 3),
  4. edge-betweenness centrality (weighted shortest paths through each edge) is
     the SPM score -- the most-traversed edges/nodes are the allosteric map.

Outputs, cached to ``<cache-dir>/<cond>/cluster<c>.npz`` (C, dmean, edge
betweenness matrix, node betweenness) and figures to ``<out-dir>/``:
  * per-cluster SPM edge-betweenness maps, nomon vs wmon vs difference,
  * per-residue node-betweenness (SPM node importance), nomon vs wmon,
  * ChimeraX ``.defattr`` files colouring the structure by node betweenness.

    python scripts/spm_analysis.py [--data-dir data] [--out-dir plots/spm]
        [--cache-dir spm] [--contact 6.0] [--top-frac 0.05] [--recompute]
"""
import os
import sys
import argparse
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402
from df_analysis import EXIT_ORDER, REGIONS, add_region_ticks  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
import networkx as nx  # noqa: E402

CONDS = ("nomon", "wmon")


def kabsch_to_mean(coords, iters=2):
    """Ca-superpose all frames onto their (iteratively refined) mean structure.

    ``coords`` (n, R, 3) in Angstrom.  Returns the aligned coordinates.  Fully
    batched: per-frame optimal rotations come from a stacked 3x3 SVD.
    """
    X = coords - coords.mean(1, keepdims=True)          # remove translation
    ref = X[0]
    for _ in range(iters):
        ref = ref - ref.mean(0)
        H = np.einsum("nri,rj->nij", X, ref)            # (n,3,3) covariance
        U, _S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(np.einsum("nij,njk->nik",
                                            Vt.transpose(0, 2, 1),
                                            U.transpose(0, 2, 1))))
        D = np.zeros((len(X), 3, 3))
        D[:, 0, 0] = D[:, 1, 1] = 1.0
        D[:, 2, 2] = d
        R = np.einsum("nij,njk,nkl->nil", Vt.transpose(0, 2, 1), D,
                      U.transpose(0, 2, 1))              # (n,3,3)
        X = np.einsum("nri,nij->nrj", X, R)
        ref = X.mean(0)
    return X


def mean_distance(coords, chunk=128):
    """Mean Ca-Ca distance matrix (R, R) in Angstrom (rotation invariant)."""
    n, R, _ = coords.shape
    sumd = np.zeros((R, R), np.float64)
    for s in range(0, n, chunk):
        x = coords[s:s + chunk].astype(np.float64)
        nrm = np.einsum("mrc,mrc->mr", x, x)
        gram = x @ x.transpose(0, 2, 1)
        d2 = nrm[:, :, None] + nrm[:, None, :] - 2.0 * gram
        sumd += np.sqrt(np.maximum(d2, 0.0)).sum(0)
    return sumd / n


def displacement_correlation(aligned):
    """Paper eq 2: C_ij = <dr_i.dr_j>/sqrt(<|dr_i|^2><|dr_j|^2>)."""
    ref = aligned.mean(0)
    disp = aligned - ref                                 # (n, R, 3)
    n = len(disp)
    R = disp.shape[1]
    cross = np.zeros((R, R), np.float64)
    for c in range(3):                                   # sum over x,y,z
        Dc = disp[:, :, c]
        cross += Dc.T @ Dc
    cross /= n
    var = np.diag(cross).copy()
    denom = np.sqrt(np.outer(var, var))
    C = cross / np.where(denom > 0, denom, 1.0)
    return C


def spm_scores(C, dmean, contact, eps=1e-6):
    """Edge- and node-betweenness on the -log|C| contact graph."""
    R = C.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(R))
    iu, iv = np.triu_indices(R, k=1)
    mask = (dmean[iu, iv] < contact)                     # Ca contacts only
    for i, j in zip(iu[mask], iv[mask]):
        c = abs(C[i, j])
        length = -np.log(max(c, eps))                    # eq 3 (>= 0)
        G.add_edge(int(i), int(j), weight=length)
    ebc = nx.edge_betweenness_centrality(G, weight="weight", normalized=True)
    nbc = nx.betweenness_centrality(G, weight="weight", normalized=True)
    E = np.zeros((R, R))
    for (i, j), v in ebc.items():
        E[i, j] = E[j, i] = v
    node = np.array([nbc[i] for i in range(R)])
    return E, node, G.number_of_edges()


def compute_all(data_dir, cache_dir, contact, recompute):
    spm, res_ids, counts = {}, None, {}
    for cond in CONDS:
        h5 = f"{data_dir}/ca_coords_{cond}.h5"
        with h5py.File(h5, "r") as f:
            cl = f["cluster"][:]
            res_ids = f["res_ids"][:]
            for c in range(N_CLUSTERS):
                cache = f"{cache_dir}/{cond}/cluster{c}.npz"
                sel = cl == c
                counts[(cond, c)] = int(sel.sum())
                if os.path.exists(cache) and not recompute:
                    d = np.load(cache)
                    spm[(cond, c)] = (d["E"], d["node"])
                    continue
                if sel.sum() < 10:
                    continue
                ca = f["ca"][:][sel] * 10.0              # nm -> Angstrom
                aligned = kabsch_to_mean(ca)
                C = displacement_correlation(aligned)
                dmean = mean_distance(ca[::5])           # contacts converge fast
                E, node, nedge = spm_scores(C, dmean, contact)
                os.makedirs(os.path.dirname(cache), exist_ok=True)
                np.savez_compressed(cache, C=C, dmean=dmean, E=E, node=node)
                spm[(cond, c)] = (E, node)
                print(f"  SPM {cond} c{c}: n={sel.sum()} edges={nedge} -> {cache}",
                      flush=True)
    return spm, res_ids, counts


def write_chimerax_attr(path, res_ids, node, name):
    """ChimeraX .defattr: colour residues by SPM node betweenness."""
    with open(path, "w") as f:
        f.write(f"attribute: {name}\nmatch mode: 1-to-1\nrecipient: residues\n")
        for r, v in zip(res_ids, node):
            f.write(f"\t/A:{int(r)}\t{v:.6f}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default="plots/spm")
    ap.add_argument("--cache-dir", default="spm")
    ap.add_argument("--contact", type=float, default=6.0,
                    help="mean Ca-Ca distance (A) defining a graph edge")
    ap.add_argument("--top-frac", type=float, default=0.05,
                    help="fraction of top edges kept when reporting the SPM")
    ap.add_argument("--recompute", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    attr_dir = f"{a.out_dir}/chimerax"
    os.makedirs(attr_dir, exist_ok=True)

    spm, res_ids, counts = compute_all(a.data_dir, a.cache_dir, a.contact,
                                       a.recompute)
    if not spm:
        print("no SPM computed (missing ca_coords_*.h5?)")
        return
    extent = [res_ids[0], res_ids[-1], res_ids[-1], res_ids[0]]
    order = [c for c in EXIT_ORDER if any((cond, c) in spm for cond in CONDS)]

    allE = np.concatenate([E[E > 0].ravel() for E, _ in spm.values()])
    emax = np.percentile(allE, 99.5)

    # ---- grand grid: SPM edge-betweenness maps + drug difference ----
    fig, axes = plt.subplots(len(order), 3, figsize=(14, 4.4 * len(order)),
                             squeeze=False)
    for row, c in enumerate(order):
        for col, cond in enumerate(CONDS):
            ax = axes[row][col]
            if (cond, c) not in spm:
                ax.set_visible(False); continue
            E = spm[(cond, c)][0]
            im = ax.imshow(E, cmap="magma", vmin=0, vmax=emax, extent=extent,
                           interpolation="nearest")
            ax.set_title(f"{cond}  cluster {c}  (n={counts[(cond,c)]})", fontsize=9)
            add_region_ticks(ax, res_ids)
        axd = axes[row][2]
        if (("nomon", c) in spm) and (("wmon", c) in spm):
            D = spm[("wmon", c)][0] - spm[("nomon", c)][0]
            dmax = np.percentile(np.abs(D), 99.5) or 1.0
            imd = axd.imshow(D, cmap="RdBu_r", norm=TwoSlopeNorm(0, -dmax, dmax),
                             extent=extent, interpolation="nearest")
            axd.set_title(f"$\\Delta$SPM (wmon$-$nomon)  c{c}", fontsize=9)
            add_region_ticks(axd, res_ids)
            fig.colorbar(imd, ax=axd, fraction=0.046, label="$\\Delta$betweenness")
        else:
            axd.set_visible(False)
    fig.colorbar(im, ax=axes[:, :2], label="edge betweenness (SPM)",
                 fraction=0.02, pad=0.01)
    fig.suptitle("Eg5 ADP-exit Shortest Path Map per LPATH cluster\n"
                 "(rows = exit order 2→5→0→3→4→1; red = pathway gained with "
                 "monastrol)", y=0.997, fontsize=12)
    out = f"{a.out_dir}/spm_grand_grid.png"
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)

    # ---- per-residue node betweenness (SPM node importance) ----
    fig, axes = plt.subplots(len(order), 1, figsize=(12, 2.2 * len(order)),
                             squeeze=False, sharex=True)
    for row, c in enumerate(order):
        ax = axes[row][0]
        for cond, color in zip(CONDS, ("tab:blue", "tab:red")):
            if (cond, c) in spm:
                ax.plot(res_ids, spm[(cond, c)][1], lw=1, color=color, label=cond)
        ax.set_ylabel(f"c{c}", fontsize=9)
        for name, (lo, hi) in REGIONS.items():
            ax.axvspan(lo, hi, color="orange", alpha=0.08)
        if row == 0:
            ax.legend(fontsize=8)
    axes[-1][0].set_xlabel("residue (AMBER)")
    fig.suptitle("SPM node betweenness per residue (blue nomon, red wmon)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = f"{a.out_dir}/spm_node_betweenness.png"
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)

    # ---- ChimeraX attribute files + top-edge tables ----
    for (cond, c), (E, node) in spm.items():
        write_chimerax_attr(f"{attr_dir}/{cond}_cluster{c}_nodebetween.defattr",
                            res_ids, node, f"spmNode{cond}c{c}")
    for c in order:
        if ("nomon", c) in spm and ("wmon", c) in spm:
            D = spm[("wmon", c)][0] - spm[("nomon", c)][0]
            iu, iv = np.triu_indices(len(res_ids), k=1)
            vals = D[iu, iv]
            top = np.argsort(np.abs(vals))[::-1][:25]
            with open(f"{a.out_dir}/top_drug_edges_cluster{c}.csv", "w") as f:
                f.write("res_i,res_j,delta_betweenness\n")
                for k in top:
                    f.write(f"{int(res_ids[iu[k]])},{int(res_ids[iv[k]])},"
                            f"{vals[k]:.6e}\n")
    print(f"wrote ChimeraX attrs to {attr_dir}/ and top-edge tables to {a.out_dir}/")


if __name__ == "__main__":
    main()
