"""Sequence x sequence NMI heatmaps of the allosteric networks.

Each mdpath network is a set of per-residue-pair NMI edge weights on the shared
5 A-proximity graph, i.e. an *NMI-weighted contact map*.  This renders, per
LPATH cluster (ADP-exit stage):

  * nomon and wmon maps side by side (shared colour scale), and
  * their difference wmon - nomon (diverging), which is monastrol's allosteric
    footprint at that stage -- red = communication gained with drug, blue = lost.

A grand 6x3 grid (clusters x [nomon, wmon, delta]) gives the whole story in one
figure; per-cluster figures are also written for detail.

    python scripts/plot_nmi_heatmaps.py [--scheme pooled] [--angle-note phi]
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

NET = "networks"
PLOTS = "plots/nmi_heatmaps"


def load_nmi_matrix(cond, cluster, scheme, res_index):
    """Dense symmetric NMI matrix on the fixed residue index, or None."""
    path = f"{NET}/{cond}/cluster{cluster}_{scheme}/nmi.csv"
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    n = len(res_index)
    M = np.zeros((n, n), dtype=np.float64)
    for pair, val in zip(df["Residue Pair"], df["MI Difference"]):
        u, v = pair.split("|")
        iu = res_index.get(int(u.split()[1]))
        iv = res_index.get(int(v.split()[1]))
        if iu is not None and iv is not None:
            M[iu, iv] = val
    return M


def residue_universe(scheme):
    """Union of all residues appearing in any network (sorted), as an index map."""
    res = set()
    for cond in ("nomon", "wmon"):
        for c in range(N_CLUSTERS):
            p = f"{NET}/{cond}/cluster{c}_{scheme}/nmi.csv"
            if not os.path.exists(p):
                continue
            df = pd.read_csv(p, usecols=["Residue Pair"])
            for pair in df["Residue Pair"]:
                u, v = pair.split("|")
                res.add(int(u.split()[1]))
                res.add(int(v.split()[1]))
    res = sorted(res)
    return res, {r: i for i, r in enumerate(res)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scheme", default="pooled")
    args = ap.parse_args()
    scheme = args.scheme
    os.makedirs(PLOTS, exist_ok=True)

    reslist, ridx = residue_universe(scheme)
    extent = [reslist[0], reslist[-1], reslist[-1], reslist[0]]
    print(f"{len(reslist)} residues in the graph ({reslist[0]}-{reslist[-1]})")

    mats = {}  # (cond,cluster) -> matrix
    for cond in ("nomon", "wmon"):
        for c in range(N_CLUSTERS):
            M = load_nmi_matrix(cond, c, scheme, ridx)
            if M is not None:
                mats[(cond, c)] = M

    # shared NMI colour scale (99.5th pct of nonzero over all networks)
    allvals = np.concatenate([M[M > 0].ravel() for M in mats.values()])
    vmax = np.percentile(allvals, 99.5)
    # shared diverging scale for the differences
    diffs = {}
    for c in range(N_CLUSTERS):
        if ("nomon", c) in mats and ("wmon", c) in mats:
            diffs[c] = mats[("wmon", c)] - mats[("nomon", c)]
    dmax = np.percentile(np.abs(np.concatenate([D[D != 0].ravel()
                                                for D in diffs.values()])), 99.5)
    dnorm = TwoSlopeNorm(vmin=-dmax, vcenter=0.0, vmax=dmax)

    # ---- grand grid: rows = clusters, cols = [nomon, wmon, delta] ----
    fig, axes = plt.subplots(N_CLUSTERS, 3, figsize=(13.5, 4.2 * N_CLUSTERS))
    for c in range(N_CLUSTERS):
        for col, cond in enumerate(("nomon", "wmon")):
            ax = axes[c, col]
            M = mats.get((cond, c))
            if M is None:
                ax.set_visible(False)
                continue
            im = ax.imshow(M, cmap="magma", vmin=0, vmax=vmax,
                           extent=extent, interpolation="nearest")
            ax.set_title(f"{cond}  cluster {c}", fontsize=10)
            if col == 0:
                ax.set_ylabel("residue")
            if c == N_CLUSTERS - 1:
                ax.set_xlabel("residue")
        axd = axes[c, 2]
        if c in diffs:
            imd = axd.imshow(diffs[c], cmap="RdBu_r", norm=dnorm,
                             extent=extent, interpolation="nearest")
            axd.set_title(f"$\\Delta$ (wmon $-$ nomon)  cluster {c}", fontsize=10)
            if c == N_CLUSTERS - 1:
                axd.set_xlabel("residue")
        else:
            axd.set_visible(False)
    fig.colorbar(im, ax=axes[:, :2], label=f"NMI ({scheme})",
                 fraction=0.02, pad=0.01)
    fig.colorbar(imd, ax=axes[:, 2], label="$\\Delta$NMI", fraction=0.04, pad=0.02)
    fig.suptitle(f"Eg5 ADP-exit allosteric networks ({scheme}, phi backbone)\n"
                 f"NMI-weighted residue contact maps per LPATH cluster",
                 fontsize=13, y=0.995)
    out = f"{PLOTS}/grand_grid_{scheme}.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)

    # ---- per-cluster detail figures ----
    for c in range(N_CLUSTERS):
        if ("nomon", c) not in mats or ("wmon", c) not in mats:
            continue
        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        for k, cond in enumerate(("nomon", "wmon")):
            im = ax[k].imshow(mats[(cond, c)], cmap="magma", vmin=0, vmax=vmax,
                              extent=extent, interpolation="nearest")
            ax[k].set_title(f"{cond}  cluster {c}")
            ax[k].set_xlabel("residue"); ax[k].set_ylabel("residue")
            fig.colorbar(im, ax=ax[k], fraction=0.046, label="NMI")
        imd = ax[2].imshow(diffs[c], cmap="RdBu_r", norm=dnorm,
                           extent=extent, interpolation="nearest")
        ax[2].set_title(f"$\\Delta$ wmon$-$nomon  cluster {c}")
        ax[2].set_xlabel("residue"); ax[2].set_ylabel("residue")
        fig.colorbar(imd, ax=ax[2], fraction=0.046, label="$\\Delta$NMI")
        fig.tight_layout()
        out = f"{PLOTS}/cluster{c}_{scheme}.png"
        fig.savefig(out, dpi=170)
        plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    main()
