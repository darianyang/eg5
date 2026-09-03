"""
Per-condition pdist landscapes with the SHARED-BIN-GRID path-cluster
representative FULL-TRAJECTORY traces overlaid (max-weight/flux path per
cluster).

Difference from the source->target version: here each representative is drawn as
its complete WE walker lineage traced back to iteration 1 (wedap plot_trace),
i.e. the full trajectory from the basis state through the ADP-exit event -- not
just the successful source->target sub-path.

The shared pathway clustering on the fixed cfg-derived bin grid
(cluster_paths_shared_binned.py) separates the two conditions completely; at k=3

    cluster 1 -> 100% nomon (WT)       cluster 2 -> 100% wmon (+ monastrol)
                                       cluster 3 -> 100% wmon (+ monastrol)

Per request, the WT cluster (1) is split one level deeper -- into its two next
dendrogram branches (1a, 1b) -- so the WT panel shows two distinct exit paths.
The wmon panel keeps clusters 2 and 3.  Each rep is the highest-WE-weight (flux)
pathway within its (sub)cluster.

Outputs (into this dir), for each run and layout:
    pdist_rmsd_ene_<run>_binned.pdf/.png
    pdist_ene_mind_<run>_binned.pdf/.png
"""
import os
import csv
import pickle

import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
import wedap

# --- config -----------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
LPATH = os.path.dirname(HERE)

RUNS = {
    "nomon": {"h5": os.path.join(LPATH, "nomon", "west.h5"),
              "pickle": os.path.join(LPATH, "nomon", "succ_traj",
                                     "reassigned_binned.pickle"),
              "title": "no monastrol (WT Eg5)"},
    "wmon":  {"h5": os.path.join(LPATH, "wmon", "west.h5"),
              "pickle": os.path.join(LPATH, "wmon", "succ_traj",
                                     "reassigned_binned.pickle"),
              "title": "+ monastrol"},
}
DISTMAT = os.path.join(LPATH, "shared_paths_binned", "distmat_shared.npy")
ASSIGNMENTS = os.path.join(LPATH, "shared_paths_binned", "path_assignments.csv")
BASE_K = 3       # matches cluster_paths_shared_binned.py auto-pick

# coordinate roles: (wedap pcoord index, axis label)
COORD = {
    "rmsd": (0, r"ADP-Mg RMSD ($\mathrm{\AA}$)"),
    "ene":  (1, r"Interaction energy (kcal/mol)"),
    "mind": (2, r"min ADP-Eg5 distance ($\mathrm{\AA}$)"),
}
LAYOUTS = {
    "rmsd_ene": ("rmsd", "ene", "mind"),
    "ene_mind": ("ene", "mind", "rmsd"),
}

STYLE = "/ihome/lchong/dty7/Apps/wedap/wedap/styles/default.mplstyle"


def load_assignments():
    """Return pooled-ordered arrays (origin, weight) and n_nomon offset."""
    origin, weight = [], []
    with open(ASSIGNMENTS) as f:
        for r in csv.DictReader(f):
            origin.append(r["origin"])
            weight.append(float(r["weight"]))
    origin = np.array(origin)
    weight = np.array(weight)
    n_nomon = int((origin == "nomon").sum())
    return origin, weight, n_nomon


def final_iter_seg(pathway):
    """Final (iter, seg) of a pathway = last non-padding (iter > 0) frame."""
    a = np.asarray(pathway, dtype=float)
    real = a[a[:, 0] > 0]
    return int(real[-1, 0]), int(real[-1, 1])


def build_reps():
    """Return {run: [(label, (iter,seg), color), ...]} of full-traj reps.

    WT (all-nomon k=3 cluster) is split into its two next dendrogram branches;
    the two all-wmon k=3 clusters are kept as clusters 2 and 3.
    """
    origin, weight, n_nomon = load_assignments()
    dm = np.load(DISTMAT)
    z = sch.linkage(squareform(dm, checks=False), method="ward")
    base = sch.fcluster(z, t=BASE_K, criterion="maxclust")

    # identify the all-nomon (WT) base cluster and the wmon clusters
    wt_label = None
    wmon_labels = []
    for cl in sorted(np.unique(base)):
        members = np.where(base == cl)[0]
        origins = set(origin[members])
        if origins == {"nomon"}:
            wt_label = cl
        else:
            wmon_labels.append(cl)
    if wt_label is None:
        raise RuntimeError("no all-nomon base cluster found; check clustering")

    # split the WT cluster along its next branch: increase k until its members
    # occupy exactly two labels (that is the WT subtree's top internal merge)
    c1 = np.where(base == wt_label)[0]
    split = None
    for k in range(BASE_K + 1, len(z) + 1):
        lab = sch.fcluster(z, t=k, criterion="maxclust")
        sub = lab[c1]
        if len(np.unique(sub)) == 2:
            split = sub
            break
    if split is None:
        raise RuntimeError("WT cluster never split into two branches")

    def load_pickles():
        return {run: pickle.load(open(RUNS[run]["pickle"], "rb")) for run in RUNS}

    pk = load_pickles()

    def rep_of(pooled_members):
        """(run, local_idx, (iter,seg), weight) of the max-weight member."""
        g = pooled_members[int(np.argmax(weight[pooled_members]))]
        run = "nomon" if g < n_nomon else "wmon"
        local = g if run == "nomon" else g - n_nomon
        return run, local, final_iter_seg(pk[run][local]), float(weight[g])

    reps = {run: [] for run in RUNS}

    # WT sub-paths 1a / 1b (ordered by descending flux so 1a is the dominant one)
    sub_ids = np.unique(split)
    sub_groups = [c1[split == s] for s in sub_ids]
    sub_groups.sort(key=lambda m: -weight[m].max())
    wt_colors = ["#377eb8", "#984ea3"]          # blue, purple
    for tag, members, color in zip(["1a", "1b"], sub_groups, wt_colors):
        run, local, itseg, w = rep_of(members)
        reps[run].append((tag, itseg, color))
        print(f"WT path {tag}: origin={run} local={local} iter,seg={itseg} "
              f"weight={w:.3e} (n={len(members)})")

    # wmon clusters 2 / 3 (report in their original label order)
    wmon_colors = {wmon_labels[0]: "#d62728", wmon_labels[-1]: "#4daf4a"}
    for i, cl in enumerate(sorted(wmon_labels), start=2):
        members = np.where(base == cl)[0]
        run, local, itseg, w = rep_of(members)
        reps[run].append((str(i), itseg, wmon_colors[cl]))
        print(f"path cluster {i}: origin={run} local={local} iter,seg={itseg} "
              f"weight={w:.3e} (n={len(members)})")

    return reps


def plot(run, layout, reps):
    info = RUNS[run]
    xr, yr, zr = LAYOUTS[layout]
    xi, xlab = COORD[xr]
    yi, ylab = COORD[yr]
    zi, zlab = COORD[zr]

    wp = wedap.H5_Plot(
        h5=info["h5"], data_type="average", plot_mode="hexbin3d",
        Xname="pcoord", Xindex=xi,
        Yname="pcoord", Yindex=yi,
        Zname="pcoord", Zindex=zi,
        cmap="copper", cbar_label=zlab, hexbin_grid=120,
    )
    wp.plot()
    ax = wp.ax

    handles = []
    for tag, itseg, color in reps[run]:
        # full walker lineage back to iteration 1, on this layout's axes
        wp.plot_trace(itseg, ax=ax, color=color, linewidth=1.5,
                      mark_points=True, mp_size=45)
        label = (f"WT path {tag}" if tag in ("1a", "1b")
                 else f"path cluster {tag}")
        handles.append(plt.Line2D([0], [0], color=color, lw=2, label=label))

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"Eg5 ADP exit (full traj) -- {info['title']}")
    if handles:
        ax.legend(handles=handles, loc="best", fontsize=8, framealpha=0.85)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(HERE, f"pdist_{layout}_{run}_binned.{ext}")
        plt.savefig(out, dpi=300)
        print("wrote", out)
    plt.close("all")


def main():
    if os.path.exists(STYLE):
        try:
            plt.style.use(STYLE)
        except Exception as e:  # noqa: BLE001
            print("style skipped:", e)
    reps = build_reps()
    for run in RUNS:
        for layout in LAYOUTS:
            print(f"=== {run} / {layout} ===")
            plot(run, layout, reps)


if __name__ == "__main__":
    main()
