"""
Per-condition pdist landscapes with the shared-cluster representative path
traces overlaid, restricted to the successful source->target sub-path.

For each Eg5 WE run (nomon = no monastrol, wmon = + monastrol) we hexbin the
accessible ensemble (full west.h5) in a 2D reaction plane and color each hex by
the average of a 3rd progress coordinate, then draw the three shared-path-cluster
representatives on top.

pcoord layout (west.cfg, pcoord_ndim = 4):
    idx0  ADP-Mg RMSD (A)
    idx1  ADP+Mg / Eg5 interaction energy (kcal/mol)
    idx2  min contact distance ADP / Eg5 (A)
    idx3  min contact distance PO4 / Eg5 (A)

Two axis layouts are produced:
    A (rmsd_ene) : X = RMSD,   Y = int. energy,  color = min ADP-Eg5 dist
    B (ene_mind) : X = int. energy, Y = min dist, color = ADP-Mg RMSD

Traces: the highest-flux representative WITHIN each condition for each of the
three shared path-dendrogram clusters (extract_rep_paths.py --rep
weight-per-cond).  Unlike the earlier version, the trace is NOT wedap's full
walker-lineage trace back to iteration 1; instead we plot the pcoord values that
lpath stored for the successful pathway itself (rows with iter > 0 in
reassigned.pickle), i.e. exactly the source->target exit event.

Pathway columns (reassigned.pickle rows, verified):
    0 iter | 1 seg | 2 shared-state | 3 pcoord0 RMSD | 4 pcoord1 ene
    5 pcoord2 minADP | 6 pcoord3 minPO4 | 7 label | 8 frame | 9 weight

Outputs (into this dir), for each run and layout:
    pdist_rmsd_ene_<run>.pdf/.png
    pdist_ene_mind_<run>.pdf/.png
"""
import os
import pickle

import numpy as np
import matplotlib.pyplot as plt
import wedap

# --- config -----------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
LPATH = os.path.dirname(HERE)

RUNS = {
    "nomon": {"h5": os.path.join(LPATH, "nomon", "west.h5"),
              "pickle": os.path.join(LPATH, "nomon", "succ_traj", "reassigned.pickle"),
              "title": "no monastrol (WT Eg5)"},
    "wmon":  {"h5": os.path.join(LPATH, "wmon", "west.h5"),
              "pickle": os.path.join(LPATH, "wmon", "succ_traj", "reassigned.pickle"),
              "title": "+ monastrol"},
}
SUMMARY = os.path.join(LPATH, "shared_paths", "rep_paths_percond",
                       "rep_paths_summary.csv")

# distinct, colorblind-friendly per-cluster colors (match the WIPA poster set)
CL_COLORS = {1: "#377eb8", 2: "#d62728", 3: "#4daf4a"}

# coordinate roles: (wedap pcoord index, pathway column, axis label)
COORD = {
    "rmsd": (0, 3, r"ADP-Mg RMSD ($\mathrm{\AA}$)"),
    "ene":  (1, 4, r"Interaction energy (kcal/mol)"),
    "mind": (2, 5, r"min ADP-Eg5 distance ($\mathrm{\AA}$)"),
}
# axis layouts: name -> (x role, y role, color/z role)
LAYOUTS = {
    "rmsd_ene": ("rmsd", "ene", "mind"),
    "ene_mind": ("ene", "mind", "rmsd"),
}

STYLE = "/ihome/lchong/dty7/Apps/wedap/wedap/styles/default.mplstyle"


def load_reps():
    """Return {run: {cluster: pathway_ndarray}} for the per-cond reps."""
    reps = {r: {} for r in RUNS}
    with open(SUMMARY) as f:
        header = f.readline().strip().split(",")
        ci = {name: i for i, name in enumerate(header)}
        for line in f:
            row = line.strip().split(",")
            run = row[ci["origin"]]
            cl = int(row[ci["cluster"]])
            li = int(row[ci["local_idx"]])
            with open(RUNS[run]["pickle"], "rb") as pf:
                paths = pickle.load(pf)
            reps[run][cl] = np.asarray(paths[li], dtype=object)
    return reps


def subpath_xy(pathway, xcol, ycol):
    """(x, y) arrays of the successful source->target sub-path (iter > 0)."""
    a = pathway
    nz = a[a[:, 0].astype(float) > 0]
    x = nz[:, xcol].astype(float)
    y = nz[:, ycol].astype(float)
    return x, y


def plot(run, layout, reps):
    info = RUNS[run]
    xr, yr, zr = LAYOUTS[layout]
    xi, xcol, xlab = COORD[xr]
    yi, ycol, ylab = COORD[yr]
    zi, _, zlab = COORD[zr]

    wp = wedap.H5_Plot(
        h5=info["h5"], data_type="average", plot_mode="hexbin3d",
        Xname="pcoord", Xindex=xi,
        Yname="pcoord", Yindex=yi,
        Zname="pcoord", Zindex=zi,
        cmap="copper", cbar_label=zlab, hexbin_grid=120,
    )
    wp.plot()
    ax = wp.ax

    for cl in sorted(reps[run]):
        x, y = subpath_xy(reps[run][cl], xcol, ycol)
        c = CL_COLORS[cl]
        # black underlay + colored line, like wedap's own trace styling
        ax.plot(x, y, color="black", lw=2.5, zorder=3)
        ax.plot(x, y, color=c, lw=1.5, zorder=4)
        ax.scatter(x[0], y[0], marker="o", color=c, s=45, edgecolor="k", zorder=5)
        ax.scatter(x[-1], y[-1], marker="v", color=c, s=45, edgecolor="k", zorder=5)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(f"Eg5 ADP exit -- {info['title']}")
    handles = [plt.Line2D([0], [0], color=CL_COLORS[cl], lw=2,
                          label=f"path cluster {cl}") for cl in sorted(reps[run])]
    ax.legend(handles=handles, loc="best", fontsize=8, framealpha=0.85)

    plt.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(HERE, f"pdist_{layout}_{run}.{ext}")
        plt.savefig(out, dpi=300)
        print("wrote", out)
    plt.close("all")


def main():
    if os.path.exists(STYLE):
        try:
            plt.style.use(STYLE)
        except Exception as e:  # noqa: BLE001
            print("style skipped:", e)
    reps = load_reps()
    for run in RUNS:
        for layout in LAYOUTS:
            print(f"=== {run} / {layout} ===")
            plot(run, layout, reps)


if __name__ == "__main__":
    main()
