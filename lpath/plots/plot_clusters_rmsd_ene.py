"""
Plot the shared-space clusters on the ADP-Mg RMSD (x) vs interaction energy (y)
plane -- the same axes used for the eg5_poster WE probability distributions.

Reads pcoord (dim0 = RMSD, dim1 = Int. Ene.) and the shared cluster labels
(auxdata/labels, written by ../cluster_shared.py) directly from each west.h5.
Read-only, so it is safe to run alongside the LPATH pipeline.
"""
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

RUNS = {
    "nomon": "../nomon/west.h5",
    "wmon":  "../wmon/west.h5",
}
# axis ranges chosen to match the eg5_poster RMSD vs Int. Ene. distributions
XLIM = (0, 20)      # ADP-Mg RMSD (A)
YLIM = (-250, 400)  # interaction energy (kcal/mol)


def load(path):
    rmsd, ene, lab = [], [], []
    with h5py.File(path, "r") as f:
        for k in sorted(f["iterations"].keys()):
            pc = f["iterations/" + k + "/pcoord"][:]              # (segs,len,4)
            lb = f["iterations/" + k + "/auxdata/labels"][:]      # (segs,len,1)
            rmsd.append(pc[:, :, 0].ravel())
            ene.append(pc[:, :, 1].ravel())
            lab.append(lb[:, :, 0].ravel())
    return (np.concatenate(rmsd), np.concatenate(ene),
            np.concatenate(lab).astype(int))


data = {name: load(path) for name, path in RUNS.items()}
uniq = np.unique(np.concatenate([d[2] for d in data.values()]))
cmap = matplotlib.colormaps.get_cmap("tab10")
colors = {c: cmap(c % 10) for c in uniq}

# ---- one panel per run + a combined panel --------------------------------
panels = list(RUNS) + ["combined"]
fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 4.6),
                         sharex=True, sharey=True)

for ax, name in zip(axes, panels):
    if name == "combined":
        rmsd = np.concatenate([data[n][0] for n in RUNS])
        ene = np.concatenate([data[n][1] for n in RUNS])
        lab = np.concatenate([data[n][2] for n in RUNS])
    else:
        rmsd, ene, lab = data[name]
    for c in uniq:
        m = lab == c
        ax.scatter(rmsd[m], ene[m], s=2, color=colors[c],
                   label=f"{c}", rasterized=True)
        ax.scatter(rmsd[m].mean(), ene[m].mean(), color="black", s=55,
                   marker="X", edgecolors="white", linewidths=0.7, zorder=3)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xlabel(r"ADP-Mg RMSD ($\AA$)")
    ax.set_title({"nomon": "no monastrol (WT)",
                  "wmon": "+ monastrol",
                  "combined": "combined"}[name])

axes[0].set_ylabel("interaction energy (kcal/mol)")
axes[-1].legend(title="cluster", markerscale=4, fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig("clusters_rmsd_ene.png", dpi=250)
fig.savefig("clusters_rmsd_ene.pdf")
print("wrote plots/clusters_rmsd_ene.png and .pdf")
