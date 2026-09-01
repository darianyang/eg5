"""
Shared-space clustering for the two Eg5 ADP-unbinding WE simulations
(nomon = WT/no monastrol, wmon = with monastrol).

Both simulations share the same 4D progress coordinate:
    dim0 = ADP-Mg RMSD (A)
    dim1 = interaction energy ADP+Mg / Eg5 (kcal/mol)
    dim2 = ADP min contact distance (A)
    dim3 = PO4 min contact distance (A)

We build ONE clustering model on the combined data from both runs, so the
resulting cluster IDs mean the same thing in both conditions and the LPATH
pathways can be compared directly.  The cluster labels are written back into
each run's west.h5 as an 'auxdata/labels' dataset (matching the LPATH WE
example) for use in the extract/match steps.
"""
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tqdm.auto import trange
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn import neighbors

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
RUNS = {
    "nomon": "nomon/west.h5",
    "wmon":  "wmon/west.h5",
}
FEATURES = ["ADP-Mg RMSD ($\\AA$)", "Int. Ene. (kcal/mol)",
            "ADP min dist ($\\AA$)", "PO4 min dist ($\\AA$)"]
N_CLUSTERS   = 6       # number of shared clusters
N_TRAIN      = 5000    # size of stratified training subset for agglomerative fit
LINKAGE      = "ward"
RANDOM_SEED  = 42
rng = np.random.default_rng(RANDOM_SEED)


def load_all_frames(path):
    """Return (n_frames, 4) array of every pcoord frame in the run."""
    rows = []
    with h5py.File(path, "r") as f:
        iters = sorted(f["iterations"].keys())
        for k in iters:
            pc = f["iterations/" + k + "/pcoord"][:]      # (segs, len, 4)
            rows.append(pc.reshape(-1, 4))
    return np.concatenate(rows).astype(np.float64)


# ----------------------------------------------------------------------------
# 1. Load combined data from both runs
# ----------------------------------------------------------------------------
print("Loading pcoord data from both runs ...")
data = {}
for name, path in RUNS.items():
    data[name] = load_all_frames(path)
    print(f"  {name}: {data[name].shape[0]} frames")

combined = np.concatenate([data[n] for n in RUNS], axis=0)
origin = np.concatenate([np.full(data[n].shape[0], n) for n in RUNS])

# ----------------------------------------------------------------------------
# 2. Standardize (shared scaler) and fit clustering on a stratified subset
# ----------------------------------------------------------------------------
print("Standardizing features (shared scaler) ...")
scaler = StandardScaler().fit(combined)
combined_s = scaler.transform(combined)

# stratified training subset: equal draw from each run
train_idx = []
for name in RUNS:
    idx = np.where(origin == name)[0]
    n_take = min(N_TRAIN // len(RUNS), idx.size)
    train_idx.append(rng.choice(idx, size=n_take, replace=False))
train_idx = np.concatenate(train_idx)
X_train = combined_s[train_idx]

print(f"Agglomerative clustering ({LINKAGE}, k={N_CLUSTERS}) on "
      f"{X_train.shape[0]} training points ...")
ag = AgglomerativeClustering(n_clusters=N_CLUSTERS, linkage=LINKAGE).fit(X_train)
labels_train = ag.labels_

# ----------------------------------------------------------------------------
# 3. Propagate labels to every frame via KNN (shared model)
# ----------------------------------------------------------------------------
print("Propagating cluster labels with KNN ...")
knn = neighbors.KNeighborsClassifier(n_neighbors=5)
knn.fit(X_train, labels_train)
labels_all = knn.predict(combined_s)

# report cluster populations and centroids (in original feature units)
uniq = np.unique(labels_all)
centroids = []
print("\nShared cluster summary (centroids in original units):")
print(f"{'cid':>3} {'n_total':>9} {'n_nomon':>8} {'n_wmon':>8}  "
      + "  ".join(f"{n.split()[0]:>8}" for n in FEATURES))
for c in uniq:
    m = labels_all == c
    cen = combined[m].mean(axis=0)
    centroids.append(cen)
    n_no = np.sum(m & (origin == "nomon"))
    n_wm = np.sum(m & (origin == "wmon"))
    print(f"{c:>3} {m.sum():>9} {n_no:>8} {n_wm:>8}  "
          + "  ".join(f"{v:>8.2f}" for v in cen))
centroids = np.asarray(centroids)
np.save("centroids.npy", centroids)

# ----------------------------------------------------------------------------
# 4. Write labels back into each run's west.h5 as auxdata/labels
# ----------------------------------------------------------------------------
print("\nWriting auxdata/labels into each west.h5 ...")
for name, path in RUNS.items():
    with h5py.File(path, "a") as f:
        iters = sorted(f["iterations"].keys())
        for k in iters:
            pc = f["iterations/" + k + "/pcoord"][:]       # (segs, len, 4)
            n_seg, plen, _ = pc.shape
            Xs = scaler.transform(pc.reshape(-1, 4))
            lab = knn.predict(Xs).reshape(n_seg, plen, 1).astype(np.int32)
            ds = "iterations/" + k + "/auxdata/labels"
            if ds in f:
                f[ds][...] = lab
            else:
                f.create_dataset(ds, data=lab)
    print(f"  {name}: labels written for {len(iters)} iterations")

# ----------------------------------------------------------------------------
# 5. Plots: 2D projections colored by shared cluster
# ----------------------------------------------------------------------------
print("Making cluster plots ...")
cmap = matplotlib.colormaps.get_cmap("tab10")
colors = [cmap(i) for i in range(10)]

def scatter_proj(ax, X, labs, ix, iy, title):
    for c in uniq:
        m = labs == c
        ax.scatter(X[m, ix], X[m, iy], s=2, color=colors[c % 10],
                   label=f"{c}", rasterized=True)
        cen = X[m].mean(axis=0)
        ax.scatter(cen[ix], cen[iy], color="black", s=40, zorder=3,
                   marker="X", edgecolors="white", linewidths=0.5)
    ax.set_xlabel(FEATURES[ix])
    ax.set_ylabel(FEATURES[iy])
    ax.set_title(title)

# projection pairs: RMSD vs ADPdist, and IntEne vs ADPdist
pairs = [(0, 2), (1, 2), (0, 1)]
for name in list(RUNS) + ["combined"]:
    if name == "combined":
        X, labs = combined, labels_all
    else:
        m = origin == name
        X, labs = combined[m], labels_all[m]
    fig, axes = plt.subplots(1, len(pairs), figsize=(5 * len(pairs), 4.2))
    for ax, (ix, iy) in zip(axes, pairs):
        scatter_proj(ax, X, labs, ix, iy, f"{name}")
    axes[-1].legend(title="cluster", markerscale=4, fontsize=8,
                    loc="upper right")
    fig.tight_layout()
    fig.savefig(f"plots/clusters_{name}.png", dpi=200)
    plt.close(fig)
    print(f"  wrote plots/clusters_{name}.png")

print("\nDone. centroids.npy and plots/clusters_*.png written; "
      "auxdata/labels added to both west.h5 files.")
