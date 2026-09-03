"""Distance Fluctuation (DF) and Delta-DF allosteric analysis per LPATH cluster.

Reimplements the DF method of Rodriguez-Santos et al. (JCIM 2026) for the two
Eg5 ADP-unbinding WE simulations, resolved by the 6 shared LPATH clusters.

    DF_ij = < ( d_ij - <d_ij> )^2 >

where d_ij is the Ca-Ca distance of residues i,j in a frame and <.> averages
over that cluster's frame ensemble (uniform, unweighted -- see extract_ca_coords).
DF is rotation/translation invariant, so no alignment is needed.  A LOW DF_ij
means residues i,j move in a coordinated (allosteric) way.

Delta-DF compares two ensembles (paper convention: red/positive = LOSS of
coordination, black/negative = GAIN):
  * drug  : wmon - nomon within a cluster (monastrol's footprint at that stage),
  * stage : DF(later) - DF(earlier) between adjacent clusters along the ADP-exit
            coordinate, within a condition.

ADP-exit ordering (from lpath/centroids.npy, bound -> unbound):
    2 -> 5 -> 0 -> 3 -> 4 -> 1

    python scripts/df_analysis.py [--data-dir data] [--out-dir plots/df]
        [--cache-dir df] [--recompute]

Caches DF matrices to ``<cache-dir>/<cond>/cluster<c>.npy`` and writes figures
to ``<out-dir>/``.
"""
import os
import sys
import argparse
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

CONDS = ("nomon", "wmon")
# ADP-exit order (bound -> unbound), inferred from lpath/centroids.npy
EXIT_ORDER = [2, 5, 0, 3, 4, 1]
# motor-domain regions in this construct's AMBER numbering.  Landmark residues
# verified against the topology (Trp127, Tyr211, P-loop, SWI Ser233/Arg234,
# Glu270); helix spans (α3/α4/α6) placed by DSSP on data/ref_nomon.pdb.
REGIONS = {
    "N-term": (1, 16), "α0": (30, 36), "P-loop": (105, 112),
    "L5": (116, 134), "L8": (166, 206), "α3": (209, 226),
    "SWI": (230, 237), "SWII": (262, 270), "L11": (271, 281),
    "α4": (282, 303), "α6": (340, 356),
}


def df_matrix(ca, chunk=128):
    """Distance-fluctuation matrix (n_res, n_res) for a frame ensemble.

    ``ca`` is (n_frame, n_res, 3).  DF_ij = E[d_ij^2] - E[d_ij]^2, accumulated
    in frame chunks.  Per-frame squared distances come from the Gram matrix
    (d^2 = |r_i|^2 + |r_j|^2 - 2 r_i.r_j), which is BLAS-fast and avoids the
    (m, R, R, 3) coordinate-difference tensor.
    """
    n, R, _ = ca.shape
    sumd = np.zeros((R, R), np.float64)
    sumd2 = np.zeros((R, R), np.float64)
    for s in range(0, n, chunk):
        x = ca[s:s + chunk].astype(np.float64) * 10.0   # (m, R, 3) nm -> Angstrom
        nrm = np.einsum("mrc,mrc->mr", x, x)            # (m, R)  |r|^2
        gram = x @ x.transpose(0, 2, 1)                 # (m, R, R)  r_i . r_j
        d2 = nrm[:, :, None] + nrm[:, None, :] - 2.0 * gram
        d = np.sqrt(np.maximum(d2, 0.0))                # (m, R, R)  distances
        sumd += d.sum(0)
        sumd2 += d2.sum(0)                              # E[d^2] uses d2 directly
    mean = sumd / n
    df = sumd2 / n - mean * mean
    np.fill_diagonal(df, 0.0)
    return df


def load_cluster_ca(h5, cluster):
    with h5py.File(h5, "r") as f:
        sel = f["cluster"][:] == cluster
        ca = f["ca"][:][sel]
        res = f["res_ids"][:]
    return ca, res


def compute_all(data_dir, cache_dir, recompute):
    """Return {(cond, cluster): DF matrix}, res_ids, {(cond,cluster): n_frame}."""
    mats, counts, res_ids = {}, {}, None
    for cond in CONDS:
        h5 = f"{data_dir}/ca_coords_{cond}.h5"
        for c in range(N_CLUSTERS):
            cache = f"{cache_dir}/{cond}/cluster{c}.npy"
            ca, res = load_cluster_ca(h5, c)
            res_ids = res
            counts[(cond, c)] = len(ca)
            if os.path.exists(cache) and not recompute:
                mats[(cond, c)] = np.load(cache)
                continue
            if len(ca) < 2:
                continue
            df = df_matrix(ca)
            os.makedirs(os.path.dirname(cache), exist_ok=True)
            np.save(cache, df)
            mats[(cond, c)] = df
            print(f"  DF {cond} c{c}: n={len(ca)} -> {cache}", flush=True)
    return mats, res_ids, counts


def add_region_ticks(ax, res_ids):
    """Draw faint region boundaries + labels on both axes."""
    lo, hi = int(res_ids[0]), int(res_ids[-1])

    def to_idx(r):
        return np.searchsorted(res_ids, r)
    for name, (a, b) in REGIONS.items():
        if b < lo or a > hi:
            continue
        for r in (a, b):
            ax.axhline(to_idx(r), color="c", lw=0.3, alpha=0.4)
            ax.axvline(to_idx(r), color="c", lw=0.3, alpha=0.4)
        mid = to_idx((a + b) // 2)
        ax.text(mid, -6, name, ha="center", va="bottom", fontsize=6, color="teal")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default="plots/df")
    ap.add_argument("--cache-dir", default="df")
    ap.add_argument("--recompute", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    mats, res_ids, counts = compute_all(a.data_dir, a.cache_dir, a.recompute)
    if not mats:
        print("no DF matrices computed (missing ca_coords_*.h5?)")
        return
    n = len(res_ids)
    extent = [res_ids[0], res_ids[-1], res_ids[-1], res_ids[0]]

    # shared DF colour scale (99th pct over all matrices)
    allv = np.concatenate([M[M > 0].ravel() for M in mats.values()])
    dfmax = np.percentile(allv, 99)

    # ---- grand grid: rows = clusters in exit order, cols = [nomon, wmon, drug d] ----
    order = [c for c in EXIT_ORDER if (("nomon", c) in mats or ("wmon", c) in mats)]
    fig, axes = plt.subplots(len(order), 3, figsize=(14, 4.4 * len(order)),
                             squeeze=False)
    drug = {}
    for row, c in enumerate(order):
        for col, cond in enumerate(CONDS):
            ax = axes[row][col]
            M = mats.get((cond, c))
            if M is None:
                ax.set_visible(False); continue
            im = ax.imshow(M, cmap="viridis_r", vmin=0, vmax=dfmax,
                           extent=extent, interpolation="nearest")
            ax.set_title(f"{cond}  cluster {c}  (n={counts[(cond,c)]})", fontsize=9)
            add_region_ticks(ax, res_ids)
        axd = axes[row][2]
        if ("nomon", c) in mats and ("wmon", c) in mats:
            D = mats[("wmon", c)] - mats[("nomon", c)]
            drug[c] = D
            dmax = np.percentile(np.abs(D), 99) or 1.0
            imd = axd.imshow(D, cmap="RdGy_r", norm=TwoSlopeNorm(0, -dmax, dmax),
                             extent=extent, interpolation="nearest")
            axd.set_title(f"$\\Delta$DF drug (wmon$-$nomon)  c{c}", fontsize=9)
            add_region_ticks(axd, res_ids)
            fig.colorbar(imd, ax=axd, fraction=0.046, label="$\\Delta$DF")
        else:
            axd.set_visible(False)
    fig.colorbar(im, ax=axes[:, :2], label="DF ($\\AA^2$)", fraction=0.02, pad=0.01)
    fig.suptitle("Eg5 ADP-exit Distance Fluctuation per LPATH cluster\n"
                 "(rows = exit order 2→5→0→3→4→1; "
                 "red = coordination lost with monastrol)", y=0.997, fontsize=12)
    out = f"{a.out_dir}/df_grand_grid.png"
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)

    # ---- stage Delta-DF along the exit coordinate, per condition ----
    pairs = list(zip(EXIT_ORDER[:-1], EXIT_ORDER[1:]))
    fig, axes = plt.subplots(len(CONDS), len(pairs),
                             figsize=(3.4 * len(pairs), 3.4 * len(CONDS)),
                             squeeze=False)
    for r, cond in enumerate(CONDS):
        for k, (a0, a1) in enumerate(pairs):
            ax = axes[r][k]
            if (cond, a0) not in mats or (cond, a1) not in mats:
                ax.set_visible(False); continue
            D = mats[(cond, a1)] - mats[(cond, a0)]
            dmax = np.percentile(np.abs(D), 99) or 1.0
            im = ax.imshow(D, cmap="RdGy_r", norm=TwoSlopeNorm(0, -dmax, dmax),
                           extent=extent, interpolation="nearest")
            ax.set_title(f"{cond}: c{a0}$\\to$c{a1}", fontsize=9)
            if k == 0:
                ax.set_ylabel(cond)
    fig.suptitle("Stage $\\Delta$DF along ADP exit (DF$_{later}-$DF$_{earlier}$)  "
                 "red = coordination lost as ADP leaves", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = f"{a.out_dir}/df_stage_delta.png"
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)

    # ---- per-residue drug footprint: mean |dDF| involvement per residue ----
    if drug:
        fig, ax = plt.subplots(figsize=(12, 4))
        for c in order:
            if c in drug:
                involve = np.abs(drug[c]).mean(1)   # mean over partners
                ax.plot(res_ids, involve, lw=1, label=f"cluster {c}")
        ax.set_xlabel("residue (AMBER)")
        ax.set_ylabel("mean |$\\Delta$DF| (drug)")
        ax.set_title("Per-residue monastrol DF footprint across ADP-exit stages")
        for name, (lo, hi) in REGIONS.items():
            ax.axvspan(lo, hi, color="orange", alpha=0.08)
            ax.text((lo + hi) / 2, ax.get_ylim()[1], name, fontsize=6,
                    ha="center", va="top", color="darkorange")
        ax.legend(fontsize=8, ncol=3)
        fig.tight_layout()
        out = f"{a.out_dir}/df_drug_footprint.png"
        fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    main()
