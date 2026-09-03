"""Render the Shortest Path Map on the 3D structure, per (condition, cluster).

Two outputs from the SPM caches written by spm_analysis.py (spm/<cond>/clusterN.npz):

  1. ChimeraX ``.bild`` + ``.cxc`` (open with ``chimerax viz/spm/clusterN.cxc``):
     top-betweenness edges as cylinders and top nodes as spheres on the shared
     reference structure -- nomon, wmon, and their drug difference as separate
     models.  Radius scales with SPM betweenness; the difference map colours
     edges red where monastrol gains path traffic and blue where it loses it.

  2. A matplotlib 2D projection PNG per cluster (PCA of the reference Ca cloud),
     three panels nomon / wmon / drug-difference, so the map is viewable without
     ChimeraX -- analogous to the paper's SPM structure figures.

    python scripts/make_spm_viz.py [--ref data/ref_nomon.pdb] [--cache-dir spm]
        [--out-dir viz/spm] [--data-dir data] [--top-edges 120] [--top-nodes 25]
"""
import os
import sys
import argparse
import numpy as np
import h5py
import mdtraj as md

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402
from df_analysis import EXIT_ORDER, REGIONS  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

CONDS = ("nomon", "wmon")


def ref_ca(ref_pdb, res_ids):
    """Ca coordinates (Angstrom) ordered to match the SPM node index (res_ids)."""
    t = md.load(ref_pdb)
    ca = t.top.select("name CA")
    xyz = {t.top.atom(a).residue.resSeq: t.xyz[0, a] * 10.0 for a in ca}
    return np.array([xyz[int(r)] for r in res_ids])          # (n_res, 3)


def top_edges(M, k):
    """Indices (i<j) and values of the k largest-|value| upper-triangle edges."""
    R = M.shape[0]
    iu, iv = np.triu_indices(R, k=1)
    v = M[iu, iv]
    sel = np.argsort(np.abs(v))[::-1][:k]
    return iu[sel], iv[sel], v[sel]


def rgb(c):
    return f"{c[0]:.3f} {c[1]:.3f} {c[2]:.3f}"


def write_bild(path, coords, ie, je, ve, nodes, node_val, diverging):
    """One .bild: edges as cylinders, nodes as spheres.  For diverging (drug)
    maps red = gained / blue = lost; otherwise colour by betweenness (plasma)."""
    vmax = np.abs(ve).max() or 1.0
    plasma = plt.get_cmap("plasma")
    with open(path, "w") as f:
        for i, j, v in zip(ie, je, ve):
            if diverging:
                col = (0.85, 0.15, 0.15) if v > 0 else (0.15, 0.3, 0.85)
            else:
                col = plasma(v / vmax)[:3]
            r = 0.15 + 0.7 * (abs(v) / vmax)
            f.write(f".color {rgb(col)}\n")
            x1, y1, z1 = coords[i]; x2, y2, z2 = coords[j]
            f.write(f".cylinder {x1:.3f} {y1:.3f} {z1:.3f} "
                    f"{x2:.3f} {y2:.3f} {z2:.3f} {r:.3f}\n")
        nmax = node_val[nodes].max() or 1.0
        f.write(".color 0.1 0.1 0.5\n")
        for n in nodes:
            x, y, z = coords[n]
            r = 0.5 + 1.6 * (node_val[n] / nmax)
            f.write(f".sphere {x:.3f} {y:.3f} {z:.3f} {r:.3f}\n")


def write_cxc(path, ref_abspath, bilds):
    with open(path, "w") as f:
        f.write(f"open {ref_abspath}\n")
        f.write("hide atoms\nshow cartoon\ncolor light gray\n")
        f.write("transparency 55 target c\n")
        for name, b in bilds:
            f.write(f"open {os.path.abspath(b)}\n")
        f.write("set bgColor white\nlighting soft\n")
        f.write("# models: 2=nomon SPM, 3=wmon SPM, 4=drug difference "
                "(red gained / blue lost)\n")


def project(coords):
    """PCA of the Ca cloud -> 2D coordinates for a flat structure view."""
    X = coords - coords.mean(0)
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    return X @ Vt[:2].T                                       # (n_res, 2)


def panel(ax, P, res_ids, ie, je, ve, node_val, title, diverging):
    ax.plot(P[:, 0], P[:, 1], color="0.8", lw=0.8, zorder=1)  # backbone trace
    vmax = np.abs(ve).max() or 1.0
    segs, cols, lws = [], [], []
    for i, j, v in zip(ie, je, ve):
        segs.append([P[i], P[j]])
        if diverging:
            cols.append((0.85, 0.15, 0.15) if v > 0 else (0.15, 0.3, 0.85))
        else:
            cols.append(plt.get_cmap("plasma")(v / vmax))
        lws.append(0.4 + 3.0 * (abs(v) / vmax))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws, zorder=2))
    nmax = node_val.max() or 1.0
    sizes = 6 + 120 * (node_val / nmax)
    ax.scatter(P[:, 0], P[:, 1], s=sizes, c="0.35", zorder=3, edgecolors="none")
    # label a few landmark regions at their Ca centroid
    for name, (lo, hi) in REGIONS.items():
        m = (res_ids >= lo) & (res_ids <= hi)
        if m.any():
            c = P[m].mean(0)
            ax.text(c[0], c[1], name, fontsize=6, color="teal", ha="center",
                    zorder=4)
    ax.set_title(title, fontsize=9)
    ax.set_aspect("equal"); ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="data/ref_nomon.pdb")
    ap.add_argument("--cache-dir", default="spm")
    ap.add_argument("--out-dir", default="viz/spm")
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--top-edges", type=int, default=120)
    ap.add_argument("--top-nodes", type=int, default=25)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    plots_dir = f"{a.out_dir}/projections"
    os.makedirs(plots_dir, exist_ok=True)

    with h5py.File(f"{a.data_dir}/ca_coords_nomon.h5", "r") as f:
        res_ids = f["res_ids"][:]
    coords = ref_ca(a.ref, res_ids)
    P = project(coords)
    ref_abs = os.path.abspath(a.ref)

    def load(cond, c):
        p = f"{a.cache_dir}/{cond}/cluster{c}.npz"
        if not os.path.exists(p):
            return None
        d = np.load(p)
        return d["E"], d["node"]

    order = [c for c in EXIT_ORDER
             if load("nomon", c) is not None or load("wmon", c) is not None]

    for c in order:
        got = {cond: load(cond, c) for cond in CONDS}
        bilds = []
        # per-condition bild + projection panels
        panels = []
        for cond in CONDS:
            if got[cond] is None:
                continue
            E, node = got[cond]
            ie, je, ve = top_edges(E, a.top_edges)
            nodes = np.argsort(node)[::-1][:a.top_nodes]
            cond_dir = f"{a.out_dir}/{cond}"
            os.makedirs(cond_dir, exist_ok=True)
            b = f"{cond_dir}/cluster{c}.bild"
            write_bild(b, coords, ie, je, ve, nodes, node, diverging=False)
            bilds.append((cond, b))
            panels.append((cond, ie, je, ve, node, False))
        # drug-difference bild + panel
        if got["nomon"] is not None and got["wmon"] is not None:
            D = got["wmon"][0] - got["nomon"][0]
            dnode = got["wmon"][1] - got["nomon"][1]
            ie, je, ve = top_edges(D, a.top_edges)
            nodes = np.argsort(np.abs(dnode))[::-1][:a.top_nodes]
            ddir = f"{a.out_dir}/delta"
            os.makedirs(ddir, exist_ok=True)
            b = f"{ddir}/cluster{c}.bild"
            write_bild(b, coords, ie, je, ve, nodes, np.abs(dnode), diverging=True)
            bilds.append(("delta", b))
            panels.append(("wmon-nomon", ie, je, ve, np.abs(dnode), True))

        write_cxc(f"{a.out_dir}/cluster{c}.cxc", ref_abs, bilds)

        # projection figure
        fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 5.2),
                                 squeeze=False)
        for ax, (tag, ie, je, ve, nv, div) in zip(axes[0], panels):
            ttl = (f"$\\Delta$SPM {tag} (red gained / blue lost)" if div
                   else f"{tag} SPM  cluster {c}")
            panel(ax, P, res_ids, ie, je, ve, nv, ttl, div)
        fig.suptitle(f"Shortest Path Map on structure -- cluster {c} "
                     f"(top {a.top_edges} edges; PCA projection)", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        out = f"{plots_dir}/cluster{c}_spm_structure.png"
        fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
        print("wrote", out, "and", f"{a.out_dir}/cluster{c}.cxc")


if __name__ == "__main__":
    main()
