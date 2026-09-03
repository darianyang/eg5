"""
Build the SHARED, fixed bin grid used to discretize the LPATH intermediate
states for both Eg5 ADP-unbinding WE runs (nomon = WT, wmon = + monastrol).

Motivation
----------
The earlier LPATH run used 6 agglomerative shared-space clusters as the match
states (see cluster_shared.py / reassign_custom.py).  This replaces those 6
clusters with a much finer, *fixed* rectilinear grid whose edges are taken
directly from the WE binning structure in west.cfg, so the state space is:

  * physically meaningful (the same axes the WE actually sampled along),
  * identical between nomon and wmon (built once, here), so the resulting
    pathway dendrogram is directly comparable, exactly as with the shared
    6-cluster space.

Grid (3 coords; PO4 min-dist is dropped -- it is unconstrained in the
FULLY_UNBOUND2 source/sink definition anyway):

  * dim0  ADP-Mg RMSD (A):
        west.cfg base boundaries 0, 7.5, 10.5, inf.
        MAB mapper 1 -> 5 equidistant bins in [0, 7.5]
        MAB mapper 2 -> 5 equidistant bins in [7.5, 10.5]
        RMSD >10.5 -> 1 bin
        => 11 bins
  * dim1  ADP+Mg / Eg5 interaction energy (kcal/mol):
        west.cfg base boundary at 10 (['-inf', 10, 'inf']).  The refining MAB
        mappers each request 5 bins but over infinite ranges, so we make the
        grid finite and shared by splitting at 10 and placing 5 equidistant
        bins on each side across the pooled observed range => 10 bins.
  * dim2  ADP min contact distance (A):
        west.cfg base boundaries 0, 6, inf.
        MAB mapper 2 -> 5 equidistant bins in [0, 6]
        dist >6 -> 1 bin
        => 6 bins

NOTE on "exact west.cfg bins": the refining regions in west.cfg are
MABBinMappers (Minimal Adaptive Binning) whose edges are recomputed each WE
iteration from the walker min/max, so there is no single fixed grid to read off.
This script reconstructs a *static* grid that matches the cfg's structure and
per-region bin counts (Option A) -- fixed and reproducible for both runs.

The two w_assign macrostates (source = bound, sink = unbound) are kept as their
own dedicated anchor states; only the intermediate (w_assign "unknown", col 2 ==
2) frames are placed on the grid.

Output
------
shared_grid.pkl : dict with keys
    edges_rmsd, edges_ene, edges_dist : bin edges (interior + inf/finite ends)
    nbins        : (n_rmsd, n_ene, n_dist)
    bound_id, unbound_id, unknown_id  : reserved state ids
    dictionary   : {state_id (int) -> single unicode char}, last key = unknown
    char_base    : the chr() offset used
    cols         : which pathway columns hold each coord
shared_grid_legend.csv : human-readable state -> bin-range / centroid / counts
"""
import os
import pickle

import numpy as np

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = ["nomon", "wmon"]
PICKLE = "succ_traj/output.pickle"          # relative to each run dir

# pathway column layout (verified; see pdists/plot_pdist_traces.py):
#   0 iter | 1 seg | 2 w_assign state | 3 pcoord0 RMSD | 4 pcoord1 ene
#   5 pcoord2 ADPmin | 6 pcoord3 PO4min | 7 label | 8 frame | 9 weight
COL_STATE = 2
COL_RMSD = 3
COL_ENE = 4
COL_DIST = 5

# w_assign macrostate ids from the FULLY_UNBOUND2 scheme
WA_BOUND = 0
WA_UNBOUND = 1

# int-ene base split from west.cfg (['-inf', 10, 'inf']) and #bins per side
ENE_SPLIT = 10.0
ENE_NBINS_PER_SIDE = 5

CHAR_BASE = 0x100      # chr() offset -> unique single unicode chars, avoids '!'
OUT_GRID = os.path.join(HERE, "shared_grid.pkl")
OUT_LEGEND = os.path.join(HERE, "shared_grid_legend.csv")


# ---------------------------------------------------------------------------
# Fixed cfg-derived edges (RMSD, dist).  ene edges depend on data range.
# ---------------------------------------------------------------------------
def rmsd_edges():
    lo = list(np.linspace(0.0, 7.5, 6))            # 5 bins in [0, 7.5]
    mid = list(np.linspace(7.5, 10.5, 6))[1:]      # 5 bins in [7.5, 10.5]
    return np.array(lo + mid + [np.inf])           # + 1 bin for >10.5 -> 11 bins


def dist_edges():
    lo = list(np.linspace(0.0, 6.0, 6))            # 5 bins in [0, 6]
    return np.array(lo + [np.inf])                 # + 1 bin for >6 -> 6 bins


def ene_edges(ene_min, ene_max):
    below = list(np.linspace(ene_min, ENE_SPLIT, ENE_NBINS_PER_SIDE + 1))
    above = list(np.linspace(ENE_SPLIT, ene_max, ENE_NBINS_PER_SIDE + 1))[1:]
    return np.array(below + above)                 # 10 bins, edge at ENE_SPLIT


def bin_index(x, edges, nbins):
    """Vectorized bin index 0..nbins-1 for values x given full edges."""
    return np.clip(np.digitize(x, edges[1:-1]), 0, nbins - 1)


# ---------------------------------------------------------------------------
def load_frames(run):
    with open(os.path.join(HERE, run, PICKLE), "rb") as f:
        out = pickle.load(f)
    return np.concatenate([np.asarray(v, dtype=float) for v in out], axis=0)


def main():
    # 1. pooled frames -> shared int-ene range
    frames = {r: load_frames(r) for r in RUNS}
    pooled = np.concatenate([frames[r] for r in RUNS], axis=0)
    ene_min = float(pooled[:, COL_ENE].min())
    ene_max = float(pooled[:, COL_ENE].max())
    print(f"pooled frames: " + ", ".join(f"{r}={frames[r].shape[0]}" for r in RUNS))
    print(f"int-ene pooled range: [{ene_min:.2f}, {ene_max:.2f}] "
          f"(split at {ENE_SPLIT})")

    e_rmsd = rmsd_edges()
    e_ene = ene_edges(ene_min, ene_max)
    e_dist = dist_edges()
    n_rmsd = len(e_rmsd) - 1
    n_ene = len(e_ene) - 1
    n_dist = len(e_dist) - 1
    n_grid = n_rmsd * n_ene * n_dist
    print(f"grid: RMSD {n_rmsd} x ene {n_ene} x dist {n_dist} = {n_grid} cells")

    # 2. state-id layout: 0 = bound, 1 = unbound, 2.. = grid cells, last = unknown
    bound_id = 0
    unbound_id = 1
    grid_offset = 2
    unknown_id = grid_offset + n_grid

    dictionary = {i: chr(CHAR_BASE + i) for i in range(unknown_id)}
    dictionary[unknown_id] = "!"

    grid = {
        "edges_rmsd": e_rmsd, "edges_ene": e_ene, "edges_dist": e_dist,
        "nbins": (n_rmsd, n_ene, n_dist),
        "grid_offset": grid_offset,
        "bound_id": bound_id, "unbound_id": unbound_id, "unknown_id": unknown_id,
        "dictionary": dictionary, "char_base": CHAR_BASE,
        "cols": {"state": COL_STATE, "rmsd": COL_RMSD, "ene": COL_ENE,
                 "dist": COL_DIST},
        "wa_bound": WA_BOUND, "wa_unbound": WA_UNBOUND,
        "ene_min": ene_min, "ene_max": ene_max, "ene_split": ENE_SPLIT,
    }
    with open(OUT_GRID, "wb") as f:
        pickle.dump(grid, f)
    print(f"wrote {OUT_GRID}  ({len(dictionary)} states incl. unknown)")

    # 3. occupancy report + legend, using intermediate frames only (state 2)
    occ = {}
    for r in RUNS:
        fr = frames[r]
        inter = fr[~np.isin(fr[:, COL_STATE].astype(int), [WA_BOUND, WA_UNBOUND])]
        ri = bin_index(inter[:, COL_RMSD], e_rmsd, n_rmsd)
        ei = bin_index(inter[:, COL_ENE], e_ene, n_ene)
        di = bin_index(inter[:, COL_DIST], e_dist, n_dist)
        flat = ri * (n_ene * n_dist) + ei * n_dist + di
        vals, cnts = np.unique(flat, return_counts=True)
        occ[r] = dict(zip(vals.tolist(), cnts.tolist()))

    all_cells = sorted(set().union(*[set(occ[r]) for r in RUNS]))
    shared = set(occ[RUNS[0]])
    for r in RUNS[1:]:
        shared &= set(occ[r])
    print(f"occupied intermediate cells: total={len(all_cells)} "
          f"shared={len(shared)} "
          + " ".join(f"{r}-only={len(set(occ[r]) - shared)}" for r in RUNS))

    def cell_center(flat):
        ri, rem = divmod(flat, n_ene * n_dist)
        ei, di = divmod(rem, n_dist)

        def ctr(i, e):
            lo, hi = e[i], e[i + 1]
            if not np.isfinite(hi):
                return f">{lo:.2f}"
            return f"{lo:.2f}-{hi:.2f}"
        return ctr(ri, e_rmsd), ctr(ei, e_ene), ctr(di, e_dist)

    with open(OUT_LEGEND, "w") as f:
        f.write("state_id,char,kind,rmsd_bin,ene_bin,dist_bin,"
                + ",".join(f"n_{r}" for r in RUNS) + "\n")
        f.write(f"{bound_id},{dictionary[bound_id]!r},bound,,,,"
                + ",".join("" for _ in RUNS) + "\n")
        f.write(f"{unbound_id},{dictionary[unbound_id]!r},unbound,,,,"
                + ",".join("" for _ in RUNS) + "\n")
        for flat in all_cells:
            sid = grid_offset + flat
            rb, eb, db = cell_center(flat)
            counts = [str(occ[r].get(flat, 0)) for r in RUNS]
            f.write(f"{sid},{dictionary[sid]!r},grid,{rb},{eb},{db},"
                    + ",".join(counts) + "\n")
    print(f"wrote {OUT_LEGEND}")


if __name__ == "__main__":
    main()
