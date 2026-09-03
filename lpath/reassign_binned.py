"""
LPATH match-step reassignment onto the SHARED fixed bin grid.

Replaces reassign_custom.py (6 agglomerative shared clusters) with a finer,
fixed, west.cfg-derived rectilinear grid on (ADP-Mg RMSD, interaction energy,
ADP min contact distance).  The grid + the state->char dictionary are built
ONCE by build_shared_grid.py and saved to shared_grid.pkl, so nomon and wmon
share an identical alphabet and their pathway dendrograms stay comparable
(cluster_paths_shared_binned.py asserts the dictionaries match).

State assignment per frame:
    w_assign state (col 2) == 0 (bound)   -> dedicated 'bound'  anchor state
    w_assign state (col 2) == 1 (unbound) -> dedicated 'unbound' anchor state
    otherwise (intermediate)              -> its 3D grid-cell state

`lpath match -ra reassign_binned.reassign_binned ...`  (module must be importable
from cwd; the run script copies it into each run dir, and loads the grid from
../shared_grid.pkl).
"""
import os
import pickle

import numpy


def _load_grid():
    """Locate and load shared_grid.pkl (built by build_shared_grid.py)."""
    candidates = [
        os.environ.get("LPATH_SHARED_GRID"),
        os.path.join(os.getcwd(), "..", "shared_grid.pkl"),  # cwd = <run>/
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "shared_grid.pkl"),
        "shared_grid.pkl",
    ]
    for path in candidates:
        if path and os.path.exists(path):
            with open(path, "rb") as f:
                return pickle.load(f)
    raise FileNotFoundError(
        "shared_grid.pkl not found. Run build_shared_grid.py first "
        f"(looked in: {[c for c in candidates if c]}).")


def _bin_index(x, edges, nbins):
    return numpy.clip(numpy.digitize(x, edges[1:-1]), 0, nbins - 1)


def reassign_binned(data, pathways, dictionary, assign_file=None):
    """Reassign frame states to the shared bin grid.

    Parameters mirror lpath.match.reassign_custom: ``data`` is the list of
    successful pathways (variable length), ``pathways`` is the padded output
    array to fill, ``dictionary`` is the (empty) state->char map to return.
    """
    grid = _load_grid()
    c = grid["cols"]
    e_rmsd, e_ene, e_dist = grid["edges_rmsd"], grid["edges_ene"], grid["edges_dist"]
    n_rmsd, n_ene, n_dist = grid["nbins"]
    off = grid["grid_offset"]
    bound_id, unbound_id = grid["bound_id"], grid["unbound_id"]
    wa_bound, wa_unbound = grid["wa_bound"], grid["wa_unbound"]

    for idx, pathway in enumerate(data):
        arr = numpy.asarray(pathway, dtype=float)
        n = arr.shape[0]

        ri = _bin_index(arr[:, c["rmsd"]], e_rmsd, n_rmsd)
        ei = _bin_index(arr[:, c["ene"]], e_ene, n_ene)
        di = _bin_index(arr[:, c["dist"]], e_dist, n_dist)
        sid = off + ri * (n_ene * n_dist) + ei * n_dist + di

        wa = arr[:, c["state"]].astype(int)
        sid[wa == wa_bound] = bound_id
        sid[wa == wa_unbound] = unbound_id

        arr[:, c["state"]] = sid
        pathways[idx, :n] = arr

    # dictionary was built (with the unknown sentinel last) in build_shared_grid
    return dict(grid["dictionary"])
