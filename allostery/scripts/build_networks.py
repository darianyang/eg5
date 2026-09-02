"""Build per-(condition, cluster) allosteric networks with mdpath.

For each condition x cluster x scheme (pooled / WE-weighted) this:
  1. slices the phi movements of that cluster into an mdpath-style DataFrame,
  2. computes Normalized Mutual Information between residue pairs (pooled uses
     mdpath's NMICalculator; weighted uses a WE-weight-aware variant),
  3. builds mdpath's residue graph on the shared reference PDB and assigns the
     NMI edge weights,
  4. extracts the top max-weight shortest paths between distant residues.

Outputs, under networks/<cond>/cluster<c>_<scheme>/:
  nmi.csv        -- residue-pair NMI edge weights (the comparison currency)
  paths.pkl      -- list of (path, total_weight), sorted, top ``numpath``
  meta.json      -- n_move, ESS, params

Only phi drives the mdpath network (as in mdpath itself); chi1/chi2 movements
are carried in the movements file for later sidechain-aware extensions.
"""
import os
import sys
import json
import pickle
import argparse
import numpy as np
import pandas as pd
import h5py

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402

from mdpath.src.mutual_information import NMICalculator
from mdpath.src.structure import StructureCalculations
from mdpath.src.graph import GraphBuilder


def weighted_nmi(df, weights, num_bins=35, hist_range=(-180.0, 180.0)):
    """WE-weighted NMI edge table, matching NMICalculator's output format.

    Identical definition to mdpath's histogram NMI, but every movement (row)
    contributes its normalized WE weight to the 1-D and 2-D histograms instead
    of a unit count.  MI and entropy are scale-invariant in the histogram, so a
    weighted contingency is a valid weighted MI estimate.

    Returns a DataFrame with columns ``Residue Pair`` (tuple of "Res u","Res v")
    and ``MI Difference`` -- the shape GraphBuilder expects.
    """
    from itertools import combinations
    from scipy.stats import entropy
    from sklearn.metrics import mutual_info_score

    w = np.asarray(weights, dtype=np.float64)
    w = w / w.sum()
    cols = df.columns.tolist()
    rng = [hist_range, hist_range]

    hist1d, ent = {}, {}
    for c in cols:
        h, _ = np.histogram(df[c].to_numpy(), bins=num_bins,
                            range=hist_range, weights=w)
        hist1d[c] = h
        ent[c] = entropy(h)

    rows = {}
    for c1, c2 in combinations(cols, 2):
        hj, _, _ = np.histogram2d(df[c1].to_numpy(), df[c2].to_numpy(),
                                  bins=num_bins, range=rng, weights=w)
        mi = mutual_info_score(hist1d[c1], hist1d[c2], contingency=hj)
        denom = np.sqrt(ent[c1] * ent[c2])
        nmi = mi / denom if denom > 0 else 0.0
        rows[(c1, c2)] = nmi
        rows[(c2, c1)] = nmi
    return pd.DataFrame(rows.items(), columns=["Residue Pair", "MI Difference"])


def load_cluster_df(mov_h5, cluster, angle="phi"):
    """Return (DataFrame of movements, weights) for one cluster's phi movements."""
    with h5py.File(mov_h5, "r") as f:
        cl = f["cluster"][:]
        sel = cl == cluster
        mov = f[f"movements_{angle}"][sel]
        rid = f[f"res_ids_{angle}"][:]
        w = f["weight"][sel]
    df = pd.DataFrame(mov, columns=[f"Res {int(r)}" for r in rid])
    return df, w


def kish_ess(w):
    w = np.asarray(w, float)
    w = w / w.sum()
    return float(1.0 / np.sum(w ** 2))


def build_one(cond, cluster, scheme, mov_h5, pdb, out_root,
              num_bins=35, graphdist=5.0, fardist=12.0, numpath=500,
              min_moves=2000, ncpu=1):
    """Build and save one network.  Returns the output directory or None."""
    df, w = load_cluster_df(mov_h5, cluster)
    n = len(df)
    out_dir = f"{out_root}/{cond}/cluster{cluster}_{scheme}"
    os.makedirs(out_dir, exist_ok=True)
    ess = kish_ess(w) if n else 0.0

    if n < min_moves:
        json.dump({"cond": cond, "cluster": cluster, "scheme": scheme,
                   "n_move": n, "ess": ess, "skipped": "too_few_moves"},
                  open(f"{out_dir}/meta.json", "w"), indent=2)
        print(f"  {cond} c{cluster} {scheme}: n={n} < {min_moves}, skipped")
        return None

    if scheme == "pooled":
        nmi_df = NMICalculator(df, num_bins=num_bins, show_progress=False).nmi_df
    elif scheme == "weighted":
        nmi_df = weighted_nmi(df, w, num_bins=num_bins)
    else:
        raise ValueError(scheme)

    struct = StructureCalculations(pdb)
    gb = GraphBuilder(pdb, struct.last_res_num, nmi_df, graphdist)
    graph = gb.graph
    df_far = struct.calculate_residue_suroundings(fardist, "far")
    if ncpu > 1:
        paths = gb.collect_path_total_weights_parallel(df_far, ncpu)
    else:
        paths = gb.collect_path_total_weights(df_far)
    paths.sort(key=lambda x: x[1], reverse=True)
    paths = paths[:numpath]

    nmi_out = nmi_df.copy()
    nmi_out["Residue Pair"] = nmi_out["Residue Pair"].apply(lambda p: f"{p[0]}|{p[1]}")
    nmi_out.to_csv(f"{out_dir}/nmi.csv", index=False)
    pickle.dump(paths, open(f"{out_dir}/paths.pkl", "wb"))
    json.dump({"cond": cond, "cluster": cluster, "scheme": scheme,
               "n_move": n, "ess": ess, "n_graph_edges": graph.number_of_edges(),
               "n_paths": len(paths), "num_bins": num_bins,
               "graphdist": graphdist, "fardist": fardist},
              open(f"{out_dir}/meta.json", "w"), indent=2)
    print(f"  {cond} c{cluster} {scheme}: n={n} ess={ess:.0f} "
          f"edges={graph.number_of_edges()} paths={len(paths)} -> {out_dir}")
    return out_dir


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--conds", nargs="+", default=["nomon", "wmon"])
    p.add_argument("--schemes", nargs="+", default=["pooled", "weighted"])
    p.add_argument("--clusters", nargs="+", type=int,
                   default=list(range(N_CLUSTERS)))
    p.add_argument("--data-dir", default="data")
    p.add_argument("--pdb", default="data/ref_nomon.pdb",
                   help="shared reference PDB for the graph")
    p.add_argument("--out-root", default="networks")
    p.add_argument("--num-bins", type=int, default=35)
    p.add_argument("--numpath", type=int, default=500)
    p.add_argument("--min-moves", type=int, default=2000)
    p.add_argument("--ncpu", type=int, default=1)
    a = p.parse_args()

    for cond in a.conds:
        mov_h5 = f"{a.data_dir}/movements_{cond}.h5"
        for cluster in a.clusters:
            for scheme in a.schemes:
                build_one(cond, cluster, scheme, mov_h5, a.pdb, a.out_root,
                          num_bins=a.num_bins, numpath=a.numpath,
                          min_moves=a.min_moves, ncpu=a.ncpu)


if __name__ == "__main__":
    main()
