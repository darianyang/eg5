"""Rank the strongest-NMI residues within a cluster's allosteric network.

Reads the per-(condition, cluster, scheme) ``nmi.csv`` edge tables produced by
build_networks.py and reports, for the chosen cluster:

  * top residue *pairs* by NMI (includes trivial i,i+1 backbone neighbours),
  * top *non-local* pairs (sequence separation > ``--min-sep``) -- the
    long-range couplings that are the allosterically interesting ones, and
  * top *hub* residues by summed NMI strength (sum of a residue's edge weights).

By default it runs cluster 4 for both conditions (pooled scheme) and writes the
tables to ``<out-dir>/cluster<c>_top_nmi_<cond>.csv``.

    python scripts/nmi_top_residues.py [--cluster 4] [--conds nomon wmon]
        [--scheme pooled] [--topn 20] [--min-sep 5] [--out-dir plots/nmi_top]
"""
import os
import argparse
import numpy as np
import pandas as pd

NET = "networks"


def load_pairs(cond, cluster, scheme):
    """Undirected (u<v) residue-pair NMI table for one network, or None."""
    path = f"{NET}/{cond}/cluster{cluster}_{scheme}/nmi.csv"
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    uv = df["Residue Pair"].str.split("|", expand=True)
    u = uv[0].str.split().str[1].astype(int).to_numpy()
    v = uv[1].str.split().str[1].astype(int).to_numpy()
    val = df["MI Difference"].to_numpy()
    keep = u < v                                   # each pair is stored twice
    return pd.DataFrame({"u": u[keep], "v": v[keep], "nmi": val[keep]})


def hub_strength(pairs):
    """Summed NMI over each residue's edges, descending."""
    s = (pd.concat([pairs[["u", "nmi"]].rename(columns={"u": "res"}),
                    pairs[["v", "nmi"]].rename(columns={"v": "res"})])
         .groupby("res")["nmi"].sum().sort_values(ascending=False))
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cluster", type=int, default=4)
    ap.add_argument("--conds", nargs="+", default=["nomon", "wmon"])
    ap.add_argument("--scheme", default="pooled")
    ap.add_argument("--topn", type=int, default=20)
    ap.add_argument("--min-sep", type=int, default=5,
                    help="minimum |u-v| for a pair to count as non-local")
    ap.add_argument("--out-dir", default="plots/nmi_top")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    for cond in a.conds:
        p = load_pairs(cond, a.cluster, a.scheme)
        if p is None:
            print(f"{cond} cluster {a.cluster}: no nmi.csv, skipped")
            continue
        p["sep"] = (p["v"] - p["u"]).abs()

        print(f"\n===== {cond}  cluster {a.cluster}  ({a.scheme}) =====")
        print(f"{len(p)} edges, NMI {p.nmi.min():.4f}-{p.nmi.max():.4f}")

        top = p.sort_values("nmi", ascending=False).head(a.topn)
        nl = p[p.sep > a.min_sep].sort_values("nmi", ascending=False).head(a.topn)
        hubs = hub_strength(p).head(a.topn)

        print(f"\n-- top {a.topn} pairs (all) --")
        for _, r in top.iterrows():
            print(f"  Res {int(r.u):>3} - Res {int(r.v):>3}  "
                  f"(|d|={int(r.sep):>3})  NMI={r.nmi:.4f}")
        print(f"\n-- top {a.topn} non-local pairs (|d|>{a.min_sep}) --")
        for _, r in nl.iterrows():
            print(f"  Res {int(r.u):>3} - Res {int(r.v):>3}  "
                  f"(|d|={int(r.sep):>3})  NMI={r.nmi:.4f}")
        print(f"\n-- top {a.topn} hub residues (summed NMI) --")
        for res, val in hubs.items():
            print(f"  Res {int(res):>3}  strength={val:.3f}")

        # persist the tables
        top.to_csv(f"{a.out_dir}/cluster{a.cluster}_toppairs_{cond}.csv", index=False)
        nl.to_csv(f"{a.out_dir}/cluster{a.cluster}_nonlocal_{cond}.csv", index=False)
        hubs.rename("strength").to_csv(
            f"{a.out_dir}/cluster{a.cluster}_hubs_{cond}.csv")
    print(f"\nwrote tables to {a.out_dir}/")


if __name__ == "__main__":
    main()
