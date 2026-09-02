"""Compare the per-(condition, cluster) allosteric networks built by
build_networks.py, along the two axes the study cares about:

  * drug effect  -- wmon vs nomon within each cluster (monastrol's allosteric
    footprint at each stage of ADP exit);
  * path/stage   -- adjacent clusters within a condition (how the network
    rewires as ADP leaves).

The comparison currency is the per-edge NMI on the shared reference graph.  For
each comparison we write a ranked edge Delta-NMI table and a per-residue
"allosteric involvement change" (sum of |Delta-NMI| over a residue's edges),
which maps directly onto the structure for visualization.  Top-path edge sets
are compared by Jaccard overlap.  A summary table collects n_move and Kish ESS
per network so low-confidence (low-ESS) networks are easy to flag.

Outputs go to plots/ and networks/comparisons/.
"""
import os
import sys
import json
import glob
import itertools
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

NET_ROOT = "networks"
OUT = "networks/comparisons"
PLOTS = "plots"


def _net_dir(cond, cluster, scheme):
    return f"{NET_ROOT}/{cond}/cluster{cluster}_{scheme}"


def load_nmi(cond, cluster, scheme):
    """Return {(u,v): nmi} for one network, or None if it was skipped/missing."""
    path = f"{_net_dir(cond, cluster, scheme)}/nmi.csv"
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    out = {}
    for pair, val in zip(df["Residue Pair"], df["MI Difference"]):
        u, v = pair.split("|")
        key = tuple(sorted((int(u.split()[1]), int(v.split()[1]))))
        out[key] = val
    return out


def load_meta(cond, cluster, scheme):
    path = f"{_net_dir(cond, cluster, scheme)}/meta.json"
    return json.load(open(path)) if os.path.exists(path) else None


def load_path_edges(cond, cluster, scheme):
    """Set of undirected edges appearing in a network's top paths."""
    import pickle
    path = f"{_net_dir(cond, cluster, scheme)}/paths.pkl"
    if not os.path.exists(path):
        return None
    paths = pickle.load(open(path, "rb"))
    edges = set()
    for nodes, _ in paths:
        for a, b in zip(nodes[:-1], nodes[1:]):
            edges.add(tuple(sorted((int(a), int(b)))))
    return edges


def delta_edges(nmi_a, nmi_b, label_a, label_b, out_csv):
    """Ranked per-edge Delta = nmi_b - nmi_a and per-residue |Delta| involvement."""
    keys = sorted(set(nmi_a) | set(nmi_b))
    rows = [(u, v, nmi_a.get((u, v), 0.0), nmi_b.get((u, v), 0.0),
             nmi_b.get((u, v), 0.0) - nmi_a.get((u, v), 0.0)) for u, v in keys]
    edf = pd.DataFrame(rows, columns=["res1", "res2",
                                      f"nmi_{label_a}", f"nmi_{label_b}", "delta"])
    edf = edf.reindex(edf["delta"].abs().sort_values(ascending=False).index)
    edf.to_csv(out_csv, index=False)

    res = {}
    for _, r in edf.iterrows():
        res[r.res1] = res.get(r.res1, 0.0) + abs(r.delta)
        res[r.res2] = res.get(r.res2, 0.0) + abs(r.delta)
    rdf = (pd.DataFrame(res.items(), columns=["residue", "abs_delta_involvement"])
           .sort_values("residue"))
    rdf.to_csv(out_csv.replace(".csv", "_per_residue.csv"), index=False)
    return edf, rdf


def jaccard(a, b):
    if a is None or b is None:
        return np.nan
    return len(a & b) / len(a | b) if (a | b) else np.nan


def main(scheme_focus="pooled"):
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(PLOTS, exist_ok=True)
    conds = ["nomon", "wmon"]
    schemes = ["pooled", "weighted"]

    # summary table: n_move, ESS, built? per network
    summ = []
    for cond in conds:
        for c in range(N_CLUSTERS):
            for s in schemes:
                m = load_meta(cond, c, s)
                if m is None:
                    continue
                summ.append({"cond": cond, "cluster": c, "scheme": s,
                             "n_move": m.get("n_move"), "ess": m.get("ess"),
                             "built": "skipped" not in m,
                             "n_paths": m.get("n_paths")})
    summ = pd.DataFrame(summ)
    summ.to_csv(f"{OUT}/network_summary.csv", index=False)
    print("Network summary:\n", summ.to_string(index=False))

    # --- drug effect: wmon vs nomon within each cluster ---
    per_res_by_cluster = {}
    for c in range(N_CLUSTERS):
        a = load_nmi("nomon", c, scheme_focus)
        b = load_nmi("wmon", c, scheme_focus)
        if a is None or b is None:
            continue
        edf, rdf = delta_edges(a, b, "nomon", "wmon",
                               f"{OUT}/drug_cluster{c}_{scheme_focus}.csv")
        per_res_by_cluster[c] = rdf.set_index("residue")["abs_delta_involvement"]
        je = jaccard(load_path_edges("nomon", c, scheme_focus),
                     load_path_edges("wmon", c, scheme_focus))
        print(f"cluster {c}: drug Delta computed; top-path edge Jaccard "
              f"nomon/wmon = {je:.3f}")

    # per-residue drug-effect heatmap (residue x cluster)
    if per_res_by_cluster:
        allres = sorted(set().union(*[s.index for s in per_res_by_cluster.values()]))
        M = np.full((len(allres), N_CLUSTERS), np.nan)
        ridx = {r: i for i, r in enumerate(allres)}
        for c, s in per_res_by_cluster.items():
            for r, v in s.items():
                M[ridx[r], c] = v
        fig, ax = plt.subplots(figsize=(6, 10))
        im = ax.imshow(M, aspect="auto", cmap="magma",
                       extent=[-0.5, N_CLUSTERS - 0.5, allres[-1], allres[0]])
        ax.set_xlabel("cluster (ADP-exit stage)")
        ax.set_ylabel("residue")
        ax.set_title(f"Monastrol allosteric footprint\n"
                     f"sum |Delta-NMI| per residue ({scheme_focus})")
        fig.colorbar(im, ax=ax, label="|Delta-NMI| involvement")
        fig.tight_layout()
        fig.savefig(f"{PLOTS}/drug_footprint_{scheme_focus}.png", dpi=200)
        plt.close(fig)
        print(f"wrote {PLOTS}/drug_footprint_{scheme_focus}.png")

    # --- stage effect: adjacent clusters within each condition ---
    for cond in conds:
        for c1, c2 in itertools.pairwise(range(N_CLUSTERS)):
            a = load_nmi(cond, c1, scheme_focus)
            b = load_nmi(cond, c2, scheme_focus)
            if a is None or b is None:
                continue
            delta_edges(a, b, f"c{c1}", f"c{c2}",
                        f"{OUT}/stage_{cond}_c{c1}_c{c2}_{scheme_focus}.csv")

    # --- pooled vs weighted agreement per network (edge correlation) ---
    rows = []
    for cond in conds:
        for c in range(N_CLUSTERS):
            p = load_nmi(cond, c, "pooled")
            w = load_nmi(cond, c, "weighted")
            if p is None or w is None:
                continue
            keys = sorted(set(p) & set(w))
            pv = np.array([p[k] for k in keys])
            wv = np.array([w[k] for k in keys])
            r = np.corrcoef(pv, wv)[0, 1] if len(keys) > 2 else np.nan
            rows.append({"cond": cond, "cluster": c, "edge_pearson_pooled_vs_weighted": r})
    if rows:
        agr = pd.DataFrame(rows)
        agr.to_csv(f"{OUT}/pooled_vs_weighted_agreement.csv", index=False)
        print("\nPooled vs weighted edge agreement:\n", agr.to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "pooled")
