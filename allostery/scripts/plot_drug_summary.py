"""One-figure summary of monastrol's effect on the Eg5 ADP-exit allosteric
network, resolved by LPATH cluster (exit stage).

Panels:
  A  network divergence per stage: total |Delta-NMI| (wmon - nomon), annotated
     with n_move so the low-sample (bias-prone) stages are obvious;
  B  dominant-pathway turnover: 1 - Jaccard of the top-path edge sets
     (how much the actual shortest-path channels are rewired by the drug);
  C  confidence per stage: pooled-vs-weighted edge correlation and Kish ESS;
  D  which residues: top drug-affected residues (max |Delta-NMI| involvement
     over stages), the shortlist to map onto structure.

    python scripts/plot_drug_summary.py [--scheme pooled]
"""
import os
import sys
import json
import pickle
import argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

NET = "networks"
CMP = "networks/comparisons"
PLOTS = "plots"


def path_edges(cond, c, scheme):
    p = f"{NET}/{cond}/cluster{c}_{scheme}/paths.pkl"
    if not os.path.exists(p):
        return None
    edges = set()
    for nodes, _ in pickle.load(open(p, "rb")):
        for a, b in zip(nodes[:-1], nodes[1:]):
            edges.add(tuple(sorted((int(a), int(b)))))
    return edges


def meta(cond, c, scheme):
    p = f"{NET}/{cond}/cluster{c}_{scheme}/meta.json"
    return json.load(open(p)) if os.path.exists(p) else {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scheme", default="pooled")
    args = ap.parse_args()
    s = args.scheme

    clusters = list(range(N_CLUSTERS))
    divergence, nmove_min, jturn = [], [], []
    ess_min = []
    for c in clusters:
        # A: total |delta| over all edges
        f = f"{CMP}/drug_cluster{c}_{s}.csv"
        divergence.append(pd.read_csv(f, usecols=["delta"])["delta"].abs().sum()
                          if os.path.exists(f) else np.nan)
        mn = meta("nomon", c, s)
        mw = meta("wmon", c, s)
        nmove_min.append(min(mn.get("n_move", 0), mw.get("n_move", 0)))
        ess_min.append(min(mn.get("ess", np.nan), mw.get("ess", np.nan)))
        # B: path-edge turnover
        ea, eb = path_edges("nomon", c, s), path_edges("wmon", c, s)
        if ea and eb:
            jturn.append(1 - len(ea & eb) / len(ea | eb))
        else:
            jturn.append(np.nan)

    # C: pooled vs weighted edge correlation
    agr = pd.read_csv(f"{CMP}/pooled_vs_weighted_agreement.csv")
    corr = {(r["cond"], int(r["cluster"])): r["edge_pearson_pooled_vs_weighted"]
            for _, r in agr.iterrows()}

    # D: top drug-affected residues across stages (max involvement)
    peak = {}
    for c in clusters:
        f = f"{CMP}/drug_cluster{c}_{s}_per_residue.csv"
        if not os.path.exists(f):
            continue
        for _, r in pd.read_csv(f).iterrows():
            res = int(r["residue"]); v = r["abs_delta_involvement"]
            if v > peak.get(res, (0, 0))[0]:
                peak[res] = (v, c)
    top = sorted(peak.items(), key=lambda kv: kv[1][0], reverse=True)[:15]

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    # A
    a = ax[0, 0]
    bars = a.bar(clusters, divergence, color="#b2182b")
    a.set_xlabel("LPATH cluster (ADP-exit stage)")
    a.set_ylabel(r"$\sum |\Delta$NMI$|$ over edges")
    a.set_title("A  Network divergence per stage (wmon vs nomon)")
    for c, b in zip(clusters, bars):
        a.annotate(f"n={nmove_min[c]:,}", (c, b.get_height()),
                   ha="center", va="bottom", fontsize=7, rotation=0)
    a.text(0.5, -0.22, "n = min movements across conditions; low-n stages have "
           "inflated NMI (finite-sample bias) --- read with panel C",
           transform=a.transAxes, ha="center", fontsize=7.5, color="dimgray")

    # B
    b = ax[0, 1]
    b.bar(clusters, jturn, color="#4575b4")
    b.set_xlabel("LPATH cluster (ADP-exit stage)")
    b.set_ylabel("1 - Jaccard(top-path edges)")
    b.set_title("B  Dominant-pathway turnover with drug")
    b.set_ylim(0, 1)

    # C
    cax = ax[1, 0]
    cn = [corr.get(("nomon", c), np.nan) for c in clusters]
    cw = [corr.get(("wmon", c), np.nan) for c in clusters]
    x = np.arange(len(clusters))
    cax.bar(x - 0.2, cn, 0.4, label="nomon", color="#2166ac")
    cax.bar(x + 0.2, cw, 0.4, label="wmon", color="#d6604d")
    cax.set_xticks(x); cax.set_xticklabels(clusters)
    cax.set_xlabel("LPATH cluster (ADP-exit stage)")
    cax.set_ylabel("pooled vs weighted edge r")
    cax.set_title("C  Weighting confidence (r; low ESS -> trust pooled)")
    cax.set_ylim(0, 1); cax.legend(fontsize=8, loc="upper left")
    ax2 = cax.twinx()
    ax2.plot(clusters, ess_min, "o-", color="black", lw=1, ms=4, label="min ESS")
    ax2.set_yscale("log"); ax2.set_ylabel("min Kish ESS (log)")
    ax2.legend(fontsize=8, loc="upper right")

    # D
    d = ax[1, 1]
    labels = [f"{res} (c{cst})" for res, (v, cst) in top][::-1]
    vals = [v for _, (v, _) in top][::-1]
    stage = [cst for _, (_, cst) in top][::-1]
    cmap = plt.get_cmap("tab10")
    d.barh(range(len(vals)), vals, color=[cmap(cs) for cs in stage])
    d.set_yticks(range(len(vals))); d.set_yticklabels(labels, fontsize=8)
    d.set_xlabel(r"peak $|\Delta$NMI$|$ involvement")
    d.set_title("D  Top drug-affected residues (label: res (peak stage))")

    fig.suptitle(f"Monastrol reshapes the Eg5 ADP-exit allosteric network "
                 f"({s}, phi backbone)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = f"{PLOTS}/drug_summary_{s}.png"
    fig.savefig(out, dpi=190)
    plt.close(fig)
    print("wrote", out)
    print("top residues:", [(res, round(v, 3), f"c{cst}")
                            for res, (v, cst) in top])


if __name__ == "__main__":
    main()
