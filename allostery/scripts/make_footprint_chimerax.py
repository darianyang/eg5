"""Paint monastrol's allosteric footprint onto the Eg5 structure for ChimeraX.

For each LPATH cluster we already have, per residue, the summed |Delta-NMI| over
its edges (networks/comparisons/drug_cluster<c>_<scheme>_per_residue.csv) --- how
much that residue's allosteric involvement changes between wmon and nomon at that
exit stage.  This writes, per cluster:

  viz/footprint/cluster<c>_<scheme>.defattr  -- ChimeraX residue attribute
  viz/footprint/cluster<c>_<scheme>.cxc      -- opens ref + colours by it

and a combined defattr with every cluster's attribute plus a peak-stage column,
so `color byattribute` shows where in the structure (and at which exit stage)
the drug reshapes communication most.

    python scripts/make_footprint_chimerax.py [--scheme pooled]
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402

CMP = "networks/comparisons"
OUT = "viz/footprint"
CHAIN = "A"


def load_involvement(cluster, scheme):
    p = f"{CMP}/drug_cluster{cluster}_{scheme}_per_residue.csv"
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    return dict(zip(df["residue"].astype(int), df["abs_delta_involvement"]))


def write_defattr(path, attr_name, res_vals):
    lines = [f"attribute: {attr_name}", "match mode: any",
             "recipient: residues"]
    for r in sorted(res_vals):
        lines.append(f"\t/{CHAIN}:{int(r)}\t{res_vals[r]:.6g}")
    open(path, "w").write("\n".join(lines) + "\n")


def write_cxc(path, pdb_abs, defattr_rel, attr_name, vmax, title):
    lines = [
        f"# {title}",
        "set bgColor white",
        f"open {pdb_abs}",
        "hide atoms",
        "show cartoon",
        f"open {defattr_rel}",
        # grey->hot palette; 0 = no change, vmax = strongest drug effect
        f"color byattribute {attr_name} palette "
        f"0,#dcdcdc:{vmax*0.5:.6g},#f4a000:{vmax:.6g},#b2182b "
        "target c",
        "lighting soft",
        "graphics silhouettes true",
        "view",
        "",
    ]
    open(path, "w").write("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scheme", default="pooled")
    ap.add_argument("--pdb", default="data/ref_nomon.pdb")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    pdb_abs = os.path.abspath(args.pdb)

    all_res = set()
    per_cluster = {}
    for c in range(N_CLUSTERS):
        d = load_involvement(c, args.scheme)
        if d is None:
            continue
        per_cluster[c] = d
        all_res |= set(d)

    if not per_cluster:
        print("no per-residue comparison CSVs found; run compare_networks.py first")
        return

    # shared colour ceiling across clusters (95th pct of all involvement values)
    allvals = np.concatenate([list(d.values()) for d in per_cluster.values()])
    vmax = float(np.percentile(allvals, 95))

    for c, d in per_cluster.items():
        attr = f"dnmi_c{c}"
        da = f"{OUT}/cluster{c}_{args.scheme}.defattr"
        write_defattr(da, attr, d)
        write_cxc(f"{OUT}/cluster{c}_{args.scheme}.cxc", pdb_abs,
                  os.path.basename(da), attr, vmax,
                  f"Monastrol allosteric footprint, cluster {c} ({args.scheme})")

    # combined: one defattr with every cluster attribute + peak stage + peak value
    combined = f"{OUT}/all_clusters_{args.scheme}.defattr"
    blocks = []
    for c, d in per_cluster.items():
        b = [f"attribute: dnmi_c{c}", "match mode: any", "recipient: residues"]
        for r in sorted(all_res):
            b.append(f"\t/{CHAIN}:{int(r)}\t{d.get(r, 0.0):.6g}")
        blocks.append("\n".join(b))
    # peak involvement across stages, and which stage
    peak_val, peak_stage = {}, {}
    for r in all_res:
        vals = {c: per_cluster[c].get(r, 0.0) for c in per_cluster}
        cbest = max(vals, key=vals.get)
        peak_val[r] = vals[cbest]
        peak_stage[r] = cbest
    for name, dd in (("dnmi_peak", peak_val), ("dnmi_peakstage", peak_stage)):
        b = [f"attribute: {name}", "match mode: any", "recipient: residues"]
        for r in sorted(all_res):
            b.append(f"\t/{CHAIN}:{int(r)}\t{dd[r]:.6g}")
        blocks.append("\n".join(b))
    open(combined, "w").write("\n\n".join(blocks) + "\n")

    write_cxc(f"{OUT}/all_clusters_peak_{args.scheme}.cxc", pdb_abs,
              os.path.basename(combined), "dnmi_peak", vmax,
              f"Monastrol footprint, peak over exit stages ({args.scheme})")

    print(f"wrote {len(per_cluster)} per-cluster + combined footprint files to {OUT}")
    print(f"colour ceiling vmax (95th pct) = {vmax:.4g}")
    print("open e.g.:  chimerax viz/footprint/all_clusters_peak_"
          f"{args.scheme}.cxc")


if __name__ == "__main__":
    main()
