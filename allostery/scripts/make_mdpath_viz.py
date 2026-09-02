"""Generate mdpath 3D pathway visualizations (STL tube meshes + ChimeraX script)
for each (condition, cluster) allosteric network.

Reuses mdpath's own visualization pipeline on the paths we already saved
(networks/<cond>/cluster<c>_<scheme>/paths.pkl):

  top paths -> PatwayClustering (overlap-based path bundles) ->
  CA-coordinate backtracking -> quick JSON -> STL tube meshes (one per bundle).

For every network we then write a self-contained ChimeraX script (open.cxc) that
loads the shared reference structure and the bundle meshes in mdpath's colours,
so `chimerax networks/<cond>/cluster<c>_<scheme>/open.cxc` shows the pathways on
the protein.  Radius grows with how many shortest paths share a segment, so the
thick tubes are the dominant allosteric channels.

    python scripts/make_mdpath_viz.py --conds nomon wmon \
        --clusters 0 1 2 3 4 5 --scheme pooled --top 200 --ncpu 16
"""
import os
import sys
import json
import pickle
import argparse

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402

from mdpath.src.structure import StructureCalculations
from mdpath.src.cluster import PatwayClustering
from mdpath.src.visualization import MDPathVisualize, Colors


def rgb_str(c):
    return f"{c[0]:.4f},{c[1]:.4f},{c[2]:.4f}"


def write_cxc(out_dir, pdb_abspath, quick_json, cond, cluster, scheme):
    """ChimeraX script: protein cartoon + one STL surface per path bundle,
    each coloured with the bundle colour mdpath assigned it in the quick JSON."""
    quick = json.load(open(quick_json))
    # clusterid -> RGB colour (as emitted into the STL-driving quick JSON)
    id_color = {}
    for rec in quick:
        id_color.setdefault(rec["clusterid"], rec["color"])
    lines = [
        "# mdpath allosteric pathway bundles for "
        f"{cond} cluster {cluster} ({scheme})",
        "set bgColor white",
        f"open {pdb_abspath}",
        "hide atoms",
        "show cartoon",
        "color #1 gray(150)",
        "transparency #1 40 target c",
    ]
    for i, (cid, col) in enumerate(id_color.items()):
        model = i + 2  # #1 is the protein
        lines.append(f"open cluster_meshes/cluster_{cid}.stl")
        lines.append(f"color #{model} rgb({int(col[0]*255)},"
                     f"{int(col[1]*255)},{int(col[2]*255)})")
    lines += ["lighting soft", "graphics silhouettes true",
              "view", ""]
    with open(f"{out_dir}/open.cxc", "w") as f:
        f.write("\n".join(lines))


def process(cond, cluster, scheme, pdb, top, ncpu, closedist, out_root,
            radius_scale):
    net_dir = f"networks/{cond}/cluster{cluster}_{scheme}"
    ppath = f"{net_dir}/paths.pkl"
    if not os.path.exists(ppath):
        print(f"  skip {cond} c{cluster} {scheme}: no paths.pkl")
        return
    sorted_paths = pickle.load(open(ppath, "rb"))
    if not sorted_paths:
        print(f"  skip {cond} c{cluster} {scheme}: empty paths")
        return
    # mdpath's viz expects plain int node ids
    sorted_paths = [([int(x) for x in nodes], w) for nodes, w in sorted_paths]
    top_pathways = [p for p, _ in sorted_paths[:top]]

    struct = StructureCalculations(pdb)
    df_close = struct.calculate_residue_suroundings(closedist, "close")

    clustering = PatwayClustering(df_close, top_pathways, ncpu)
    clusters = clustering.pathways_cluster()  # {cluster_id: [pathway idx,...]}
    cluster_pathways_dict = clustering.pathway_clusters_dictionary(
        clusters, sorted_paths)

    coords = MDPathVisualize.residue_CA_coordinates(pdb, struct.last_res_num)
    updated = MDPathVisualize.apply_backtracking(cluster_pathways_dict, coords)
    formatted = MDPathVisualize.format_dict(updated)

    out_dir = f"{out_root}/{cond}/cluster{cluster}_{scheme}"
    os.makedirs(out_dir, exist_ok=True)
    quick = MDPathVisualize.precompute_cluster_properties_quick(formatted)
    # mdpath's base tube radius (0.015 A) is invisible on a ~50 A protein; scale
    # it up so the dominant (thick) channels read in ChimeraX. Relative
    # thickness -- how many shortest paths share a segment -- is preserved.
    for rec in quick:
        rec["radius"] *= radius_scale
    quick_json = f"{out_dir}/quick_precomputed_clusters_paths.json"
    json.dump(quick, open(quick_json, "w"), indent=2)

    # create_splines writes <dir>/cluster_meshes/cluster_<id>.stl next to the json
    MDPathVisualize.create_splines(quick_json)

    n_bundles = len(cluster_pathways_dict)
    pdb_abs = os.path.abspath(pdb)
    write_cxc(out_dir, pdb_abs, quick_json, cond, cluster, scheme)
    json.dump({"cond": cond, "cluster": cluster, "scheme": scheme,
               "n_paths_used": len(top_pathways), "n_bundles": n_bundles,
               "bundle_sizes": {int(k): len(v)
                                for k, v in cluster_pathways_dict.items()}},
              open(f"{out_dir}/viz_meta.json", "w"), indent=2)
    print(f"  {cond} c{cluster} {scheme}: {n_bundles} bundles from "
          f"{len(top_pathways)} paths -> {out_dir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conds", nargs="+", default=["nomon", "wmon"])
    ap.add_argument("--clusters", nargs="+", type=int,
                    default=list(range(N_CLUSTERS)))
    ap.add_argument("--scheme", default="pooled")
    ap.add_argument("--pdb", default="data/ref_nomon.pdb")
    ap.add_argument("--top", type=int, default=200,
                    help="top shortest paths fed to bundle clustering")
    ap.add_argument("--closedist", type=float, default=12.0)
    ap.add_argument("--ncpu", type=int, default=8)
    ap.add_argument("--radius-scale", type=float, default=10.0,
                    help="multiply mdpath tube radii for ChimeraX visibility")
    ap.add_argument("--out-root", default="viz")
    a = ap.parse_args()
    for cond in a.conds:
        for c in a.clusters:
            try:
                process(cond, c, a.scheme, a.pdb, a.top, a.ncpu,
                        a.closedist, a.out_root, a.radius_scale)
            except Exception as e:  # keep going across networks
                print(f"  ERROR {cond} c{c} {a.scheme}: {e}")


if __name__ == "__main__":
    main()
