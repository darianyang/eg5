"""
Extract a representative continuous trajectory for each SHARED path-dendrogram
cluster (see cluster_paths_shared.py) and drop the coordinates in a separate dir.

The shared pathway clustering pools the successful pathways of both Eg5 runs
(nomon = no monastrol, wmon = + monastrol) into one LPATH distance matrix and
cuts one merged dendrogram.  Here, for each path cluster we pick a representative
pathway and rebuild its full atomistic trajectory from the WE segment files.

Representative selection (``--rep``):
    medoid  (default) : pathway with the smallest summed shared-space distance to
                        the other members of its cluster -- the most central path
                        by the very metric that defines the cluster.
    weight            : highest-WE-weight pathway (LPATH's native select_rep).

A pathway is a parent->child WE lineage: its rows are (iter, seg, ...) with one
segment per iteration (verified contiguous).  We concatenate each lineage
segment's ``seg.nc`` in ascending-iteration order into one continuous trajectory
(the combine_trace.py / concat_iter method), optionally superposing on protein
CA so global tumbling is removed for viewing.

The representative may come from EITHER condition; whichever is most central for
that shared cluster is used, and its origin is recorded.

Outputs, under ``<out_dir>/`` (default shared_paths/rep_paths/):
    cluster<k>_rep.nc     continuous AMBER NetCDF trajectory of the representative
    cluster<k>_rep.pdb     frame 0 (topology reference / quick view)
    cluster<k>_rep.csv     per-segment provenance (iter, seg, pathway weight)
    rep_paths_summary.csv  one row per cluster (origin, medoid index, weight, ...)

Usage:
    python extract_rep_paths.py [--k 3] [--rep medoid|weight] [--no-align]
        [--out-dir shared_paths/rep_paths] [--clusters 1 2 3]
"""
import os
import sys
import argparse
import pickle

import numpy as np
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
import mdtraj as md

# reuse the project's condition config (prmtop + traj_segs + west.h5 paths)
sys.path.insert(0, "/ix/lchong/dty7/eg5/allostery/scripts")
from eg5_allostery import CONDITIONS                       # noqa: E402

RUNS = ["nomon", "wmon"]
DISTMAT = "shared_paths/distmat_shared.npy"
REASSIGNED = "{run}/succ_traj/reassigned.pickle"


def load_pooled():
    """Load reassigned pathways for both runs in RUNS order.

    Returns paths (dict run->ndarray), offsets (cumulative counts) so a pooled
    global index maps back to (run, local index) exactly as the pooled distance
    matrix was built in cluster_paths_shared.py.
    """
    paths, counts = {}, []
    for run in RUNS:
        with open(REASSIGNED.format(run=run), "rb") as f:
            paths[run] = np.asarray(pickle.load(f))
        counts.append(len(paths[run]))
    offsets = np.cumsum([0] + counts)
    return paths, offsets


def locate(gidx, offsets):
    """Map pooled global index -> (run, local index)."""
    for r, run in enumerate(RUNS):
        if gidx < offsets[r + 1]:
            return run, int(gidx - offsets[r])
    raise IndexError(gidx)


def pathway_weight(pathway):
    """WE weight of a pathway: weight of its last non-unknown frame, matching
    lpath.match.gen_dist_matrix.  (dictionary unknown id is the max col-2 value;
    here col-2 in {0..3} real, 4 = unknown.)"""
    p = np.asarray(pathway)
    nz = p[p[:, 2] < p[:, 2].max()] if (p[:, 2] == 4).any() else p
    return float(nz[-1][-1])


def lineage_segments(pathway):
    """Ordered unique (iter, seg) lineage of a pathway, ascending in iteration.

    Rows are (iter, seg, ...); one seg per iteration for a WE lineage.  Zero-iter
    padding rows (if any) are dropped.
    """
    itseg = np.asarray(pathway)[:, :2].astype(int)
    seen, uniq = set(), []
    for it, sg in itseg:
        if it != 0 and (it, sg) not in seen:
            seen.add((it, sg))
            uniq.append((int(it), int(sg)))
    return sorted(uniq)


def build_trajectory(run, segs, align):
    """Concatenate each lineage segment's seg.nc into one continuous trajectory."""
    info = CONDITIONS[run]
    top = md.load_prmtop(info["prmtop"])
    chunks = []
    for it, sg in segs:
        seg_nc = f"{info['trajseg']}/{it:06d}/{sg:06d}/seg.nc"
        if not os.path.exists(seg_nc):
            print(f"    WARNING missing {seg_nc}", flush=True)
            continue
        chunks.append(md.load(seg_nc, top=top))
    traj = md.join(chunks) if len(chunks) > 1 else chunks[0]
    if align:
        ca = traj.top.select("protein and name CA")
        traj.superpose(traj, frame=0, atom_indices=ca)
    return traj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=3,
                    help="number of shared path clusters (must match the "
                         "dendrogram cut you chose; default 3)")
    ap.add_argument("--rep", choices=["medoid", "weight", "weight-per-cond"],
                    default="medoid",
                    help="representative selection: medoid (most central, "
                         "default); weight (single highest-WE-weight path, "
                         "pooled -- biased to the higher-weight condition); "
                         "weight-per-cond (highest-flux path within EACH "
                         "condition -> one nomon + one wmon rep per cluster)")
    ap.add_argument("--clusters", nargs="+", type=int, default=None,
                    help="which cluster ids to extract (default: all)")
    ap.add_argument("--out-dir", default="shared_paths/rep_paths")
    ap.add_argument("--no-align", action="store_true",
                    help="skip protein-CA superposition of output frames")
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    paths, offsets = load_pooled()

    dm = np.load(DISTMAT)
    z = sch.linkage(squareform(dm, checks=False), method="ward")
    labels = sch.fcluster(z, t=a.k, criterion="maxclust")     # 1..k
    all_cl = sorted(np.unique(labels))
    clusters = a.clusters if a.clusters is not None else all_cl

    def max_weight_member(members):
        """Pooled index of the highest-WE-weight pathway among ``members``."""
        w = [pathway_weight(paths[locate(g, offsets)[0]][locate(g, offsets)[1]])
             for g in members]
        return int(members[int(np.argmax(w))])

    # build the list of (cluster, rep_pooled_idx, tag) extraction jobs.  tag is
    # "" for a single rep per cluster, or the condition name in weight-per-cond
    # (which yields one nomon + one wmon rep per cluster).
    jobs = []
    for cl in clusters:
        members = np.where(labels == cl)[0]
        if len(members) == 0:
            print(f"cluster {cl}: empty, skipped", flush=True)
            continue
        if a.rep == "medoid":
            sub = dm[np.ix_(members, members)]
            jobs.append((cl, int(members[np.argmin(sub.sum(axis=1))]), ""))
        elif a.rep == "weight":
            jobs.append((cl, max_weight_member(members), ""))
        else:  # weight-per-cond: highest-flux path within each condition
            for run in RUNS:
                run_members = [g for g in members
                               if locate(g, offsets)[0] == run]
                if not run_members:
                    print(f"cluster {cl}: no {run} paths, skipped", flush=True)
                    continue
                jobs.append((cl, max_weight_member(run_members), run))

    summary = []
    for cl, rep, tag in jobs:
        run, li = locate(rep, offsets)
        pathway = paths[run][li]
        w = pathway_weight(pathway)
        segs = lineage_segments(pathway)
        suffix = f"_{tag}" if tag else ""
        print(f"cluster {cl}{suffix}: rep={a.rep} pooled_idx={rep} origin={run} "
              f"local={li} | {len(segs)} segments, iters "
              f"{segs[0][0]}..{segs[-1][0]}, weight={w:.3e}", flush=True)

        traj = build_trajectory(run, segs, align=not a.no_align)
        base = f"{a.out_dir}/cluster{cl}{suffix}_rep"
        traj.save_netcdf(f"{base}.nc")
        traj[0].save_pdb(f"{base}.pdb")

        # per-segment provenance
        with open(f"{base}.csv", "w") as f:
            f.write("order,iter,seg\n")
            for i, (it, sg) in enumerate(segs):
                f.write(f"{i},{it},{sg}\n")

        print(f"    -> {base}.nc ({traj.n_frames} frames, "
              f"{traj.n_atoms} atoms), {base}.pdb, {base}.csv", flush=True)
        summary.append((cl, run, rep, li, w, len(segs),
                        segs[0][0], segs[-1][0], traj.n_frames))

    with open(f"{a.out_dir}/rep_paths_summary.csv", "w") as f:
        f.write("cluster,origin,pooled_idx,local_idx,weight,n_segments,"
                "iter_start,iter_end,n_frames\n")
        for row in summary:
            f.write(",".join(str(x) for x in row) + "\n")
    print(f"\nDone. Representatives + provenance written to {a.out_dir}/")


if __name__ == "__main__":
    main()
