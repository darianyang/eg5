"""Extract representative frames from each LPATH cluster for visualization.

Each WE frame carries a shared 6-cluster label (0..5) in ``auxdata/labels`` --
the same labelling that defines the per-cluster allosteric networks.  For each
(condition, cluster) this pulls ``--nframes`` example conformations out of the
trajectory so the cluster (an ADP-exit stage) can be looked at in ChimeraX/VMD.

Frame -> label alignment mirrors eg5_allostery.py: ``seg.nc`` holds the 5
dynamics frames (pcoord points 1..len-1), so dynamics frame ``i`` of a segment
carries label ``labels[seg, i+1]``.  A frame belongs to cluster ``c`` iff its
own label is ``c`` (no endpoint/movement pairing here -- we want single
conformations, not transitions).

Selection samples *distinct segments* (one frame each), so no two frames come
from the same short seg.nc -- maximizing conformational diversity for viewing.
``--mode uniform`` (default) draws segments uniformly, giving an even spread
across the whole cluster; this matters for the late/rare ADP-exit clusters,
where WE weight concentrates on a handful of walkers so ``--mode weighted``
(probability ~ WE weight) collapses onto a narrow, near-duplicate set of frames.
Use ``weighted`` only when you want equilibrium-probability-representative
frames rather than a diverse visual survey.  ``--allow-repeat-seg`` samples
individual frames instead of one per segment.

Outputs, under ``<out_root>/<cond>/``:
  cluster<c>_rep.pdb  -- multi-model PDB, one MODEL per representative frame,
                         optionally superposed on protein CA (``--align``)
  cluster<c>_rep.csv  -- provenance: model, iter, seg, frame, weight, label

    python scripts/extract_rep_frames.py [--conds nomon wmon]
        [--clusters 0 1 2 3 4 5] [--nframes 10] [--mode weighted|uniform]
        [--seed 0] [--out-root rep_frames] [--no-align] [--allow-repeat-seg]
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
import mdtraj as md

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import CONDITIONS, N_CLUSTERS, _load_iter_meta  # noqa: E402


def cluster_frame_pool(westh5):
    """Enumerate every dynamics frame with its cluster label and WE weight.

    Returns a DataFrame with columns iter, seg, frame, weight, label where
    ``frame`` is the 0-based dynamics-frame index within that segment's seg.nc.
    """
    meta = _load_iter_meta(westh5)
    recs = []
    for k in sorted(meta):
        it = int(k.split("_")[-1])
        w, lab = meta[k]                 # w:(n_seg,)  lab:(n_seg, pcoord_len)
        nseg, plen = lab.shape
        nfr = plen - 1                   # dynamics frames in seg.nc
        # dynamics frame i  <->  label column i+1
        fr_lab = lab[:, 1:]              # (n_seg, nfr)
        seg_idx = np.repeat(np.arange(nseg), nfr)
        fr_idx = np.tile(np.arange(nfr), nseg)
        recs.append(pd.DataFrame({
            "iter": it,
            "seg": seg_idx,
            "frame": fr_idx,
            "weight": np.repeat(w, nfr),
            "label": fr_lab.reshape(-1),
        }))
    return pd.concat(recs, ignore_index=True)


def _weights(rows, mode):
    """Sampling probabilities for ``rows`` under ``mode`` (or None = uniform)."""
    if mode == "uniform":
        return None
    if mode != "weighted":
        raise ValueError(mode)
    p = rows["weight"].to_numpy(float)
    return p / p.sum()


def choose(pool_c, nframes, mode, allow_repeat_seg, rng):
    """Pick up to ``nframes`` rows from one cluster's frame pool."""
    if len(pool_c) == 0:
        return pool_c
    if allow_repeat_seg:
        cand = pool_c.reset_index(drop=True)
        n = min(nframes, len(cand))
        idx = rng.choice(len(cand), size=n, replace=False, p=_weights(cand, mode))
        return cand.iloc[idx].reset_index(drop=True)

    # one frame per distinct segment: sample segments (weighted by WE weight),
    # then a random cluster frame within each chosen segment
    segs = (pool_c.drop_duplicates(["iter", "seg"])[["iter", "seg", "weight"]]
            .reset_index(drop=True))
    n = min(nframes, len(segs))
    sidx = rng.choice(len(segs), size=n, replace=False, p=_weights(segs, mode))
    picks = []
    for _, s in segs.iloc[sidx].iterrows():
        g = pool_c[(pool_c["iter"] == s["iter"]) & (pool_c["seg"] == s["seg"])]
        picks.append(g.iloc[rng.integers(len(g))])
    return pd.DataFrame(picks).reset_index(drop=True)


def load_frames(chosen, cond):
    """Load the chosen frames (in row order) into one mdtraj Trajectory."""
    info = CONDITIONS[cond]
    top = md.load_prmtop(info["prmtop"])
    frames = []
    # load each seg.nc once, then pull the needed frame(s)
    for (it, seg), grp in chosen.groupby(["iter", "seg"], sort=False):
        seg_nc = f"{info['trajseg']}/{int(it):06d}/{int(seg):06d}/seg.nc"
        if not os.path.exists(seg_nc):
            print(f"    WARNING missing {seg_nc}", flush=True)
            continue
        tr = md.load(seg_nc, top=top)
        for _, r in grp.iterrows():
            frames.append((r["_order"], tr[int(r["frame"])]))
    frames.sort(key=lambda x: x[0])
    traj = frames[0][1]
    for _, fr in frames[1:]:
        traj = traj.join(fr)
    return traj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conds", nargs="+", default=["nomon", "wmon"])
    ap.add_argument("--clusters", nargs="+", type=int,
                    default=list(range(N_CLUSTERS)))
    ap.add_argument("--nframes", type=int, default=10)
    ap.add_argument("--mode", choices=["weighted", "uniform"], default="uniform")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-root", default="rep_frames")
    ap.add_argument("--no-align", action="store_true",
                    help="skip protein-CA superposition of the output frames")
    ap.add_argument("--allow-repeat-seg", action="store_true",
                    help="sample individual frames (may repeat a seg.nc) "
                         "instead of one frame per distinct segment")
    a = ap.parse_args()

    for cond in a.conds:
        info = CONDITIONS[cond]
        print(f"[{cond}] enumerating frames from {info['westh5']}", flush=True)
        pool = cluster_frame_pool(info["westh5"])
        rng = np.random.default_rng(a.seed)
        out_dir = f"{a.out_root}/{cond}"
        os.makedirs(out_dir, exist_ok=True)

        for c in a.clusters:
            pool_c = pool[pool["label"] == c]
            chosen = choose(pool_c, a.nframes, a.mode, a.allow_repeat_seg, rng)
            if len(chosen) == 0:
                print(f"  cluster {c}: no frames, skipped", flush=True)
                continue
            chosen = chosen.copy()
            chosen["_order"] = np.arange(len(chosen))

            traj = load_frames(chosen, cond)
            if not a.no_align:
                ca = traj.top.select("protein and name CA")
                traj.superpose(traj, frame=0, atom_indices=ca)

            pdb = f"{out_dir}/cluster{c}_rep.pdb"
            traj.save_pdb(pdb)
            man = chosen[["_order", "iter", "seg", "frame", "weight", "label"]]
            man = man.rename(columns={"_order": "model"})
            man["model"] += 1                      # PDB MODEL numbers are 1-based
            man.to_csv(f"{out_dir}/cluster{c}_rep.csv", index=False)
            print(f"  cluster {c}: {len(chosen)} frames "
                  f"(of {len(pool_c)} in cluster) -> {pdb}", flush=True)


if __name__ == "__main__":
    main()
