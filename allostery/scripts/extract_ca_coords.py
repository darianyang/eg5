"""Extract per-cluster Ca coordinate ensembles for DF / SPM allosteric analysis.

Distance Fluctuation (DF) and Shortest Path Map (SPM) -- the two methods from
Rodriguez-Santos et al. (JCIM 2026) -- operate on Cartesian Ca coordinates of a
conformational ensemble, not on the dihedral *movements* used by the mdpath NMI
networks.  This pulls a large uniform sample of each LPATH cluster's frames and
stores their Ca coordinates so DF and SPM can be built per (condition, cluster)
and compared between nomon (WT) and wmon (+monastrol).

Frame -> label alignment mirrors eg5_allostery.py: ``seg.nc`` holds the 5
dynamics frames (pcoord points 1..len-1), so dynamics frame ``i`` of a segment
carries label ``labels[seg, i+1]``.  A frame joins cluster ``c`` iff its own
label is ``c``.

Sampling is *uniform* over each cluster's frames (WE weights are deliberately
NOT used here -- they are not converged -- but the per-frame weight is stored so
a weighted variant can be built later).  Coordinates are stored UNALIGNED; DF is
alignment-free, and SPM does its own Ca superposition against a per-cluster
reference at analysis time.

Writes one HDF5 per condition with:
    ca       (n_frame, n_res, 3) float32 -- Ca coordinates (nm, mdtraj units)
    res_ids  (n_res,)            int32   -- AMBER residue id per Ca column
    cluster  (n_frame,)          int8    -- LPATH cluster id of each frame
    weight   (n_frame,)          float64 -- WE walker weight (stored, unused here)

    python scripts/extract_ca_coords.py COND OUT [--max-frames 20000]
        [--seed 0] [--ncpu 16] [--clusters 0 1 2 3 4 5]
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
import h5py
import mdtraj as md

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import CONDITIONS, N_CLUSTERS, _load_iter_meta  # noqa: E402


def cluster_frame_pool(westh5):
    """Enumerate every dynamics frame: iter, seg, frame, weight, label."""
    meta = _load_iter_meta(westh5)
    recs = []
    for k in sorted(meta):
        it = int(k.split("_")[-1])
        w, lab = meta[k]                        # w:(n_seg,)  lab:(n_seg, pcoord_len)
        nseg, plen = lab.shape
        nfr = plen - 1                          # dynamics frames in seg.nc
        fr_lab = lab[:, 1:]                     # frame i <-> label column i+1
        recs.append(pd.DataFrame({
            "iter": it,
            "seg": np.repeat(np.arange(nseg), nfr),
            "frame": np.tile(np.arange(nfr), nseg),
            "weight": np.repeat(w, nfr),
            "label": fr_lab.reshape(-1),
        }))
    return pd.concat(recs, ignore_index=True)


def sample_pool(pool, clusters, max_frames, rng):
    """Uniformly sample up to ``max_frames`` rows per cluster."""
    picks = []
    for c in clusters:
        pc = pool[pool["label"] == c]
        n = min(max_frames, len(pc))
        if n == 0:
            print(f"  cluster {c}: no frames", flush=True)
            continue
        idx = rng.choice(len(pc), size=n, replace=False)
        picks.append(pc.iloc[idx])
        print(f"  cluster {c}: sampled {n} of {len(pc)} frames", flush=True)
    return pd.concat(picks, ignore_index=True)


_WORKER = {}


def _init_worker(prmtop):
    top = md.load_prmtop(prmtop)
    ca = top.select("protein and name CA")
    _WORKER["top"] = top
    _WORKER["ca"] = ca
    _WORKER["res_ids"] = np.array([top.atom(a).residue.resSeq + 1 for a in ca],
                                  dtype=np.int32)


def _worker_segment(task):
    """(coords[k,n_res,3], labels[k], weights[k], order[k]) for one segment."""
    seg_nc, weight, rows = task                 # rows: list of (frame, label, order)
    if not os.path.exists(seg_nc):
        return None
    try:
        tr = md.load(seg_nc, top=_WORKER["top"], atom_indices=_WORKER["ca"])
    except Exception:
        return None
    fr = np.array([r[0] for r in rows], dtype=int)
    if fr.max() >= tr.n_frames:                 # guard against short segments
        keep = fr < tr.n_frames
        rows = [r for r, k in zip(rows, keep) if k]
        fr = fr[keep]
        if fr.size == 0:
            return None
    coords = tr.xyz[fr].astype(np.float32)      # (k, n_res, 3)
    labels = np.array([r[1] for r in rows], dtype=np.int8)
    order = np.array([r[2] for r in rows], dtype=np.int64)
    weights = np.full(len(rows), weight, dtype=np.float64)
    return coords, labels, weights, order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cond", choices=list(CONDITIONS))
    ap.add_argument("out")
    ap.add_argument("--max-frames", type=int, default=20000,
                    help="max frames sampled per cluster (uniform)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ncpu", type=int, default=1)
    ap.add_argument("--clusters", nargs="+", type=int,
                    default=list(range(N_CLUSTERS)))
    a = ap.parse_args()

    info = CONDITIONS[a.cond]
    print(f"[{a.cond}] enumerating frames from {info['westh5']}", flush=True)
    pool = cluster_frame_pool(info["westh5"])
    rng = np.random.default_rng(a.seed)
    chosen = sample_pool(pool, a.clusters, a.max_frames, rng)
    chosen = chosen.reset_index(drop=True)
    chosen["order"] = np.arange(len(chosen))
    print(f"[{a.cond}] {len(chosen)} frames sampled total; loading Ca coords "
          f"on {a.ncpu} cpu(s)", flush=True)

    # group chosen frames by (iter, seg) so each seg.nc is read once
    tasks = []
    for (it, seg), g in chosen.groupby(["iter", "seg"], sort=False):
        seg_nc = f"{info['trajseg']}/{int(it):06d}/{int(seg):06d}/seg.nc"
        rows = list(zip(g["frame"].astype(int), g["label"].astype(int),
                        g["order"].astype(int)))
        tasks.append((seg_nc, float(g["weight"].iloc[0]), rows))
    print(f"[{a.cond}] {len(tasks)} segments to read", flush=True)

    coords_all, lab_all, wt_all, ord_all = [], [], [], []
    n_missing = 0

    def _collect(i, r):
        nonlocal n_missing
        if r is None:
            n_missing += 1
            return
        c, l, w, o = r
        coords_all.append(c); lab_all.append(l); wt_all.append(w); ord_all.append(o)
        if i % 2000 == 0:
            print(f"  {i}/{len(tasks)} segments read", flush=True)

    if a.ncpu > 1:
        from multiprocessing import Pool
        with Pool(a.ncpu, initializer=_init_worker,
                  initargs=(info["prmtop"],)) as pool_:
            for i, r in enumerate(pool_.imap(_worker_segment, tasks, chunksize=16)):
                _collect(i, r)
        _init_worker(info["prmtop"])            # for res_ids in the parent
    else:
        _init_worker(info["prmtop"])
        for i, t in enumerate(tasks):
            _collect(i, _worker_segment(t))

    coords = np.concatenate(coords_all)
    labels = np.concatenate(lab_all)
    weights = np.concatenate(wt_all)
    order = np.concatenate(ord_all)
    srt = np.argsort(order)                     # deterministic, sample order
    coords, labels, weights = coords[srt], labels[srt], weights[srt]

    with h5py.File(a.out, "w") as f:
        f.create_dataset("ca", data=coords, compression="gzip", compression_opts=4)
        f.create_dataset("res_ids", data=_WORKER["res_ids"])
        f.create_dataset("cluster", data=labels)
        f.create_dataset("weight", data=weights)
        f.attrs["cond"] = a.cond
        f.attrs["max_frames"] = a.max_frames
        f.attrs["seed"] = a.seed
        f.attrs["n_missing_segs"] = n_missing
    uniq, cnts = np.unique(labels, return_counts=True)
    print(f"[{a.cond}] wrote {coords.shape} to {a.out} "
          f"({n_missing} missing segs); per-cluster: "
          + ", ".join(f"c{u}={n}" for u, n in zip(uniq, cnts)), flush=True)


if __name__ == "__main__":
    main()
