"""Shared helpers for LPATH-cluster-resolved allosteric network analysis of the
two Eg5 ADP-unbinding WE simulations (nomon = WT, wmon = +monastrol).

The idea: each WE frame already carries a shared 6-cluster label (0..5) in
``auxdata/labels`` (see ../../lpath/cluster_shared.py).  Those clusters are
stages along the ADP-exit coordinate, so a per-cluster mdpath network gives a
path/reaction-coordinate-resolved allostery comparison.

mdpath's signal is the *movement* of a residue dihedral between consecutive
frames.  In a WE run the only genuinely continuous timeseries is the set of
frames inside a single ``seg.nc`` (5 frames, uniform 10 ps spacing), so we
compute wrapped dihedral differences ONLY within a segment and never across a
segment/file boundary.  We then pool the resulting *movements* (not raw frames)
per cluster, exactly as mdpath pools movements across replicas.

A movement is assigned to cluster ``c`` only when BOTH of its endpoint frames
carry label ``c`` -- i.e. it is a genuine within-state transition.  Each movement
inherits its segment's WE weight, so a WE-weighted network is a weighted
histogram over movements (weights normalized within the cluster at build time).

Angles: we extract backbone ``phi`` (what mdpath uses today) plus sidechain
``chi1`` and ``chi2`` movements, stored as separate angle-type datasets that
share the same rows (movements), weights and cluster vector.  phi drives the
current mdpath networks; chi1/chi2 are kept so sidechain dynamics can be folded
into the allostery analysis later.

Residue numbering: residues are labelled in AMBER numbering
(``amber = mdtraj_resSeq + 1``), so protein residues are 1..368 and phi is
defined for 2..368, matching the ``:1-368`` selection in the WE runseg scripts.
"""
import os
import numpy as np
import h5py
import mdtraj as md

# keep BLAS single-threaded; we parallelize over segments ourselves
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

ROOT = "/ix/lchong/dty7/eg5"

# per-condition inputs: labelled west.h5 (shared cluster labels live here), the
# water-stripped topology that matches seg.nc, and the traj_segs tree.
CONDITIONS = {
    "nomon": {
        "westh5":  f"{ROOT}/lpath/nomon/west.h5",
        "prmtop":  f"{ROOT}/h2p-we/multi-mab_nomon_v01/common_files/eg5_2022_dry.prmtop",
        "trajseg": f"{ROOT}/h2p-we/multi-mab_nomon_v01/traj_segs",
    },
    "wmon": {
        "westh5":  f"{ROOT}/lpath/wmon/west.h5",
        "prmtop":  f"{ROOT}/h2p-we/multi-mab_wmon_v00/common_files/1x88_dry.prmtop",
        "trajseg": f"{ROOT}/h2p-we/multi-mab_wmon_v00/traj_segs",
    },
}

N_CLUSTERS = 6
RAD2DEG = 180.0 / np.pi

# angle types to extract.  ``fn`` is the mdtraj compute function; ``res_atom`` is
# the position within each returned atom-index quartet whose residue the angle is
# assigned to (phi spans residues i-1..i, its CA is atom 2 -> residue i; chi
# angles lie within one residue, so atom 0 suffices).
ANGLE_SPECS = {
    "phi":  {"fn": md.compute_phi,  "res_atom": 2},
    "chi1": {"fn": md.compute_chi1, "res_atom": 0},
    "chi2": {"fn": md.compute_chi2, "res_atom": 0},
}


def angle_residue_ids(prmtop, angle):
    """AMBER residue ids (resSeq + 1) each angle of ``angle`` type belongs to.

    Aligned with the columns returned by the corresponding mdtraj compute
    function.  chi angles are defined only for residues that have that sidechain
    torsion, so the returned list is a subset of the protein residues.
    """
    top = md.load_prmtop(prmtop)
    spec = ANGLE_SPECS[angle]
    idx, _ = spec["fn"](md.Trajectory(np.zeros((1, top.n_atoms, 3)), top))
    ra = spec["res_atom"]
    return [top.atom(a[ra]).residue.resSeq + 1 for a in idx]


def _wrapped_deg_diff(ang_rad):
    """Consecutive wrapped movements in degrees for one segment.

    ``ang_rad`` is (n_frames, n_angles).  Returns (n_frames-1, n_angles) with
    each difference wrapped onto [-180, 180) (dihedrals are periodic; an
    unwrapped step across the branch cut would fake a ~360 deg jump).  Matches
    mdpath's DihedralAngles wrapping convention.
    """
    deg = ang_rad * RAD2DEG
    raw = np.diff(deg, axis=0)
    return (raw + 180.0) % 360.0 - 180.0


def _load_iter_meta(westh5):
    """Return {iter_key: (weights[n_seg], labels[n_seg, pcoord_len])} for a run."""
    meta = {}
    with h5py.File(westh5, "r") as f:
        for k in sorted(f["iterations"].keys()):
            w = f[f"iterations/{k}/seg_index"]["weight"][:]
            lab = f[f"iterations/{k}/auxdata/labels"][:][..., 0].astype(np.int16)
            meta[k] = (w, lab)
    return meta


def _segment_movements(seg_nc, top, labels_seg, angles):
    """Per-angle movements and endpoint cluster labels for one segment.

    Returns ``(mov_by_angle, lab_lo, lab_hi)`` where ``mov_by_angle`` maps each
    angle name to an (n_move, n_angles) array of wrapped movements (deg), and
    ``lab_lo``/``lab_hi`` are the endpoint cluster labels of each movement.
    seg.nc holds the dynamics frames (pcoord points 1..len-1), so labels are
    aligned to the last ``n_frames`` entries.
    """
    tr = md.load(seg_nc, top=top)
    if tr.n_frames < 2:
        return None
    mov_by_angle = {}
    for a in angles:
        _, ang = ANGLE_SPECS[a]["fn"](tr)
        mov_by_angle[a] = _wrapped_deg_diff(ang).astype(np.float32)
    lab = labels_seg[-tr.n_frames:]
    return mov_by_angle, lab[:-1], lab[1:]


_WORKER = {}


def _init_worker(prmtop, angles):
    """Pool initializer: load the topology once per worker."""
    _WORKER["top"] = md.load_prmtop(prmtop)
    _WORKER["angles"] = angles


def _worker_segment(task):
    """Pool task: (mov_by_angle, weight_vec, cluster_vec) for one segment."""
    seg_nc, labels_seg, weight = task
    if not os.path.exists(seg_nc):
        return None
    try:
        res = _segment_movements(seg_nc, _WORKER["top"], labels_seg, _WORKER["angles"])
    except Exception:
        return None
    if res is None:
        return None
    mov_by_angle, lo, hi = res
    n = lo.shape[0]
    cl = np.where(lo == hi, lo, -1).astype(np.int8)
    wt = np.full(n, weight, dtype=np.float64)
    return mov_by_angle, wt, cl


def extract_shard(cond, iter_lo, iter_hi, out_path, angles=("phi", "chi1", "chi2"),
                  ncpu=1, verbose=True):
    """Extract within-segment dihedral movements for a range of iterations.

    Writes an HDF5 shard with, per angle type ``A``:
        movements_A (n_move, n_A) float32 -- wrapped movements (deg)
        res_ids_A   (n_A,)        int32   -- AMBER residue id per column
    and shared:
        weight  (n_move,) float64 -- WE walker weight of the source segment
        cluster (n_move,) int8    -- shared cluster id, or -1 for a boundary
                                     movement whose endpoint frames differ
    All angle movement arrays share the same rows, ``weight`` and ``cluster``.
    """
    angles = list(angles)
    info = CONDITIONS[cond]
    res_ids = {a: np.asarray(angle_residue_ids(info["prmtop"], a), dtype=np.int32)
               for a in angles}
    meta = _load_iter_meta(info["westh5"])

    keys = sorted(k for k in meta if iter_lo <= int(k.split("_")[-1]) <= iter_hi)

    tasks = []   # (seg_nc, labels_seg, weight)
    for k in keys:
        it = int(k.split("_")[-1])
        w, lab = meta[k]
        for seg in range(w.shape[0]):
            seg_nc = f"{info['trajseg']}/{it:06d}/{seg:06d}/seg.nc"
            tasks.append((seg_nc, lab[seg], float(w[seg])))

    mov_all = {a: [] for a in angles}
    wt_all, cl_all = [], []
    n_missing = 0

    def _collect(i, r):
        nonlocal n_missing
        if r is None:
            n_missing += 1
            return
        mov_by_angle, wt, cl = r
        for a in angles:
            mov_all[a].append(mov_by_angle[a])
        wt_all.append(wt); cl_all.append(cl)
        if verbose and i % 5000 == 0:
            print(f"  {cond} {iter_lo}-{iter_hi}: {i}/{len(tasks)} segs", flush=True)

    if ncpu > 1:
        from multiprocessing import Pool
        with Pool(ncpu, initializer=_init_worker,
                  initargs=(info["prmtop"], angles)) as pool:
            for i, r in enumerate(pool.imap(_worker_segment, tasks, chunksize=32)):
                _collect(i, r)
    else:
        _init_worker(info["prmtop"], angles)
        for i, task in enumerate(tasks):
            _collect(i, _worker_segment(task))

    weight = np.concatenate(wt_all) if wt_all else np.zeros(0, np.float64)
    cluster = np.concatenate(cl_all) if cl_all else np.zeros(0, np.int8)

    with h5py.File(out_path, "w") as f:
        for a in angles:
            m = (np.concatenate(mov_all[a]) if mov_all[a]
                 else np.zeros((0, res_ids[a].size), np.float32))
            f.create_dataset(f"movements_{a}", data=m,
                             compression="gzip", compression_opts=4)
            f.create_dataset(f"res_ids_{a}", data=res_ids[a])
        f.create_dataset("weight", data=weight)
        f.create_dataset("cluster", data=cluster)
        f.attrs["cond"] = cond
        f.attrs["angles"] = ",".join(angles)
        f.attrs["iter_lo"] = iter_lo
        f.attrs["iter_hi"] = iter_hi
        f.attrs["n_missing_segs"] = n_missing
        f.attrs["n_segments"] = len(tasks)
    if verbose:
        kept = int((cluster >= 0).sum())
        print(f"  {cond} {iter_lo}-{iter_hi}: {cluster.size} movements "
              f"({kept} in-cluster, {cluster.size-kept} boundary), "
              f"{n_missing} missing segs -> {out_path}", flush=True)
    return out_path


def build_reference_pdb(cond, out_pdb):
    """Write a protein-only reference PDB (AMBER residue numbering) for the graph.

    Uses the first dynamics frame of iteration 1, segment 0 -- a bound-state
    conformation.  A single shared reference graph (same residues, same 5 A
    proximity edges) is used for every cluster/condition so that only the NMI
    edge weights differ between networks.
    """
    info = CONDITIONS[cond]
    top = md.load_prmtop(info["prmtop"])
    seg_nc = f"{info['trajseg']}/000001/000000/seg.nc"
    tr = md.load(seg_nc, top=top)[0]
    prot = tr.atom_slice(tr.top.select("protein"))
    for res in prot.top.residues:              # renumber to AMBER (resSeq + 1)
        res.resSeq = res.resSeq + 1
    prot.save_pdb(out_pdb)
    return out_pdb


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Extract one WE dihedral-movement shard.")
    p.add_argument("cond", choices=list(CONDITIONS))
    p.add_argument("iter_lo", type=int)
    p.add_argument("iter_hi", type=int)
    p.add_argument("out_path")
    p.add_argument("--angles", default="phi,chi1,chi2")
    p.add_argument("--ncpu", type=int, default=1)
    a = p.parse_args()
    extract_shard(a.cond, a.iter_lo, a.iter_hi, a.out_path,
                  angles=tuple(a.angles.split(",")), ncpu=a.ncpu)
