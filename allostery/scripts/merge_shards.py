"""Merge per-iteration extraction shards into one movements file per condition.

Concatenates ``data/shards/<cond>_*.h5`` (written by submit_extract.sh) into
``data/movements_<cond>.h5`` with the same dataset layout, and prints a per
-cluster summary: movement counts, boundary losses, and the Kish effective
sample size (ESS = (sum w)^2 / sum w^2) of the within-cluster-normalized WE
weights -- the key diagnostic for whether a WE-weighted network is trustworthy
or dominated by a few high-weight walkers.
"""
import os
import sys
import glob
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(__file__))
from eg5_allostery import N_CLUSTERS  # noqa: E402


def merge(cond, shard_dir="data/shards", out_dir="data"):
    shards = sorted(glob.glob(f"{shard_dir}/{cond}_*.h5"))
    if not shards:
        raise SystemExit(f"no shards for {cond} in {shard_dir}")

    with h5py.File(shards[0], "r") as f0:
        angles = f0.attrs["angles"].split(",")
        res_ids = {a: f0[f"res_ids_{a}"][:] for a in angles}

    mov = {a: [] for a in angles}
    weight, cluster = [], []
    n_missing = n_segments = 0
    for s in shards:
        with h5py.File(s, "r") as f:
            for a in angles:
                mov[a].append(f[f"movements_{a}"][:])
            weight.append(f["weight"][:])
            cluster.append(f["cluster"][:])
            n_missing += int(f.attrs["n_missing_segs"])
            n_segments += int(f.attrs["n_segments"])

    weight = np.concatenate(weight)
    cluster = np.concatenate(cluster)

    out = f"{out_dir}/movements_{cond}.h5"
    with h5py.File(out, "w") as f:
        for a in angles:
            f.create_dataset(f"movements_{a}", data=np.concatenate(mov[a]),
                             compression="gzip", compression_opts=4)
            f.create_dataset(f"res_ids_{a}", data=res_ids[a])
        f.create_dataset("weight", data=weight)
        f.create_dataset("cluster", data=cluster)
        f.attrs["cond"] = cond
        f.attrs["angles"] = ",".join(angles)
        f.attrs["n_missing_segs"] = n_missing
        f.attrs["n_segments"] = n_segments

    kept = cluster >= 0
    print(f"\n{cond}: {cluster.size} movements from {n_segments} segments "
          f"({n_missing} missing), {int(kept.sum())} in-cluster, "
          f"{int((~kept).sum())} boundary ({(~kept).mean():.1%})")
    print(f"  angles: {angles}  cols: " +
          ", ".join(f"{a}={res_ids[a].size}" for a in angles))
    print(f"  {'clust':>5} {'n_move':>9} {'ESS':>9} {'ESS/n':>7} {'wmax_frac':>9}")
    for c in range(N_CLUSTERS):
        m = cluster == c
        n = int(m.sum())
        if n == 0:
            print(f"  {c:>5} {0:>9}")
            continue
        w = weight[m]
        w = w / w.sum()
        ess = 1.0 / np.sum(w ** 2)
        print(f"  {c:>5} {n:>9} {ess:>9.1f} {ess / n:>7.3f} {w.max():>9.3f}")
    return out


if __name__ == "__main__":
    conds = sys.argv[1:] or ["nomon", "wmon"]
    for c in conds:
        merge(c)
