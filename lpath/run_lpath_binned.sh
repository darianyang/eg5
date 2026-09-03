#!/bin/bash
# LPATH match step re-run with the SHARED FIXED BIN GRID discretization
# (build_shared_grid.py + reassign_binned.py) instead of the 6 agglomerative
# shared clusters (cluster_shared.py + reassign_custom.py).
#
#   nomon = WT Eg5 (no monastrol)   wmon = Eg5 + monastrol
#
# The slow discretize (w_assign) and extract steps are UNCHANGED and reused from
# the previous run: this reads each run's existing succ_traj/output.pickle (which
# already carries the pcoords and w_assign source/sink states) and only re-runs
# `match` with -ra reassign_binned.  Outputs use a *_binned suffix so the earlier
# 6-cluster results are preserved for comparison.
#
# Source/sink (FULLY_UNBOUND2, unchanged): state 0 = bound, state 1 = unbound.
# Intermediate frames are placed on the shared grid; bound/unbound stay as
# dedicated anchor states.  Run cluster_paths_shared_binned.py afterwards for the
# comparable shared pathway dendrogram.
set -e

export MPLBACKEND=Agg
export HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

N_CLUSTERS=${1:-6}      # per-run dendrogram cut (shared analysis re-cuts anyway)

# 0. build the shared fixed grid + dictionary once (reads both output.pickles)
echo "== building shared grid =="
python build_shared_grid.py

for RUN in nomon wmon; do
    echo "=================================================================="
    echo " LPATH match (binned grid): $RUN"
    echo "=================================================================="
    pushd "$RUN" > /dev/null

    cp ../reassign_binned.py .          # must be importable from cwd
    mkdir -p plots_binned

    if [ ! -f succ_traj/output.pickle ]; then
        echo "  ERROR: succ_traj/output.pickle missing -- run the extract step first" >&2
        exit 1
    fi

    printf 'n\n' | lpath match -we -ra reassign_binned.reassign_binned \
        -ip succ_traj/output.pickle \
        -op succ_traj/reassigned_binned.pickle \
        -dF succ_traj/distmat_binned.npy \
        -co succ_traj/cluster_labels_binned.npy \
        --condense 1 --n-clusters "$N_CLUSTERS" --timeout 5 \
        --plots-hide --plot-out-path plots_binned

    popd > /dev/null
done

echo "== shared pathway dendrogram (binned) =="
python cluster_paths_shared_binned.py

echo "Done. See {nomon,wmon}/plots_binned and shared_paths_binned/."
