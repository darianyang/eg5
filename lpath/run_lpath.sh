#!/bin/bash
# LPATH pathway analysis for the two Eg5 ADP-unbinding WE simulations, using a
# SHARED clustering space (see cluster_shared.py, run that first).
#
#   nomon = WT Eg5 (no monastrol)     -> multi-mab_nomon_v01
#   wmon  = Eg5 + monastrol           -> multi-mab_wmon_v00
#
# Source/target states come from the eg5_poster FULLY_UNBOUND2 w_assign scheme
# in each run's west.cfg:  state 0 = bound (source), state 1 = unbound (target).
# The shared cluster labels (auxdata/labels) are used as the states for matching.
set -e

# headless plotting
export MPLBACKEND=Agg
# /ix is NFS; HDF5 file locking is unreliable there
export HDF5_USE_FILE_LOCKING=FALSE

SCHEME=FULLY_UNBOUND2

for RUN in nomon wmon; do
    echo "=================================================================="
    echo " LPATH: $RUN"
    echo "=================================================================="
    pushd "$RUN" > /dev/null

    # reassign_custom must be importable from the cwd
    cp ../reassign_custom.py .
    mkdir -p plots

    # w_assign writes into the analysis directory named in west.cfg (WIPA)
    ASSIGN=./WIPA/$SCHEME/assign.h5

    # 1. discretize: assign source/target states with w_assign + FULLY_UNBOUND2
    #    (skip if the assign.h5 already exists -- w_assign is the slow step)
    if [ ! -f "$ASSIGN" ]; then
        lpath discretize -we -W ./west.h5 \
            --assign-arguments="--config-from-file --scheme $SCHEME -W west.h5"
    else
        echo "  $ASSIGN exists, skipping discretize"
    fi

    # 2. extract successful pathways bound(0) -> unbound(1); carry shared labels
    lpath extract -we -W ./west.h5 -A "$ASSIGN" \
        -ss 0 -ts 1 -p -a labels

    # 3. match; states = shared cluster labels. condense=1 strips only
    #    consecutive-duplicate dwell frames. --n-clusters bypasses the
    #    interactive "how many clusters?" prompt; `printf 'n\n'` answers the
    #    "regenerate dendrogram?" prompt (avoids timedinput EOFError when
    #    stdin is not a TTY).
    printf 'n\n' | lpath match -we -ra reassign_custom.reassign_custom \
        -op succ_traj/reassigned.pickle --condense 1 --n-clusters 6 --timeout 5 \
        --plots-hide --plot-out-path plots

    popd > /dev/null
done

echo "Done. See nomon/plots and wmon/plots for dendrograms."
