#!/bin/bash
# LPATH pathway analysis for the two Eg5 ADP-unbinding WE simulations, using a
# SHARED clustering space (run cluster_shared.py first to write auxdata/labels).
#
#   nomon = WT Eg5 (no monastrol)  -> multi-mab_nomon_v01
#   wmon  = Eg5 + monastrol        -> multi-mab_wmon_v00
#
# Source/target come from the eg5_poster FULLY_UNBOUND2 w_assign scheme in each
# west.cfg: state 0 = bound (source), state 1 = unbound (target). The shared
# cluster labels are used as the states for matching (reassign_custom).
#
#     sbatch submit_lpath.sh
#SBATCH --job-name=eg5_lpath
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --account=lchong
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/lpath_%j.out
#SBATCH --error=logs/lpath_%j.out

set -e
cd "$SLURM_SUBMIT_DIR"

# headless plotting; /ix is NFS so disable HDF5 file locking; keep BLAS single-threaded
export MPLBACKEND=Agg
export HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

source "$HOME/.bashrc"
conda activate eg5

SCHEME=FULLY_UNBOUND2

for RUN in nomon wmon; do
    echo "=================================================================="
    echo " LPATH: $RUN   ($(date))"
    echo "=================================================================="
    pushd "$RUN" > /dev/null

    cp ../reassign_custom.py .            # must be importable from cwd
    mkdir -p plots
    ASSIGN=./WIPA/$SCHEME/assign.h5       # w_assign writes into west.cfg analysis dir

    # 1. discretize: assign source/target states with w_assign + FULLY_UNBOUND2
    #    (skip if assign.h5 already exists -- w_assign is the slow step)
    if [ ! -f "$ASSIGN" ]; then
        lpath discretize -we -W ./west.h5 \
            --assign-arguments="--config-from-file --scheme $SCHEME -W west.h5"
    else
        echo "  $ASSIGN exists, skipping discretize"
    fi

    # 2. extract successful pathways bound(0) -> unbound(1); carry shared labels
    #    (skip if output.pickle already exists)
    if [ ! -f succ_traj/output.pickle ]; then
        lpath extract -we -W ./west.h5 -A "$ASSIGN" -ss 0 -ts 1 -p -a labels
    else
        echo "  succ_traj/output.pickle exists, skipping extract"
    fi

    # 3. match; states = shared cluster labels. condense=1 strips only
    #    consecutive-duplicate dwell frames (condense 2 over-flattened the
    #    paths into ~3 strings). --n-clusters bypasses the interactive
    #    "how many clusters?" prompt; `printf 'n\n'` answers the remaining
    #    "regenerate dendrogram?" prompt (batch stdin is /dev/null, which
    #    otherwise makes timedinput raise EOFError and crash). Do NOT pipe
    #    `yes n` without --n-clusters: it feeds "n" into the number prompt,
    #    which loops forever on Invalid input.
    printf 'n\n' | lpath match -we -ra reassign_custom.reassign_custom \
        -op succ_traj/reassigned.pickle --condense 1 --n-clusters 6 --timeout 5 \
        --plots-hide --plot-out-path plots

    popd > /dev/null
done

echo "Done ($(date)). See nomon/plots and wmon/plots for dendrograms."
