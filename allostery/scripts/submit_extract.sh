#!/bin/bash
# Extract WE dihedral movements (phi, chi1, chi2) as HDF5 shards, one SLURM
# array task per contiguous block of iterations.  Submit once per condition:
#
#     sbatch --export=COND=nomon submit_extract.sh
#     sbatch --export=COND=wmon  submit_extract.sh
#
# then merge the shards with merge_shards.py.
#SBATCH --job-name=eg5_allo_extract
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --account=lchong
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH --array=0-15
#SBATCH --output=logs/extract_%x_%A_%a.out
#SBATCH --error=logs/extract_%x_%A_%a.out

set -e
cd "$SLURM_SUBMIT_DIR"

: "${COND:?set COND=nomon or COND=wmon via --export}"

# total iterations per condition (nomon: 586, wmon: 421)
case "$COND" in
    nomon) NITER=586 ;;
    wmon)  NITER=421 ;;
    *) echo "unknown COND=$COND"; exit 1 ;;
esac

NTASK=16
TID=${SLURM_ARRAY_TASK_ID:-0}

# contiguous iteration block [lo, hi] for this task (1-indexed, inclusive)
CHUNK=$(( (NITER + NTASK - 1) / NTASK ))
LO=$(( TID * CHUNK + 1 ))
HI=$(( LO + CHUNK - 1 ))
if [ "$HI" -gt "$NITER" ]; then HI=$NITER; fi
if [ "$LO" -gt "$NITER" ]; then echo "task $TID: no iterations, exiting"; exit 0; fi

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE

source "$HOME/.bashrc"
conda activate eg5

mkdir -p data/shards logs
OUT="data/shards/${COND}_$(printf '%02d' "$TID").h5"

echo "COND=$COND task=$TID iters ${LO}-${HI} -> $OUT on $(hostname), ${SLURM_CPUS_PER_TASK} cpus"
python scripts/eg5_allostery.py "$COND" "$LO" "$HI" "$OUT" \
    --angles phi,chi1,chi2 --ncpu "${SLURM_CPUS_PER_TASK:-16}"
echo "done task $TID"
