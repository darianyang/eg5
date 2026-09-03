#!/bin/bash
# Extract per-cluster Ca coordinate ensembles for DF / SPM analysis, one job per
# condition.  Submit once per condition:
#
#     sbatch --export=COND=nomon scripts/submit_ca_coords.sh
#     sbatch --export=COND=wmon  scripts/submit_ca_coords.sh
#
# Output: data/ca_coords_<cond>.h5  (ca, res_ids, cluster, weight)
#SBATCH --job-name=eg5_ca_coords
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --account=lchong
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ca_coords_%x_%j.out
#SBATCH --error=logs/ca_coords_%x_%j.out

set -e
cd "$SLURM_SUBMIT_DIR"

: "${COND:?set COND=nomon or COND=wmon via --export}"
: "${MAXFRAMES:=20000}"

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE

source "$HOME/.bashrc"
conda activate eg5

mkdir -p data logs
OUT="data/ca_coords_${COND}.h5"

echo "COND=$COND max_frames=$MAXFRAMES -> $OUT on $(hostname), ${SLURM_CPUS_PER_TASK} cpus"
python scripts/extract_ca_coords.py "$COND" "$OUT" \
    --max-frames "$MAXFRAMES" --ncpu "${SLURM_CPUS_PER_TASK:-16}"
echo "done $COND"
