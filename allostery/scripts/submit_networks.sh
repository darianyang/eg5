#!/bin/bash
# Build all per-(condition, cluster, scheme) allosteric networks from the merged
# movement files, then run the cross-condition / cross-cluster comparison.
# Run after submit_extract.sh + merge_shards.py have produced
# data/movements_<cond>.h5.
#
#     sbatch submit_networks.sh
#SBATCH --job-name=eg5_allo_networks
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --account=lchong
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=logs/networks_%j.out
#SBATCH --error=logs/networks_%j.out

set -e
cd "$SLURM_SUBMIT_DIR"

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE

source "$HOME/.bashrc"
conda activate eg5

python scripts/build_networks.py \
    --conds nomon wmon \
    --schemes pooled weighted \
    --pdb data/ref_nomon.pdb \
    --ncpu "${SLURM_CPUS_PER_TASK:-16}"

python scripts/compare_networks.py
echo "done"
