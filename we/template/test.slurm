#!/bin/bash
#SBATCH --job-name=eg5_holo
#SBATCH --output=job_logs/slurm.out
#SBATCH --error=job_logs/slurm.err
#SBATCH --mail-user=atb43@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=GPU
#SBATCH --gpus=8
#SBATCH --time=00:00:59

echo "hello" > test.out

