#!/bin/bash
#SBATCH --job-name=trace_traj_gen_v00_1x88
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=2000
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 24:00:00
#SBATCH --ntasks-per-node=8
#SBATCH -A bio210074p

# source for conda init
source ~/.bashrc 
# note that md analysis and amber parmed conflict 
# so can't source amber for this to work
conda activate md
python combine_trace.py
