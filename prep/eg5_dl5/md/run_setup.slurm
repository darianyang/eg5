#!/bin/bash
#SBATCH --job-name=setup_equil_Eg5_DL5
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=4
#SBATCH --cluster=invest
#SBATCH --partition=lchong
#SBATCH --time=0:04:59 
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_setup.out 
#SBATCH --error=slurm_setup.err 

bash scripts/multi_set_and_run.sh
