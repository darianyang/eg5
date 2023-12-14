#!/bin/bash
#SBATCH --job-name=Solv_Eg5_DL5
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=4
#SBATCH --cluster=invest
#SBATCH --partition=lchong
#SBATCH --time=23:59:59 
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_solv.out 
#SBATCH --error=slurm_solv.err 

conda activate we
tleap -f tleap1.in > tleap1.out
