#!/bin/bash
#SBATCH --job-name=collect_coords
#SBATCH --output=collect.out
#SBATCH --error=collect.err
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cluster=invest
#SBATCH --partition=lchong
##SBATCH --gres=gpu:4
#SBATCH --time=47:59:59


source ~/.bashrc
conda activate hamsm
python collectCoords.py west.h5 1x88_dry.pdb ../
