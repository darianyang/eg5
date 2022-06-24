#!/bin/bash
#SBATCH --job-name=1ii6_PREP
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8
#SBATCH --cluster=invest
#SBATCH --partition=lchong
#SBATCH --time=23:59:59 
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_prep.out 
#SBATCH --error=slurm_prep.err 

# load appropriate modules.  First purge all modules then load in intel (a prereq for loading in amber) and then amber 
module purge
module load intel/2017.3.196
module load amber/18

# echo commands to stdout
set -x

# set variables
export DO_PARALLEL="mpirun -np 8"
SANDER=pmemd.MPI

###################################################
###################################################
###################################################
PDB="1ii6_solv" 
###################################################
###################################################
###################################################

# 10000 step minimization (500 steepest descent)
printf "\nstart min\n\n" 
$DO_PARALLEL \
    $SANDER -O \
            -i 02_min.in -o 02_min.out \
            -p ${PDB}.prmtop -c ${PDB}.inpcrd -r 02_min.rst &&
printf "\ndone min\n\n" 
