#!/bin/bash
#SBATCH --job-name=1ii6_dl5_PREP_vv00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=gtx1080
#SBATCH --time=4:59:59
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_prep.out 
#SBATCH --error=slurm_prep.err 

# load appropriate modules.  First purge all modules then load in intel (a prereq for loading in amber) and then amber 
module purge
#module load intel/2017.3.196
#module load amber/18
module load gcc/8.2.0
module load openmpi/4.0.3
module load amber/22

# echo commands to stdout
set -x

# set variables
#export DO_PARALLEL="mpirun -np 8"
#SANDER=pmemd.MPI
SANDER=pmemd.cuda

###################################################
###################################################
###################################################
PDB="1ii6_dl5_solv" 
###################################################
###################################################
###################################################

## 10000 step minimization (500 steepest descent)
#printf "\nstart min\n\n" 
#$DO_PARALLEL \
#    $SANDER -O \
#            -i 02_min.in -o 02_min.out \
#            -p ${PDB}.prmtop -c ${PDB}.inpcrd -r 02_min.rst &&
#printf "\ndone min\n\n" 
#
## 20ps restrained NVT equilibration with Langevin thermostat (298K)
#printf "\nstart restrained system heat in\n\n" 
#$DO_PARALLEL \
#    $SANDER -O \
#            -i 03_eq1.in -o 03_eq1.out \
#	    -x 03_eq1.nc \
#            -p ${PDB}.prmtop -c 02_min.rst -r 03_eq1.rst -ref 02_min.rst &&
#printf "\ndone restrained system heat in\n\n" 

# 1 ns restrained NPT equilibration Langevin thermostat (298K) MC barostat (1atm)
#printf "\nstart restrained equil\n\n" 
##$DO_PARALLEL \
#    $SANDER  -O \
#            -i 04_eq2.in -o 04_eq2.out \
#	    -x 04_eq2.nc \
#            -p ${PDB}.prmtop -c 03_eq1.rst -r 04_eq2.rst -ref 03_eq1.rst &&
#printf "\ndone restrained equil\n\n" 

# 1 ns unrestrained NPT equilibration Langevin thermostat (298K) MC barostat (1atm)
printf "\nstart unrestrained equil\n\n" 
#$DO_PARALLEL \
    $SANDER -O \
            -i 05_eq3.in -o 05_eq3.out \
	    -x 05_eq3.nc \
            -p ${PDB}.prmtop -c 04_eq2.rst -r 05_eq3.ncrst &&
printf "\ndone unrestrained equil\n\n"

printf "\nSUBMITTING PROD RUN\n" &&

#sbatch h2p_1gpu_prod_06.slurm
