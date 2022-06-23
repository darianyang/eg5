#!/bin/bash
#SBATCH --job-name=vv00_1x88_apo_PROD_Multi
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cluster=gpu 
#SBATCH --gres=gpu:1
#SBATCH --partition=gtx1080 
#SBATCH --time=144:00:00 
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_prod_06.out 
#SBATCH --error=slurm_prod_06.err 

# load appropriate modules.  First purge all modules then load in intel (a prereq for loading in amber) and then amber 
module purge
module load intel/2017.3.196
module load amber/18

# echo commands to stdout
set -x

# set variables
#export DO_PARALLEL="mpirun -np 16"
SANDER=pmemd.cuda

###################################################
###################################################
###################################################
PDB="1x88_apo_solv" 
###################################################
###################################################
###################################################

for PROD in {07..10} ; do
    cp -v 06_prod.in ${PROD}_prod.in
    # keep leading zero but prevent octal interpretation, keep base10
    PREV=$(printf "%02d" $((10#$PROD - 1)))
    # 1 ns pcs restrained NPT production Langevin thermostat (298K) MC barostat (1atm)
    printf "\nstart $PROD prod\n\n" 
    #$DO_PARALLEL \
        $SANDER -O \
                -i ${PROD}_prod.in -o ${PROD}_prod.out \
                -p ${PDB}.prmtop -c ${PREV}_prod.rst \
                -r ${PROD}_prod.rst -x ${PROD}_prod.nc \
                -inf ${PROD}_prod.nfo &&
    printf "\ndone $PROD prod\n\n"
done
