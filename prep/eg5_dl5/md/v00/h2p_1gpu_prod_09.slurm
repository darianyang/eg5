#!/bin/bash
#SBATCH --job-name=vv00_1ii6_dl5_PROD_09
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cluster=gpu 
#SBATCH --gres=gpu:1
#SBATCH --partition=gtx1080 
#SBATCH --time=144:00:00 
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_prod_09.out 
#SBATCH --error=slurm_prod_09.err 

# load appropriate modules.  First purge all modules then load in intel (a prereq for loading in amber) and then amber 
module purge
module load intel/2017.3.196
module load amber/18

# echo commands to stdout
set -x

# Executable
SANDER=pmemd.cuda

# Amber input files and output name
INP=09_prod.in
TOP=1ii6_dl5_solv.prmtop
CRD=08_prod.rst
OUT=09_prod

$SANDER  -O     -i   $INP   -p   $TOP   -c   $CRD   -r   $OUT.rst \
        -o   $OUT.out   -e   $OUT.ene   -v   $OUT.vel   -inf $OUT.nfo   -x   $OUT.nc &&

# make cpptraj input file
echo "parm 1ii6_dl5_solv.prmtop"                            > strip_09.in
echo "trajin 09_prod.nc"                                    >> strip_09.in
echo "autoimage"                                            >> strip_09.in
echo "strip :WAT,Cl-,Na+"                                   >> strip_09.in
echo "trajout 09_prod_dry.nc"                               >> strip_09.in
echo "go"                                                   >> strip_09.in
echo "quit"                                                 >> strip_09.in

# run strip_09 to strip water
strip_09 -i strip_09.in > strip_09.out

# remove solvated nc file
if [ -f "09_prod_dry.nc" ] ; then
    rm 09_prod.nc &&
    mv 09_prod_dry.nc 09_prod.nc
fi

sbatch h2p_1gpu_prod_10.slurm
