#!/bin/bash
#SBATCH --job-name=vv00_1x88_PROD_10
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cluster=gpu 
#SBATCH --gres=gpu:1
#SBATCH --partition=gtx1080 
#SBATCH --time=144:00:00 
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_prod_10.out 
#SBATCH --error=slurm_prod_10.err 

# load appropriate modules.  First purge all modules then load in intel (a prereq for loading in amber) and then amber 
module purge
module load intel/2017.3.196
module load amber/18

# echo commands to stdout
set -x

# Executable
SANDER=pmemd.cuda

# Amber input files and output name
INP=10_prod.in
TOP=1x88_solv.prmtop
CRD=09_prod.rst
OUT=10_prod

$SANDER  -O     -i   $INP   -p   $TOP   -c   $CRD   -r   $OUT.rst \
        -o   $OUT.out   -e   $OUT.ene   -v   $OUT.vel   -inf $OUT.nfo   -x   $OUT.nc &&

# make cpptraj input file
echo "parm 1x88_solv.prmtop"                            > strip_10.in
echo "trajin 10_prod.nc"                                    >> strip_10.in
echo "autoimage"                                            >> strip_10.in
echo "strip :WAT,Cl-,Na+"                                   >> strip_10.in
echo "trajout 10_prod_dry.nc"                               >> strip_10.in
echo "go"                                                   >> strip_10.in
echo "quit"                                                 >> strip_10.in

# run strip_10 to strip water
cpptraj -i strip_10.in > strip_10.out

# remove solvated nc file
if [ -f "10_prod_dry.nc" ] ; then
    rm 10_prod.nc &&
    mv 10_prod_dry.nc 10_prod.nc
fi

