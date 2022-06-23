#!/bin/bash
#SBATCH --job-name=vv00_1x88_PROD_07
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cluster=gpu 
#SBATCH --gres=gpu:1
#SBATCH --partition=gtx1080 
#SBATCH --time=144:00:00 
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_prod_07.out 
#SBATCH --error=slurm_prod_07.err 

# load appropriate modules.  First purge all modules then load in intel (a prereq for loading in amber) and then amber 
module purge
module load intel/2017.3.196
module load amber/18

# echo commands to stdout
set -x

# Executable
SANDER=pmemd.cuda

# Amber input files and output name
INP=07_prod.in
TOP=1x88_solv.prmtop
CRD=06_prod.rst
OUT=07_prod

$SANDER  -O     -i   $INP   -p   $TOP   -c   $CRD   -r   $OUT.rst \
        -o   $OUT.out   -e   $OUT.ene   -v   $OUT.vel   -inf $OUT.nfo   -x   $OUT.nc &&

# make strip_07.input file
echo "parm 1x88_solv.prmtop"                            > strip_07.in
echo "trajin 07_prod.nc"                                    >> strip_07.in
echo "autoimage"                                            >> strip_07.in
echo "strip :WAT,Cl-,Na+"                                   >> strip_07.in
echo "trajout 07_prod_dry.nc"                               >> strip_07.in
echo "go"                                                   >> strip_07.in
echo "quit"                                                 >> strip_07.in

# run cpptraj to strip water
cpptraj -i strip_07.in > strip_07.out

# remove solvated nc file
if [ -f "07_prod_dry.nc" ] ; then
    rm 07_prod.nc &&
    mv 07_prod_dry.nc 07_prod.nc
fi

sbatch h2p_1gpu_prod_08.slurm
