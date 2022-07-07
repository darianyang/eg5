#!/bin/bash
#SBATCH --job-name=vVER_PDB_TEMP_PROD_08
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --cluster=gpu 
#SBATCH --gres=gpu:1
#SBATCH --partition=gtx1080 
#SBATCH --time=144:00:00 
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_prod_08.out 
#SBATCH --error=slurm_prod_08.err 

# load appropriate modules.  First purge all modules then load in intel (a prereq for loading in amber) and then amber 
module purge
module load intel/2017.3.196
module load amber/18

# echo commands to stdout
set -x

# Executable
SANDER=pmemd.cuda

# Amber input files and output name
INP=08_prod.in
TOP=PDB_TEMP_solv.prmtop
CRD=07_prod.rst
OUT=08_prod

$SANDER  -O     -i   $INP   -p   $TOP   -c   $CRD   -r   $OUT.rst \
        -o   $OUT.out   -e   $OUT.ene   -v   $OUT.vel   -inf $OUT.nfo   -x   $OUT.nc &&

# make cpptraj input file
echo "parm PDB_TEMP_solv.prmtop"                            > strip_08.in
echo "trajin 08_prod.nc"                                    >> strip_08.in
echo "autoimage"                                            >> strip_08.in
echo "strip :WAT,Cl-,Na+"                                   >> strip_08.in
echo "trajout 08_prod_dry.nc"                               >> strip_08.in
echo "go"                                                   >> strip_08.in
echo "quit"                                                 >> strip_08.in

# run strip_08 to strip water
strip_08 -i strip_08.in > strip_08.out

# remove solvated nc file
if [ -f "08_prod_dry.nc" ] ; then
    rm 08_prod.nc &&
    mv 08_prod_dry.nc 08_prod.nc
fi

sbatch h2p_1gpu_prod_09.slurm
