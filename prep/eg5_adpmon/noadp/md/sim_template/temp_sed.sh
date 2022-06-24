#!/bin/bash
# temp_sed.sh
# replace PDB_TEMP and VER in the std sim scripts dir with target PDB

# Parameters
# ----------
# arg 1 = pdb system prefix         :   e.g. 2kod_wt
# arg 2 = version of the replicate  :   e.g. v00

PDB=$1
VER=$2

# apply globally to all files in current directory
sed -i "s/PDB_TEMP/${PDB}/g" *
sed -i "s/VER/${VER}/" * 
