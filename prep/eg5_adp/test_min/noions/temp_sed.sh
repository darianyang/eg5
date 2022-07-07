#!/bin/bash
# temp_sed.sh
# replace 1ii6 and v00 in the std sim scripts dir with target PDB

# Parameters
# ----------
# arg 1 = pdb system prefix         :   e.g. 2kod_wt
# arg 2 = version of the replicate  :   e.g. v00

PDB=$1
v00=$2

# apply globally to all files in current directory
sed -i "s/1ii6/${PDB}/g" *
sed -i "s/v00/${VER}/" * 
