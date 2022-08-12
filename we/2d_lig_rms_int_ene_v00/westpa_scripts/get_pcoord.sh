#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT

COMMAND="           parm $WEST_SIM_ROOT/common_files/1x88_solv.prmtop \n"
COMMAND="$COMMAND trajin $WEST_STRUCT_DATA_REF \n"
COMMAND="$COMMAND reference $WEST_SIM_ROOT/bstates/05_eq3.ncrst \n"
COMMAND="$COMMAND autoimage \n"
COMMAND="$COMMAND rmsd :18-364&!@H= reference \n"
COMMAND="$COMMAND rmsd adp-mg-rmsd :369-370&!@H= reference nofit mass out pcoord.dat \n"
COMMAND="$COMMAND distance adp-mg-proteincut-dist :369-370&!@H= :18-364&!@H= out pcoord.dat \n"
COMMAND="$COMMAND rmsd base-rmsd :370@C4*,O4*,C1*,N9,C8,N7,C5,C6,N6,N1,C2,N3,C4,C3*,O3*,C2*,O2* reference nofit out pcoord.dat \n"
COMMAND="$COMMAND rmsd phos-rmsd :370@O1B,PB,O2B,O3B,O3A,PA,O1A,O2A,O5*,C5* reference nofit out pcoord.dat \n"
COMMAND="$COMMAND nativecontacts mindist :1-368 :370 out pcoord.dat \n"

# testing distance selection syntax
#COMMAND="$COMMAND nativecontacts byresidue :370<:8&!(:370) :369 resout test.dat \n"

COMMAND="$COMMAND energy adpcpx-lt-8A-ene :370<:8 out pcoord_ene.dat \n"
COMMAND="$COMMAND energy rec-lt-8A-ene :370<:8&!(:370) out pcoord_ene.dat \n"
COMMAND="$COMMAND energy adp-ene :370 out pcoord_ene.dat \n"
COMMAND="$COMMAND go"

echo -e "$COMMAND" | $CPPTRAJ

# calc interaction energy
bash $WEST_SIM_ROOT/common_files/get_energies.sh pcoord_ene.dat ene.txt
python $WEST_SIM_ROOT/common_files/get_energies.py ene.txt ene_proc.txt

paste <(cat pcoord.dat | tail -n +2 | awk '{print $2}') <(cat ene_proc.txt) <(cat pcoord.dat | tail -n +2 | awk '{print $8}') > $WEST_PCOORD_RETURN

#rm pcoord.dat auxdata_rmsd.dat auxdata_ene.dat auxdata_sasa.dat
#rm pcoord.dat

cp $WEST_SIM_ROOT/common_files/1x88_solv.prmtop $WEST_TRAJECTORY_RETURN
cp $WEST_STRUCT_DATA_REF $WEST_TRAJECTORY_RETURN

cp $WEST_SIM_ROOT/common_files/1x88_solv.prmtop $WEST_RESTART_RETURN
cp $WEST_STRUCT_DATA_REF $WEST_RESTART_RETURN/parent.ncrst

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
