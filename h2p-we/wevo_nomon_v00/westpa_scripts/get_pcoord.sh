#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT

COMMAND="           parm $WEST_SIM_ROOT/common_files/eg5_2022.prmtop \n"
COMMAND="${COMMAND} trajin $WEST_STRUCT_DATA_REF \n"
COMMAND="${COMMAND} reference $WEST_SIM_ROOT/bstates/05_eq3.ncrst \n"
COMMAND="${COMMAND} autoimage \n"
COMMAND="${COMMAND} rmsd :18-364&!@H= reference \n"
COMMAND="${COMMAND} rmsd adp-mg-rmsd :368-369&!@H= reference nofit out pcoord.dat \n"
COMMAND="${COMMAND} rmsd base-rmsd :369@C4*,O4*,C1*,N9,C8,N7,C5,C6,N6,N1,C2,N3,C4,C3*,O3*,C2*,O2* reference nofit out pcoord.dat \n"
COMMAND="${COMMAND} rmsd phos-rmsd :369@O1B,PB,O2B,O3B,O3A,PA,O1A,O2A,O5*,C5* reference nofit out pcoord.dat \n"
COMMAND="${COMMAND} nativecontacts mindist :1-367 :369 out pcoord.dat \n"
COMMAND="${COMMAND} go"

echo -e "${COMMAND}" | $CPPTRAJ

paste <(cat pcoord.dat | tail -n +2 | awk '{print $2}') <(cat pcoord.dat | tail -n +2 | awk '{print $7}') > $WEST_PCOORD_RETURN

#rm pcoord.dat auxdata_rmsd.dat auxdata_ene.dat auxdata_sasa.dat
#rm pcoord.dat

cp $WEST_SIM_ROOT/common_files/eg5_2022.prmtop $WEST_TRAJECTORY_RETURN
cp $WEST_STRUCT_DATA_REF $WEST_TRAJECTORY_RETURN

cp $WEST_SIM_ROOT/common_files/eg5_2022.prmtop $WEST_RESTART_RETURN
cp $WEST_STRUCT_DATA_REF $WEST_RESTART_RETURN/parent.ncrst

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
