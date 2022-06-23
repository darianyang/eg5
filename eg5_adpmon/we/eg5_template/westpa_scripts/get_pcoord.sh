#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT

COMMAND="           parm $WEST_SIM_ROOT/common_files/eg5_holo.prmtop \n"
COMMAND="${COMMAND} trajin $WEST_STRUCT_DATA_REF \n"
COMMAND="${COMMAND} reference $WEST_SIM_ROOT/bstates/bstate.ncrst \n"
COMMAND="${COMMAND} autoimage \n"
COMMAND="${COMMAND} rmsd :18-364&!@H= reference \n"
COMMAND="${COMMAND} rmsd base-rmsd :370@C4*,O4*,C1*,N9,C8,N7,C5,C6,N6,N1,C2,N3,C4,C3*,O3*,C2*,O2* reference nofit out pcoord.dat \n"
COMMAND="${COMMAND} rmsd phos-rmsd :370@O1B,PB,O2B,O3B,O3A,PA,O1A,O2A,O5*,C5* reference nofit out pcoord.dat \n"
COMMAND="${COMMAND} nativecontacts mindist :1-368 :370 out pcoord.dat \n"
COMMAND="${COMMAND} rmsd pro-rmsd :1-368&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="${COMMAND} rmsd rec-rmsd :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="${COMMAND} rmsd lig-rmsd :370&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="${COMMAND} rmsd cpx-rmsd :1-370&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="${COMMAND} rmsd loop5-rmsd :116-133&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="${COMMAND} energy rec-ene :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338 out auxdata_ene.dat \n"
COMMAND="${COMMAND} energy lig-ene :370 out auxdata_ene.dat \n"
COMMAND="${COMMAND} energy cpx-ene :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338,370 out auxdata_ene.dat \n"
COMMAND="${COMMAND} surf rec-sasa :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338 out auxdata_sasa.dat \n"
COMMAND="${COMMAND} surf lig-sasa :370 out auxdata_sasa.dat \n"
COMMAND="${COMMAND} go"

echo -e "${COMMAND}" | $CPPTRAJ


paste <(cat pcoord.dat | tail -n +2 | awk '{print $2}') <(cat pcoord.dat | tail -n +2 | awk {'print $3'}) <(cat pcoord.dat | tail -n +2 | awk '{print $6}') > $WEST_PCOORD_RETURN


#rm pcoord.dat auxdata_rmsd.dat auxdata_ene.dat auxdata_sasa.dat rec_sasa.dat lig_sasa.dat

cp $WEST_SIM_ROOT/common_files/eg5_holo.prmtop $WEST_TRAJECTORY_RETURN
cp $WEST_STRUCT_DATA_REF $WEST_TRAJECTORY_RETURN

cp $WEST_SIM_ROOT/common_files/eg5_holo.prmtop $WEST_RESTART_RETURN
cp $WEST_STRUCT_DATA_REF $WEST_RESTART_RETURN/parent.ncrst

#python $WEST_SIM_ROOT/common_files/analyze_sasa.py
#cat rec_sasa.dat
#cat lig_sasa.dat

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
