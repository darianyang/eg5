#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/eg5_holo.prmtop .
ln -sv $WEST_SIM_ROOT/bstates/bstate.ncrst .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF/seg.ncrst ./parent.ncrst
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/md.in > md.in
  ln -sv $WEST_PARENT_DATA_REF ./parent.ncrst
fi

export CUDA_DEVICES=(`echo $CUDA_VISIBLE_DEVICES_ALLOCATED | tr , ' '`)
export CUDA_VISIBLE_DEVICES=${CUDA_DEVICES[$WM_PROCESS_INDEX]}

echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES_ALLOCATED = " $CUDA_VISIBLE_DEVICES_ALLOCATED
echo "RUNSEG.SH: WM_PROCESS_INDEX = " $WM_PROCESS_INDEX
echo "RUNSEG.SH: CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES

$PMEMD -O -i md.in   -p eg5_holo.prmtop  -c parent.ncrst \
          -r seg.ncrst -x seg.nc      -o seg.log    -inf seg.nfo

COMMAND="         parm eg5_holo.prmtop\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/parent.ncrst\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
COMMAND="$COMMAND reference bstate.ncrst\n"
COMMAND="$COMMAND autoimage \n"
COMMAND="$COMMAND rmsd :18-364&!@H= reference \n"
COMMAND="$COMMAND rmsd base-rmsd :370@C4*,O4*,C1*,N9,C8,N7,C5,C6,N6,N1,C2,N3,C4,C3*,O3*,C2*,O2* reference nofit out pcoord.dat \n"
COMMAND="$COMMAND rmsd phos-rmsd :370@O1B,PB,O2B,O3B,O3A,PA,O1A,O2A,O5*,C5* reference nofit out pcoord.dat \n"
COMMAND="$COMMAND nativecontacts mindist :1-368 :370 out pcoord.dat \n"
COMMAND="$COMMAND rmsd pro-rmsd :1-368&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd rec-rmsd :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd lig-rmsd :370&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd cpx-rmsd :1-370&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd loop5-rmsd :116-133&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND energy rec-ene :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338 out auxdata_ene.dat \n"
COMMAND="$COMMAND energy lig-ene :370 out auxdata_ene.dat \n"
COMMAND="$COMMAND energy cpx-ene :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338,370 out auxdata_ene.dat \n"
COMMAND="$COMMAND surf rec-sasa :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338 out auxdata_sasa.dat \n"
COMMAND="$COMMAND surf lig-sasa :370 out auxdata_sasa.dat \n"
COMMAND="$COMMAND go\n"

echo -e $COMMAND | $CPPTRAJ

python $WEST_SIM_ROOT/common_files/analyze_sasa.py $WEST_CURRENT_SEG_DATA_REF

paste <(cat pcoord.dat | tail -n +2 | awk '{print $2}') <(cat pcoord.dat | tail -n +2 | awk {'print $3'}) <(cat pcoord.dat | tail -n +2 | awk '{print $6}') > $WEST_PCOORD_RETURN

cat auxdata_rmsd.dat | tail -n +2 | awk {'print $2'} > $WEST_PROTEIN_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $3'} > $WEST_RECEPT_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $4'} > $WEST_LIGAND_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $5'} > $WEST_COMPLEX_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $6'} > $WEST_LOOP5_RMSD_RETURN
cat auxdata_ene.dat | tail -n +2 | awk {'print $10'} > $WEST_RECEPT_ENE_RETURN
cat auxdata_ene.dat | tail -n +2 | awk {'print $19'} > $WEST_LIGAND_ENE_RETURN
cat auxdata_ene.dat | tail -n +2 | awk {'print $28'} > $WEST_COMPLEX_ENE_RETURN
cat auxdata_sasa.dat | tail -n +2 | awk {'print $2'} > $WEST_RECEPT_SASA_CPT_RETURN
cat auxdata_sasa.dat | tail -n +2 | awk {'print $3'} > $WEST_LIGAND_SASA_CPT_RETURN

cat rec_sasa.dat > $WEST_RECEPT_SASA_MDT_RETURN
cat lig_sasa.dat > $WEST_LIGAND_SASA_MDT_RETURN

cp eg5_holo.prmtop $WEST_TRAJECTORY_RETURN
cp seg.nc $WEST_TRAJECTORY_RETURN

cp eg5_holo.prmtop $WEST_RESTART_RETURN
cp seg.ncrst $WEST_RESTART_RETURN/parent.ncrst

cp seg.log $WEST_LOG_RETURN

# Clean up
rm -f md.in parent.ncrst seg.nfo seg.pdb eg5_holo.prmtop bstate.ncrst *.dat
