#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/1x88_solv.prmtop .
ln -sv $WEST_SIM_ROOT/bstates/05_eq3.ncrst .

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

$PMEMD -O -i md.in   -p 1x88_solv.prmtop  -c parent.ncrst \
          -r seg.ncrst -x seg.nc      -o seg.log    -inf seg.nfo

COMMAND="         parm 1x88_solv.prmtop\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/parent.ncrst\n"
COMMAND="$COMMAND trajin $WEST_CURRENT_SEG_DATA_REF/seg.nc\n"
COMMAND="$COMMAND reference 05_eq3.ncrst\n"
COMMAND="$COMMAND autoimage \n"
COMMAND="$COMMAND rmsd :18-364&!@H= reference \n"
# pcoord calc
COMMAND="${COMMAND} rmsd adp-mg-rmsd :369-370&!@H= reference nofit out pcoord.dat \n"
COMMAND="$COMMAND rmsd base-rmsd :370@C4*,O4*,C1*,N9,C8,N7,C5,C6,N6,N1,C2,N3,C4,C3*,O3*,C2*,O2* reference nofit out pcoord.dat \n"
COMMAND="$COMMAND rmsd phos-rmsd :370@O1B,PB,O2B,O3B,O3A,PA,O1A,O2A,O5*,C5* reference nofit out pcoord.dat \n"
COMMAND="$COMMAND nativecontacts mindist :1-368 :370 out pcoord.dat \n"
# begin auxdata_rmsd.dat
COMMAND="$COMMAND rmsd prot-rmsd :1-368&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd rec-rmsd :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd adp-rmsd :370&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd adpcpx-rmsd :1-370&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd loop5-rmsd :116-133&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd protcut-rmsd :18-364&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd mon-rmsd :371&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd moncpx-rmsd :1-368,371&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd adpmoncpx-rmsd :1-371&!@H= reference nofit out auxdata_rmsd.dat \n"
# begin auxdata_ene.dat
COMMAND="$COMMAND energy rec-ene :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338 out auxdata_ene.dat \n"
COMMAND="$COMMAND energy adp-ene :370 out auxdata_ene.dat \n"
COMMAND="$COMMAND energy adpcpx-ene :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338,370 out auxdata_ene.dat \n"
COMMAND="$COMMAND energy mon-ene :371 out auxdata_ene.dat \n"
# begin auxdata_sasa.dat
COMMAND="$COMMAND surf rec-sasa :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338 out auxdata_sasa.dat \n"
COMMAND="$COMMAND surf adp-sasa :370 out auxdata_sasa.dat \n"
COMMAND="$COMMAND surf mon-sasa :371 out auxdata_sasa.dat \n"
COMMAND="$COMMAND surf prot-sasa :1-368 out auxdata_sasa.dat \n"
# pcoord ene calc: needs to be in order of: complex, receptor, ligand
COMMAND="$COMMAND energy adpcpx-lt-8A-ene :370<:8 out pcoord-ene.dat \n"
COMMAND="$COMMAND energy rec-lt-8A-ene :370<:8&!(:370) out pcoord-ene.dat \n"
COMMAND="$COMMAND energy adp-ene-pcoord :370 out pcoord-ene.dat \n"
COMMAND="$COMMAND go\n"

echo -e $COMMAND | $CPPTRAJ

python $WEST_SIM_ROOT/common_files/analyze_sasa.py $WEST_CURRENT_SEG_DATA_REF
# calc interaction energy
bash $WEST_SIM_ROOT/common_files/get_energies.sh pcoord-ene.dat int_ene.txt
python $WEST_SIM_ROOT/common_files/get_energies.py int_ene.txt int_ene_proc.txt

paste <(cat pcoord.dat | tail -n +2 | awk '{print $2}') <(cat int_ene_proc.txt) <(cat pcoord.dat | tail -n +2 | awk '{print $7}') > $WEST_PCOORD_RETURN

cat pcoord.dat | tail -n +2 | awk {'print $3'} > $WEST_ADP_BASE_RMSD_RETURN
cat pcoord.dat | tail -n +2 | awk {'print $4'} > $WEST_ADP_PHOS_RMSD_RETURN

cat auxdata_rmsd.dat | tail -n +2 | awk {'print $2'} > $WEST_PROT_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $3'} > $WEST_RECEPT_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $4'} > $WEST_ADP_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $5'} > $WEST_ADPCOMPLEX_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $6'} > $WEST_LOOP5_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $7'} > $WEST_PROTCUT_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $8'} > $WEST_MON_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $9'} > $WEST_MONCOMPLEX_RMSD_RETURN
cat auxdata_rmsd.dat | tail -n +2 | awk {'print $10'} > $WEST_ADPMONCOMPLEX_RMSD_RETURN

cat auxdata_ene.dat | tail -n +2 | awk {'print $10'} > $WEST_RECEPT_ENE_RETURN
cat auxdata_ene.dat | tail -n +2 | awk {'print $19'} > $WEST_ADP_ENE_RETURN
cat auxdata_ene.dat | tail -n +2 | awk {'print $28'} > $WEST_ADPCOMPLEX_ENE_RETURN
cat auxdata_ene.dat | tail -n +2 | awk {'print $37'} > $WEST_MON_ENE_RETURN

cat auxdata_sasa.dat | tail -n +2 | awk {'print $2'} > $WEST_RECEPT_SASA_CPT_RETURN
cat auxdata_sasa.dat | tail -n +2 | awk {'print $3'} > $WEST_ADP_SASA_CPT_RETURN
cat auxdata_sasa.dat | tail -n +2 | awk {'print $4'} > $WEST_MON_SASA_CPT_RETURN
cat auxdata_sasa.dat | tail -n +2 | awk {'print $5'} > $WEST_PROT_SASA_CPT_RETURN

cat prot_sasa.dat > $WEST_PROT_SASA_MDT_RETURN
cat recept_sasa.dat > $WEST_RECEPT_SASA_MDT_RETURN
cat adp_sasa.dat > $WEST_ADP_SASA_MDT_RETURN
cat mon_sasa.dat > $WEST_MON_SASA_MDT_RETURN

cp 1x88_solv.prmtop $WEST_TRAJECTORY_RETURN
cp seg.nc $WEST_TRAJECTORY_RETURN

cp 1x88_solv.prmtop $WEST_RESTART_RETURN
cp seg.ncrst $WEST_RESTART_RETURN/parent.ncrst

cp seg.log $WEST_LOG_RETURN

# strip water
CMD="     parm $WEST_SIM_ROOT/common_files/1x88_solv.prmtop \n"
CMD="$CMD trajin $WEST_CURRENT_SEG_DATA_REF/seg.nc \n"
CMD="$CMD autoimage \n"
# strip and replace solv nc file, make seperate stripped rst file but keep solved
CMD="$CMD strip :WAT parmout $WEST_CURRENT_SEG_DATA_REF/seg-nowat.rst \n"
CMD="$CMD trajout $WEST_CURRENT_SEG_DATA_REF/seg-nowat.nc \n"
CMD="$CMD go \n"

echo -e "$CMD" | $CPPTRAJ

# remove and replace solvated segment trajectory file
if [ -f "seg-nowat.nc" ]; then
    rm seg.nc &&
    mv seg-nowat.nc seg.nc
fi

# Clean up
rm -f md.in parent.ncrst seg.nfo seg.pdb 1x88_solv.prmtop 05_eq3.ncrst *.dat *.txt
