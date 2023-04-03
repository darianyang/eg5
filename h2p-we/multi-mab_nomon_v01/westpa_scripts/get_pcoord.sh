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

# pcoord 1: ADP+Mg RMSD
COMMAND="${COMMAND} rmsd adp-mg-rmsd :369-370&!@H= reference nofit out pcoord.dat \n"
#COMMAND="${COMMAND} rmsd base-rmsd :370@C4*,O4*,C1*,N9,C8,N7,C5,C6,N6,N1,C2,N3,C4,C3*,O3*,C2*,O2* reference nofit out pcoord.dat \n"
#COMMAND="${COMMAND} rmsd phos-rmsd :370@O1B,PB,O2B,O3B,O3A,PA,O1A,O2A,O5*,C5* reference nofit out pcoord.dat \n"

# pcoord 2: receptor and ADP + Mg interaction energy
# needs to be in order of: complex, receptor, ligand
COMMAND="$COMMAND energy cpx-ene :1-370 out pcoord-ene.dat \n"
COMMAND="$COMMAND energy rec-ene :1-368 out pcoord-ene.dat \n"
COMMAND="$COMMAND energy adp-ene :369,370 out pcoord-ene.dat \n"

# pcoord 3: min contact distance of ADP and Eg5
COMMAND="${COMMAND} nativecontacts mindist :1-368 :370 out pcoord.dat \n"

# pcoord 4: min contact distance of PO4 and eg5
COMMAND="${COMMAND} nativecontacts mindist :1-368 :370@O1B,PB,O2B,O3B,O3A,PA,O1A,O2A,O5*,C5* out pcoord.dat \n"

# done
COMMAND="$COMMAND go\n"
#echo -e "${COMMAND}" | $CPPTRAJ
echo -e "${COMMAND}" | cpptraj

# calc interaction energy
bash $WEST_SIM_ROOT/common_files/get_energies.sh pcoord-ene.dat int_ene.txt
python $WEST_SIM_ROOT/common_files/get_energies.py int_ene.txt int_ene_proc.txt

paste <(cat pcoord.dat | tail -n +2 | awk '{print $2}') <(cat int_ene_proc.txt) <(cat pcoord.dat | tail -n +2 | awk '{print $5}') <(cat pcoord.dat | tail -n +2 | awk '{print $8}')  > $WEST_PCOORD_RETURN
paste <(cat pcoord.dat | tail -n +2 | awk '{print $2}') <(cat int_ene_proc.txt) <(cat pcoord.dat | tail -n +2 | awk '{print $5}') <(cat pcoord.dat | tail -n +2 | awk '{print $8}')  > pcoord_init.dat

rm pcoord.dat pcoord-ene.dat int_ene.txt int_ene_proc.txt
#rm pcoord.dat

#cp $WEST_SIM_ROOT/common_files/eg5_2022.prmtop $WEST_TRAJECTORY_RETURN
#cp $WEST_STRUCT_DATA_REF $WEST_TRAJECTORY_RETURN

#cp $WEST_SIM_ROOT/common_files/eg5_2022.prmtop $WEST_RESTART_RETURN
#cp $WEST_STRUCT_DATA_REF $WEST_RESTART_RETURN/parent.ncrst

if [ -n "$SEG_DEBUG" ] ; then
  head -v $WEST_PCOORD_RETURN
fi
