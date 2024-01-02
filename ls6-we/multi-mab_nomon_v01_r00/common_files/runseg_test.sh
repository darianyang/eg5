#!/bin/bash


WEST_SIM_ROOT="/ix/lchong/dty7/eg5/h2p-we/multi-mab_nomon_v00"

COMMAND="         parm $WEST_SIM_ROOT/common_files/eg5_2022.prmtop\n"
COMMAND="$COMMAND trajin seg.ncrst\n"
COMMAND="$COMMAND reference $WEST_SIM_ROOT/bstates/05_eq3.ncrst\n"
COMMAND="$COMMAND autoimage \n"
COMMAND="$COMMAND rmsd :18-364&!@H= reference \n"
# pcoord calc
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
COMMAND="$COMMAND nativecontacts mindist :1-368 :370@O1B,PB,O2B,O3B,O3A,PA,O1A,O2A,O5*,C5* out pcoord.dat \n"

COMMAND="$COMMAND rmsd base-rmsd :370@C4*,O4*,C1*,N9,C8,N7,C5,C6,N6,N1,C2,N3,C4,C3*,O3*,C2*,O2* reference nofit out adp_rmsd.dat \n"
COMMAND="$COMMAND rmsd phos-rmsd :370@O1B,PB,O2B,O3B,O3A,PA,O1A,O2A,O5*,C5* reference nofit out adp_rmsd.dat \n"
# begin auxdata_rmsd.dat
COMMAND="$COMMAND rmsd prot-rmsd :1-368&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd rec-rmsd :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd adp-rmsd :370&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd adpcpx-rmsd :1-370&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd loop5-rmsd :116-133&!@H= reference nofit out auxdata_rmsd.dat \n"
COMMAND="$COMMAND rmsd protcut-rmsd :18-364&!@H= reference nofit out auxdata_rmsd.dat \n"
#COMMAND="$COMMAND rmsd mon-rmsd :371&!@H= reference nofit out auxdata_rmsd.dat \n"
#COMMAND="$COMMAND rmsd moncpx-rmsd :1-368,371&!@H= reference nofit out auxdata_rmsd.dat \n"
#COMMAND="$COMMAND rmsd adpmoncpx-rmsd :1-371&!@H= reference nofit out auxdata_rmsd.dat \n"
#COMMAND="$COMMAND energy rec-ene :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338 out auxdata_ene.dat \n"
#COMMAND="$COMMAND energy adp-ene :370 out auxdata_ene.dat \n"
#COMMAND="$COMMAND energy adpcpx-ene :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338,370 out auxdata_ene.dat \n"
#COMMAND="$COMMAND energy mon-ene :371 out auxdata_ene.dat \n"
# begin auxdata_sasa.dat
COMMAND="$COMMAND surf rec-sasa :24-30,72-78,81,103-118,132,232-233,265-270,274,333,335-338 out auxdata_sasa.dat \n"
COMMAND="$COMMAND surf adp-sasa :370 out auxdata_sasa.dat \n"
#COMMAND="$COMMAND surf mon-sasa :371 out auxdata_sasa.dat \n"
COMMAND="$COMMAND surf prot-sasa :1-368 out auxdata_sasa.dat \n"
COMMAND="$COMMAND go\n"

echo -e $COMMAND | cpptraj

python $WEST_SIM_ROOT/common_files/analyze_sasa.py .

# calc interaction energy
bash $WEST_SIM_ROOT/common_files/get_energies.sh pcoord-ene.dat int_ene.txt
python $WEST_SIM_ROOT/common_files/get_energies.py int_ene.txt int_ene_proc.txt

paste <(cat pcoord.dat | tail -n +2 | awk '{print $2}') <(cat int_ene_proc.txt) <(cat pcoord.dat | tail -n +2 | awk '{print $5}') <(cat pcoord.dat | tail -n +2 | awk '{print $8}')  > pcoord.out
#paste <(cat pcoord.dat | tail -n +2 | awk '{print $2}') <(cat pcoord.dat | tail -n +2 | awk '{print $7}') > $WEST_PCOORD_RETURN

