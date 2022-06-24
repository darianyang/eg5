#!/bin/bash
# analyze multiple (5) replicates of each system

SYSTEMS=(w4f w5f w6f w7f wt)
#SYSTEMS=(w4f)
#SYSTEMS=(wt)
# for iRED: no need for wt
#SYSTEMS=(w5f w6f w7f)
OUT_ROOT=1us_noion
PDB=3k0n
FF="ipq"
IRED_V=1816 # starting at W4F CE3
#IRED_V=1814 # W5F
CPPTRAJ=cpptraj.MPI
export DO_PARALLEL="mpirun -np 8"

cd $FF

# go to target directory
for SYS in ${SYSTEMS[@]} ; do
    cd $SYS &&
    echo "RUNNING SYSTEM : $SYS"

    # make the 5 replicate sub-directories
    for V in {00..04..1} ; do
        echo "GENERATING REPLICA V$V"
        # go to sub-directory
        cd v$V
        mkdir $OUT_ROOT

#        # Commands 00: Strip water and/or ions
#        C00="     parm ${PDB}_${SYS}_solv.prmtop \n"
#        C00="$C00 trajin 06_prod.nc \n"
#        C00="$C00 strip :WAT \n"
##        C00="$C00 trajout ${PDB}_${SYS}_w_ion.pdb pdb \n"
#        C00="$C00 trajout 06_prod_w_ion.nc \n"
#        echo -e "$C00" | $CPPTRAJ \
#         > >(tee $OUT_ROOT/cpptraj_0.out) \
#        2> >(tee $OUT_ROOT/cpptraj_0.err >&2)

        # Commands 0: Load in topology and trajectory
        #if [ -f "${PDB}_${SYS}_dry_noion.prmtop" ] ; then
        #    C0="    parm ${PDB}_${SYS}_dry_noion.prmtop \n"
        #else
        #    C0="    parm ${PDB}_${SYS}_dry.prmtop \n"
        #fi 
        C0="    parm ../${PDB}_${SYS}_dry_noion.prmtop \n"
        #C0="    parm ../${PDB}_${SYS}_dry.prmtop \n"
        # no water but retains ions
        C0="$C0 trajin 06_prod_dry.nc 1 last 10 \n" 
        C0="$C0 trajin 07_prod_dry.nc 1 last 10 \n" 
        C0="$C0 trajin 08_prod_dry.nc 1 last 10 \n" 
        C0="$C0 reference ${PDB}_${SYS}_w_ion.pdb :* [REF] \n"
        
        # Commands 1: Dihedrals, radius of gyration, RMSD, DSSP
        C1="    multidihedral dihedral_trp121 phi psi resrange 121-121"
        C1="$C1               out $OUT_ROOT/dihedral_trp121.dat \n"
        C1="$C1 rms heavyRMSD :1-165&!@H= out $OUT_ROOT/rmsd_1-165_heavy.dat ref [REF] mass time 1 \n"
        C1="$C1 rms bbRMSD :1-165@CA,C,O,N out $OUT_ROOT/rmsd_1-165_bb.dat ref [REF] mass time 1 \n"

        if [ $SYS == "wt" ] ; then
            # COM for 4 TRP WT H positions : prob can use a loop of H atom ids
            C1="$C1 distance 4F_to_BB_CO @1819,1820 @1817 out $OUT_ROOT/dist_4F_to_BB_CO.dat \n"
            C1="$C1 distance 5F_to_BB_CO @1819,1820 @1815 out $OUT_ROOT/dist_5F_to_BB_CO.dat \n"
            C1="$C1 distance 6F_to_BB_CO @1819,1820 @1813 out $OUT_ROOT/dist_6F_to_BB_CO.dat \n"
            C1="$C1 distance 7F_to_BB_CO @1819,1820 @1811 out $OUT_ROOT/dist_7F_to_BB_CO.dat \n"
        else
            # COM between C=O and 19F
            C1="$C1 distance F_to_BB_CO @1819,1820 @$(($IRED_V + 1)) out $OUT_ROOT/dist_F_to_BB_CO.dat \n"
        fi

        #C1="$C1 radgyr :1-165"
        #C1="$C1        out $OUT_ROOT/radgyr.dat \n"
        #C1="$C1 rms rms_3K0N_1-165_BB :1-165@N,CA,C ref [REF] mass"
        #C1="$C1     perres    perresmask @N,CA,C"
        #C1="$C1     perresout $OUT_ROOT/rmsd_3K0N_1-165_BB_perres.dat"
        #C1="$C1     perresavg $OUT_ROOT/rmsd_3K0N_1-165_BB_perres_avg.dat"
        #C1="$C1     out       $OUT_ROOT/rmsd_3K0N_1-165_BB.dat \n"
        #C1="$C1 secstruct ss :*"
        #C1="$C1           out       $OUT_ROOT/ss.dat"
        #C1="$C1           sumout    $OUT_ROOT/ss_sum.dat"
        #C1="$C1           assignout $OUT_ROOT/ss_assign.dat \n"

        echo -e "${C0}${C1}run" > $OUT_ROOT/cpp_c0c1.in
        $DO_PARALLEL $CPPTRAJ -i $OUT_ROOT/cpp_c0c1.in > cpp_c0c1.out
        
#        echo -e "$C0$C1" | $CPPTRAJ \
#         > >(tee $OUT_ROOT/cpptraj_1.out) \
#        2> >(tee $OUT_ROOT/cpptraj_1.err >&2)

#        # Commands 4: iRED
#        C4=""
#        C4="$C4 vector CF @${IRED_V} ired @$(($IRED_V + 1)) \n"
#        C4="$C4 matrix ired name ired_matrix order 2 \n"
#        C4="$C4 analyze matrix ired_matrix name ired_vectors"
#        C4="$C4         vecs 1"
#        C4="$C4         out $OUT_ROOT/ired_vectors.dat \n"
#        
#        C4="$C4 ired relax freq 600 NHdist 1.41 order 2 norm"
#        C4="$C4      tstep 1 tcorr 41000" # 5 * 8.2ns tc of 19F CypA
#        C4="$C4      modes ired_vectors"
#        C4="$C4      out            $OUT_ROOT/ired_${CUTOFF}_corr.dat"
#        C4="$C4      noefile        $OUT_ROOT/ired_${CUTOFF}_relax.dat"
#        C4="$C4      orderparamfile $OUT_ROOT/ired_${CUTOFF}_order.dat \n"
#        
#        echo -e "$C0$C4" | $CPPTRAJ \
#         > >(tee $OUT_ROOT/cpptraj_4.out) \
#        2> >(tee $OUT_ROOT/cpptraj_4.err >&2)

        
        cd ..
    done
    let "IRED_V-=2"
    cd ..
    echo "FINISHED RUNNING SYSTEM : $SYS"
done
