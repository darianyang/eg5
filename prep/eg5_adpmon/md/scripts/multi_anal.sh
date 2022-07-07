#!/bin/bash
# multi_anal.sh
# analyze multiple (5) replicates of each system

# for iRED: no need for wt
SYSTEMS=(w4f w5f)
#SYSTEMS=(wt)
OUT=200ns-nocl
PDB=3k0n
FF="ipq"
#IRED_V=1816 # starting at W4F CE3
IRED_V=1812 # starting at W6F

# go to target directory
for SYS in ${SYSTEMS[@]} ; do
    cd $FF/$SYS &&
    echo "RUNNING SYSTEM : $SYS"

    # make the 5 replicate sub-directories
    for V in {00..04..1} ; do
        echo "GENERATING REPLICA V$V"
        # go to sub-directory
        cd v$V
        mkdir $OUT

#        echo "parm 3k0n_${SYS}_solv.prmtop" > ${OUT}/strip.cpp
#        echo "trajin 06_prod.nc" >> ${OUT}/strip.cpp
#        echo "autoimage" >> ${OUT}/strip.cpp
#        echo "strip :WAT,Cl- parmout 3k0n_${SYS}_dry_noion.prmtop" >> ${OUT}/strip.cpp
#        echo "trajout 06_prod_dry.nc" >> ${OUT}/strip.cpp
#        echo "go" >> ${OUT}/strip.cpp
#        echo "quit" >> ${OUT}/strip.cpp
#        cpptraj -i ${OUT}/strip.cpp > ${OUT}/strip.out

        # later rm 06_prod.nc and replace with 06_prod_dry.nc

        REF=${PDB}_${SYS}_solv.pdb

        # calc RMSD
        echo "parm ${PDB}_${SYS}_dry_noion.prmtop" > ${OUT}/calc.cpp
        echo "trajin 06_prod_dry.nc" >> ${OUT}/calc.cpp
        echo "reference ${REF} name <st0>" >> ${OUT}/calc.cpp
        echo "rms heavyRMSD :*&!@H= out ${OUT}/${PDB}_rms_heavy.dat ref <st0> mass time 1" >> ${OUT}/calc.cpp
        echo "rms trpRMSD :121&!@H= out ${OUT}/${PDB}_rms_trp.dat ref <st0> mass time 1" >> ${OUT}/calc.cpp

#        # ired analysis: define CF vector, diagonalize matrix, run ired using C-F bond distance
#        echo "vector CF @${IRED_V} ired @$(($IRED_V + 1))" >> ${OUT}/calc.cpp
#        echo "matrix ired name matired order 2" >> ${OUT}/calc.cpp
#        echo "diagmatrix matired vecs 1 out ${OUT}/ired.vec name ired.vec" >> ${OUT}/calc.cpp
#        # tstep is the crd saving steps, tcorr is 5x tc of CypA (5 * 8.2ns = 41ns)
#        echo "ired relax NHdist 1.41 freq 600.0 tstep 1 tcorr 41000 out ${OUT}/CF.out noefile ${OUT}/noe order 2 modes ired.vec" >> ${OUT}/calc.cpp

        echo "run" >> ${OUT}/calc.cpp
        echo "quit" >> ${OUT}/calc.cpp

        # run cpptraj
        cpptraj -i ${OUT}/calc.cpp > ${OUT}/calc.out

        cd ..
    done

    let "IRED_V-=2"
    echo "FINISHED SYSTEM : $SYS"
    cd ../..
done        

