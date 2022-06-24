#!/usr/bin/env bash
############################ MODULES AND SETTINGS #############################
source /opt/modules/sgefix.sh
set -e
module purge
module load anaconda
module load ktd3/amber/cpptraj
#set -x
############################ ENVIRONMENT VARIABLES ############################
PROJ_ROOT=$WORK/AMBER15IPQ
CLASS=GB3
CLASS_ROOT=$PROJ_ROOT/protein/$CLASS

FORCE_FIELDS=(ff15ipq10.3_spceb)
REF=$CLASS_ROOT/ff15ipq10.3_spceb/1_leap/${CLASS}_solute.rst
INDEXES=($(seq 2 1 56))
FREQUENCIES=(400 500 600 700 800 900)
CUTOFFS=(1000 2000 5000 10000 17000 20000 50000)

CPPTRAJ="$(which cpptraj)"
H5="$(which python) /share/home/ktd3/scripts/moldynplot/cpptraj2hdf5.py"
#################################### MAIN #####################################
for FORCE_FIELD in ${FORCE_FIELDS[@]}; do
    SIM_ROOT=$CLASS_ROOT/$FORCE_FIELD
    TOP=$SIM_ROOT/1_leap/${CLASS}_solute.parm7
    for NVE_RUN in $(ls $SIM_ROOT | grep NVE); do
        RUN_ROOT=$SIM_ROOT/$NVE_RUN
        echo $RUN_ROOT
        OUT_ROOT=$RUN_ROOT/analysis

        [ -d $OUT_ROOT ] ||
            mkdir -v $OUT_ROOT

        # Commands 0: Load in topology and trajectory
        C0="    parm ${TOP} \n"
        for SEG in $(ls $RUN_ROOT/[0-9][0-9][0-9][0-9]/*_solute.nc ); do
            C0="$C0 trajin ${SEG} 1 last \n"
        done
        C0="$C0 reference ${REF} :* [1P7E] \n"

#        # Commands 1: Dihedrals, radius of gyration, RMSD, DSSP
#        C1="    multidihedral dihedral phi psi"
#        C1="$C1               out $OUT_ROOT/dihedral.dat \n"
#        C1="$C1 radgyr !:CL,NA"
#        C1="$C1        out $OUT_ROOT/radgyr.dat \n"
#        C1="$C1 rms rms_1P7E_3-54_BB :3-54@N,CA,C ref [1P7E] mass"
#        C1="$C1     perres    perresmask @N,CA,C"
#        C1="$C1     perresout $OUT_ROOT/rmsd_1P7E_3-54_BB_perres.dat"
#        C1="$C1     perresavg $OUT_ROOT/rmsd_1P7E_3-54_BB_perres_avg.dat"
#        C1="$C1     out       $OUT_ROOT/rmsd_1P7E_3-54_BB.dat \n"
#        C1="$C1 secstruct ss :*"
#        C1="$C1           out       $OUT_ROOT/ss.dat"
#        C1="$C1           sumout    $OUT_ROOT/ss_sum.dat"
#        C1="$C1           assignout $OUT_ROOT/ss_assign.dat \n"
#
#        echo -e "$C0$C1" | $CPPTRAJ \
#         > >(tee $OUT_ROOT/cpptraj_1.out) \
#        2> >(tee $OUT_ROOT/cpptraj_1.err >&2)
#        $H5 dihedral   $OUT_ROOT/dihedral.dat                 $OUT_ROOT/dihedral.h5
#        rm -v $OUT_ROOT/dihedral.dat
#        $H5 perresrmsd $OUT_ROOT/rmsd_1P7E_3-54_BB_perres.dat $OUT_ROOT/rmsd_1P7E_3-54_BB_perres.h5
#        rm -v $OUT_ROOT/rmsd_1P7E_3-54_BB_perres.dat
#        $H5 secstruct  $OUT_ROOT/ss.dat                       $OUT_ROOT/ss.h5
#        rm -v $OUT_ROOT/ss.dat

#        # Commands 2: Rotational diffusion
#        C2="    rms rms_1P7E_3-54_BB :3-54@N,CA,C ref [1P7E] mass"
#        C2="$C2     savematrices \n"
#        C2="$C2 rotdif"
#        C2="$C2        usefft"
#        C2="$C2        rmatrix rms_1P7E_3-54_BB[RM] nvecs 1000 rseed -1"
#        C2="$C2        order 2 ncorr 1000"
#        C2="$C2        tf 1 dt 0.001"
#        C2="$C2        outfile $OUT_ROOT/rotdif_fft_ncorr-01000.dat \n"
#        C2="$C2 rotdif"
#        C2="$C2        rmatrix rms_1P7E_3-54_BB[RM] nvecs 1000 rseed -1"
#        C2="$C2        order 2 ncorr 1000"
#        C2="$C2        tf 1 dt 0.001"
#        C2="$C2        outfile $OUT_ROOT/rotdif_ncorr-01000_tf-00001.dat \n"
#
#        echo -e "$C0$C2" | $CPPTRAJ \
#         > >(tee $OUT_ROOT/cpptraj_2.out) \
#        2> >(tee $OUT_ROOT/cpptraj_2.err >&2)

#        # Commands 3: N-H vector autocorrelation functions
#        C3=""
#        for INDEX in ${INDEXES[@]}; do
#            C3="$C3 vector   v$INDEX :$INDEX@N :$INDEX@H\n"
#            C3="$C3 timecorr vec1 v$INDEX"
#            C3="$C3          tstep 1 tcorr 1000 order 2"
#            C3="$C3          name ac01000_$INDEX"
#            C3="$C3          norm"
#            C3="$C3          out $OUT_ROOT/relax_ac_01000.dat \n"
#            C3="$C3 timecorr vec1 v$INDEX"
#            C3="$C3          tstep 1 tcorr 2000 order 2"
#            C3="$C3          name ac02000_$INDEX"
#            C3="$C3          norm"
#            C3="$C3          out $OUT_ROOT/relax_ac_02000.dat \n"
#            C3="$C3 timecorr vec1 v$INDEX"
#            C3="$C3          tstep 1 tcorr 5000 order 2"
#            C3="$C3          name ac05000_$INDEX"
#            C3="$C3          norm"
#            C3="$C3          out $OUT_ROOT/relax_ac_05000.dat \n"
#            C3="$C3 timecorr vec1 v$INDEX"
#            C3="$C3          tstep 1 tcorr 10000 order 2"
#            C3="$C3          name ac10000_$INDEX"
#            C3="$C3          norm"
#            C3="$C3          out $OUT_ROOT/relax_ac_10000.dat \n"
#            C3="$C3 timecorr vec1 v$INDEX"
#            C3="$C3          tstep 1 tcorr 20000 order 2"
#            C3="$C3          name ac20000_$INDEX"
#            C3="$C3          norm"
#            C3="$C3          out $OUT_ROOT/relax_ac_20000.dat \n"
#            C3="$C3 timecorr vec1 v$INDEX"
#            C3="$C3          tstep 1 tcorr 50000 order 2"
#            C3="$C3          name ac50000_$INDEX"
#            C3="$C3          norm"
#            C3="$C3          out $OUT_ROOT/relax_ac_50000.dat \n"
#            C3="$C3 timecorr vec1 v$INDEX"
#            C3="$C3          tstep 1 tcorr 17000 order 2"
#            C3="$C3          name ac17000_$INDEX"
#            C3="$C3          norm"
#            C3="$C3          out $OUT_ROOT/relax_ac_17000.dat \n"
#        done
#
#        echo -e "$C0$C3" | $CPPTRAJ \
#         > >(tee $OUT_ROOT/cpptraj_3.out) \
#        2> >(tee $OUT_ROOT/cpptraj_3.err >&2)

        # Commands 4: iRED
        C4=""
        for INDEX in ${INDEXES[@]}; do
            C4="$C4 vector v$INDEX ired :$INDEX@N :$INDEX@H \n"
        done
        C4="$C4 matrix ired name ired_matrix order 2 \n"
        C4="$C4 analyze matrix ired_matrix name ired_vectors"
        C4="$C4         vecs ${#INDEXES[@]}"
        C4="$C4         out $OUT_ROOT/ired_vectors.dat \n"
        for FREQ in ${FREQUENCIES[@]}; do
            for CUTOFF in ${CUTOFFS[@]}; do
                C4="$C4 ired relax freq ${FREQ} NHdist 1.02 order 2 norm"
                C4="$C4      tstep 1 tcorr $CUTOFF"
                C4="$C4      modes ired_vectors"
                printf -v CUTOFF "%05d" $CUTOFF
                C4="$C4      out            $OUT_ROOT/ired_${FREQ}_${CUTOFF}_corr.dat"
                C4="$C4      noefile        $OUT_ROOT/ired_${FREQ}_${CUTOFF}_relax.dat"
                C4="$C4      orderparamfile $OUT_ROOT/ired_${FREQ}_${CUTOFF}_order.dat \n"
            done
        done

        echo -e "$C0$C4" | $CPPTRAJ \
         > >(tee $OUT_ROOT/cpptraj_4.out) \
        2> >(tee $OUT_ROOT/cpptraj_4.err >&2)

    done
done
