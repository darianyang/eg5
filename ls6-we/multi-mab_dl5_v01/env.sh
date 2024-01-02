# Modify to taste

# Load Modules
ml reset
echo  "BEGIN INSIDE ENV.SH"
#export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/lib64:/opt/apps/gcc/9.4.0/lib64/libstdc++.so.6:$LD_LIBRARY_PATH

ml intel/19.1.1
ml impi/19.0.9
ml cuda/11.3
ml -python3

#base=/work/00770/milfeld/ls6/benchmark/nick2
#base=/work/09067/jml230/ls6/apps/module_build_hdf_speedup
base=/work/09416/dty7/ls6/module_build_we_amber/
ml use $base/modules/modulefiles

#ml conda/4.12.0
#ml envwestpa/1.0-C4120_1

#ml fftw/3.3.10
#ml openmm/7.7.0-EW10_1-F3310_1
#ml swig/4.1.1

#ml boost/1.72.0-EW10_1
#ml amber/22-B1720_1-EW10_1-F3310_1
#ml westpa/A22_1-EW10_1-O770_1

ml conda/4.12.0_a
ml envwestpa/1.0_b

ml fftw/3.3.10_a
ml openmm/7.7.0_a
ml swig/4.1.1_b

ml boost/1.72.0_a
ml amber/22_a
ml westpa/22.06_b

ml
which python
#xport LD_LIBRARY_PATH=/usr/lib64:/opt/apps/gcc/9.4.0/lib64/:$LD_LIBRARY_PATH
#xport LD_LIBRARY_PATH=/work/00770/milfeld/ls6/benchmark/nick2/modules/envwestpa/1.0-C4120_1/lib/:$LD_LIBRARY_PATH
# if alone breaks due to GLIBexport LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

#xport LD_LIBRARY_PATH=/usr/lib64:/work/00770/milfeld/ls6/benchmark/nick2/modules/envwestpa/1.0-C4120_1/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/work/00770/milfeld/ls6/benchmark/nick2/modules/envwestpa/1.0-C4120_1/lib:/usr/lib64:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/work/009067/jml230/ls6/apps/module_build_hdf_speedup/modules/envwestpa/1.0_a/lib:/usr/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/work/09416/dty7/ls6/apps/module_build_hdf_speedup/modules/envwestpa/1.0_b/lib:/usr/lib64:$LD_LIBRARY_PATH

echo $LD_LIBRARY_PATH > see.LD_LIBRARY_PATH_new_westpa

export CUDA_HOME=${TACC_CUDA_DIR}
#export CUDA_VISIBLE_DEVICES=0,1,2,3
export CUDA_VISIBLE_DEVICES=0,1,2

export WEST_SIM_ROOT="$PWD"
export SIM_NAME=$(basename $WEST_SIM_ROOT)
export AMBERHOME=${TACC_AMBER_DIR}
export WEST_ROOT=${TACC_WESTPA_DIR}

# This is our local scratch, where we'll store files during the dynamics.
export NODELOC=$WORK
export USE_LOCAL_SCRATCH=1

# Explicitly name our simulation root directory
if [[ -z "$WEST_SIM_ROOT" ]]; then
    # The way we're calling this, it's $SLURM_SUBMIT_DIR, which is fine.
    export WEST_SIM_ROOT="$PWD"
fi

# Export environment variables for the ZMQ work manager.
export WM_ZMQ_MASTER_HEARTBEAT=100
export WM_ZMQ_WORKER_HEARTBEAT=100
export   WM_ZMQ_TIMEOUT_FACTOR=300
export     BASH=$SWROOT/bin/bash
export     PERL=$SWROOT/usr/bin/perl
export      ZSH=$SWROOT/bin/zsh
#xport IFCONFIG=$SWROOT/bin/ifconfig  #not valid on ls6
export IFCONFIG=$SWROOT/usr/sbin/ifconfig
export      CUT=$SWROOT/usr/bin/cut
export       TR=$SWROOT/usr/bin/tr
export       LN=$SWROOT/bin/ln
export       CP=$SWROOT/bin/cp
export       RM=$SWROOT/bin/rm
export      SED=$SWROOT/bin/sed
export      CAT=$SWROOT/bin/cat
export     HEAD=$SWROOT/bin/head
export      TAR=$SWROOT/bin/tar
export      AWK=$SWROOT/usr/bin/awk
export    PASTE=$SWROOT/usr/bin/paste
export     GREP=$SWROOT/bin/grep
export     SORT=$SWROOT/usr/bin/sort
export     UNIQ=$SWROOT/usr/bin/uniq
export     HEAD=$SWROOT/usr/bin/head
export    MKDIR=$SWROOT/bin/mkdir
export     ECHO=$SWROOT/bin/echo
export     DATE=$SWROOT/bin/date
export   SANDER=$AMBERHOME/bin/sander
export    PMEMD=$AMBERHOME/bin/pmemd.cuda
export  CPPTRAJ=$AMBERHOME/bin/cpptraj
echo  "END   INSIDE ENV.SH"
