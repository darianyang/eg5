#!/bin/bash -l

#set -x

#cd $WEST_SIM_ROOT
cd $1; shift
source env.sh
export WEST_JOBID=$1; shift
export SLURM_NODENAME=$1; shift
export CUDA_VISIBLE_DEVICES_ALLOCATED=$1; shift
# export LOCAL=/local/$WEST_JOBID
echo "starting WEST client processes on: "; hostname
echo "current directory is $PWD"
echo "environment is: "
env | sort
#env | sort > ENV.node.sh.out$$

echo "CUDA_VISIBLE_DEVICES = " $CUDA_VISIBLE_DEVICES


#comment out-- /usr/lib64 should not be put before ...westpaenv...lib
# so that  GLIBCXX doesn't occur
#export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH

echo "w_run   $@"   >    client_w_run_cmd

#     w_run  "$@" |& tee client_w_run_cmd.log
#     w_run  "$@" &> west-$SLURM_NODENAME-node.log

#      w_run  "$@"
w_run "$@" &> west-$SLURM_NODENAME-node.log

echo "Shutting down.  Hopefully this was on purpose?"
