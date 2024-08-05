#!/bin/bash
#SBATCH --job-name=v01_mon_eg5_multi_mab
#SBATCH --output=job_logs/slurm.out
#SBATCH --error=job_logs/slurm.err
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --gres=gpu:8
#SBATCH --cluster=gpu
#SBATCH --partition=a100_nvlink
#SBATCH --time=71:59:59

set -x
cd $SLURM_SUBMIT_DIR
source env.sh || exit 1

env | sort

cd $WEST_SIM_ROOT
SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$SLURM_JOBID.json

which w_run

# start server
w_run --work-manager=zmq --n-workers=0 --zmq-mode=master --zmq-write-host-info=$SERVER_INFO --zmq-comm-mode=tcp &> west-$SLURM_JOBID.log &

# wait on host info file up to one minute
for ((n=0; n<60; n++)); do
    if [ -e $SERVER_INFO ] ; then
        echo "== server info file $SERVER_INFO =="
        cat $SERVER_INFO
        break
    fi
    sleep 1
done

# exit if host info file doesn't appear in one minute
if ! [ -e $SERVER_INFO ] ; then
    echo 'server failed to start'
    exit 1
fi

# start clients, with the proper number of cores on each

for node in $(scontrol show hostname $SLURM_NODELIST); do
    ssh -o StrictHostKeyChecking=no $node $PWD/node.sh $SLURM_SUBMIT_DIR $SLURM_JOBID $node $CUDA_VISIBLE_DEVICES --work-manager=zmq --n-workers=$SLURM_NTASKS_PER_NODE --zmq-mode=client --zmq-read-host-info=$SERVER_INFO --zmq-comm-mode=tcp &
done


wait
