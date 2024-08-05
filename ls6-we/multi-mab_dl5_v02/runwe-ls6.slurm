#!/bin/bash
#SBATCH -J v02_dl5_eg5_multi_mab
#SBATCH -o job_logs/slurm_%j.out
#SBATCH -e job_logs/slurm_%j.err
#SBATCH -t 48:00:00
#SBATCH -N 3
#SBATCH --ntasks-per-node 3
#SBATCH -p gpu-a100
#SBATCH -A MCB23004
#SBATCH --mail-user=dty7@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL
###SBATCH --reservation=root_45

echo "START `date +"%Y"-%m-%dT%T` `date +"%s"`" 
echo "JobID:    ${SLURM_JOB_ID}"
echo "User:     ${USER}"
echo "Hostlist: ${SLURM_NODELIST}"
#-------------------------

# WE INIT & RUN

  source env.sh || exit 1

  cd $WEST_SIM_ROOT

  SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$SLURM_JOBID.json

  if [ ! -f "west.h5" ] ; then
      source init.sh
  fi

# start server
  #_run --debug --work-manager=zmq \
#  w_run --verbose --work-manager=zmq \
#      --n-workers=0 \
#      --zmq-mode=master \
#      --zmq-write-host-info=$SERVER_INFO \
#      --zmq-comm-mode=ipc &> west-$SLURM_JOBID.log &
w_run --zmq-startup-timeout 3600 \
    --zmq-shutdown-timeout 360 \
    --zmq-timeout-factor 240 \
    --zmq-master-heartbeat 3 \
    --zmq-worker-heartbeat 120 \
    --work-manager=zmq \
    --n-workers=0 \
    --zmq-mode=master \
    --zmq-write-host-info=$SERVER_INFO \
    --zmq-comm-mode=ipc &> west-$SLURM_JOBID.log &
#    --zmq-comm-mode=ipc > west.log


  echo " -> Finished launching master."

# wait on master to create SERVER_INFO file
  for ((n=0; n<60; n++)); do
    if [ -e $SERVER_INFO ] ; then
        echo "== server info file $SERVER_INFO =="
        cat $SERVER_INFO
        break
    fi
    sleep 1
  done

  echo " -> Found SERVER_INFO file, launching clients."

# exit if host info file doesn't appear in one minute
  if ! [ -e $SERVER_INFO ] ; then
    echo 'server failed to start'
    exit 1
  fi

  echo " -> Found SERVER_INFO file, launching clients."

  scontrol show hostname $SLURM_NODELIST >& SLURM_NODELIST.log

# start clients on nodes, with the proper number of cores on each
  for node in $(scontrol show hostname $SLURM_NODELIST); do
    echo " -> Launching client on $node."
    export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
    ssh -o StrictHostKeyChecking=no $node $PWD/node.sh \
        $PWD \
        $SLURM_JOBID \
        $node \
        $CUDA_VISIBLE_DEVICES \
          --verbose \
          --work-manager=zmq \
	  --n-workers=$SLURM_NTASKS_PER_NODE \
          --zmq-mode=client \
          --zmq-read-host-info=$SERVER_INFO \
          --zmq-comm-mode=ipc &> job_logs/node.log$$ &
  done
  echo " -> Launching Completed... ."

  wait

#-------------------------

echo "END `date +"%Y"-%m-%dT%T` `date +"%s"`" 

exit
