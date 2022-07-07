#!/bin/bash
# multi_set_and_run.sh
# set up and run multiple (5) replicates of each system

SYSTEMS=(v00)
PDB=1x88

# go to target directory
for SYS in ${SYSTEMS[@]} ; do
    # make parent dir if needed
    if [ ! -d $SYS ] ; then
        mkdir $SYS
    fi
    cd $SYS
    echo -e "\nRUNNING SYSTEM : $SYS"

    cp -v ../../${PDB}_leap.pdb .
    #cp -v ../../${PDB}_solv.* .
    cp ../sim_template/* .

    # formatting
    bash temp_sed.sh ${PDB} v00
    # make the inital parm and crd files
    tleap -f tleap.in > tleap.out
    # submit the prep run, which submits the prod after finishing
    sbatch prep_mpi.slurm

    echo "FINISHED SYSTEM : $SYS"
    cd ..
done        

