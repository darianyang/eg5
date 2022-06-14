#!/bin/bash
# multi_set_and_run.sh
# set up and run multiple (5) replicates of each system

SYSTEMS=(noions addions1 addions2 ipq-noions ipq-addions1 ipq-addions2)
PDB=1ii6

# go to target directory
for SYS in ${SYSTEMS[@]} ; do
    # make parent dir if needed
    if [ ! -d $SYS ] ; then
        mkdir $SYS
    fi
    cd $SYS
    echo -e "\nRUNNING SYSTEM : $SYS"


    cp -v ../../${PDB}_leap.pdb .
    cp ../sim_template/* .

    # formatting
    bash temp_sed.sh ${PDB} v00
    # make the inital parm and crd files
    tleap -f ${SYS}.in > tleap.out
    # submit the prep run, which submits the prod after finishing
    #sbatch prep_mpi.slurm
    sbatch min.slurm

    echo "FINISHED SYSTEM : $SYS"
    cd ..
done        

