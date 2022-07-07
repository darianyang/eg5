#!/bin/bash
# multi_extend.sh
# extend prod sim for multiple (5) replicates of each system

#SYSTEMS=(wt w4f w5f w6f w7f)
SYSTEMS=(w4f)
FF="ipq"
PDB="3k0n"
cd $FF &&

SLURM=h2p_1gpu_prod_07.slurm

# go to target directory
for SYS in ${SYSTEMS[@]} ; do
    cd $SYS
    echo "EXTENDING SYSTEM : $SYS"

    # make the 5 replicate sub-directories
    for V in {00..02..1} ; do
        cd v$V
        echo "GENERATING REPLICA V$V"

        # set up sub-directory slurm file
        cp -v ../../../std_sim_scripts/$SLURM .

        # apply globally to slurm file
        sed -i "s/PDB_TEMP/${PDB}_${SYS}/g" $SLURM
        sed -i "s/VER/${V}/" $SLURM

        # submit the extension prod run
        sbatch $SLURM
        cd ..
    done

    echo "FINISHED SYSTEM : $SYS"
    cd ..
done        

