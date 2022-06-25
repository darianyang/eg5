#!/bin/bash
# setup and run multiple we runs from template
# or extend multiple runs

# replicates
for V in {00..04} ; do
    # create new dir if needed
    if [ ! -d v$V ] ; then
        cp -r template v$V &&
        sed -i "s/VER/v$V/" v$V/runwe.slurm
    fi
    cd v$V

    bash init.sh &&
    sbatch runwe.slurm

    echo "FINISHED v$V"

    cd ..
done
