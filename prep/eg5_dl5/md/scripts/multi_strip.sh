#!/bin/bash
# strip water and ions from multiple (5) replicates of each system

#SYSTEMS=(wt w4f)
SYSTEMS=(w6f w7f)
PDB=3k0n
FF="ipq"
CPPTRAJ=cpptraj.MPI
export DO_PARALLEL="mpirun -np 8"

cd $FF

# go to target directory
for SYS in ${SYSTEMS[@]} ; do
    cd $SYS &&
    echo "RUNNING SYSTEM : $SYS"

    # make the 5 replicate sub-directories
    for V in {00..00..1} ; do
        echo "GENERATING REPLICA V$V"
        # go to sub-directory
        cd v$V

        for TRAJ in {07..08} ; do
            # Commands 00: Create a water stripped but with ions pdb for ref
            #C00="     parm ${PDB}_${SYS}_solv.prmtop \n"
            C00="     parm ../${PDB}_${SYS}_dry.prmtop \n"
            C00="$C00 trajin ${TRAJ}_prod.nc \n"
            C00="$C00 strip :WAT,Cl- \n"
            C00="$C00 trajout ${TRAJ}_prod_dry.nc \n"
            echo -e "$C00" > cpp_strip.in 
            $DO_PARALLEL $CPPTRAJ -i cpp_strip.in
    
        done        

        cd ..
    done
    cd ..
    echo "FINISHED RUNNING SYSTEM : $SYS"
done
