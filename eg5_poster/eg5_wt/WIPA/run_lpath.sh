#!/bin/bash

# assign source and target states using the phi/psi angles from WE simulations
#lpath discretize -we -W ./multi.h5 --assign-arguments="--config-from-file --scheme C7_EQ -W multi.h5"

# extract all successful pathways connecting source and target states
# also extract the cluster labels for matching in the next step
#lpath extract -we -W west_labeled.h5 -A ./FULLY_UNBOUND2/assign.h5 -ss 0 -ts 1 -p -a labels --stride 5 --threads 4 #--last-iter 456 --trace-basis

# perform matching with condensing repeat pairs
# uses the cluster labels as states through reassign_custom
lpath match -we -ra reassign_custom.reassign_custom -op succ_traj/reassigned.pickle --condense 2 --stats

#lpath match -we --input-pickle succ_traj/output.pickle --output-pickle succ_traj/match-output.pickle --cluster-label-output succ_traj/cluster_labels.npy --export-h5 --file-pattern "west_succ_c{}.h5" --reassign-method "reassign_segid" --substring

#lpath plot --plot-input succ_traj/reassigned.pickle
