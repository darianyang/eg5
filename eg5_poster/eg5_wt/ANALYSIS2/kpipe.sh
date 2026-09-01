#!/bin/bash
# kinetics_pipeline.sh
# Reproduces the results of w_ipa -ao, but manually with self defined aux values
# w_ipa is essentially: w_assign and w_direct
# w_direct all is: w_kinetics, w_kinavg, and w_stateprobs
# Pipeline of commands to get kinetics information from west.h5 file.
# Run in a directory with target west.h5 file
# currently configured to run kinetics for a 2D scheme with 2 aux values

SCHEME=ANALYSIS2 
WEST=~/we_data/eg5/west.h5.i585

# make scheme dir and fill with current pipeline run file
mkdir -pv $SCHEME
cp -v kpipe.sh $SCHEME
cd $SCHEME

# define bins and states with yaml files
cat << EOF > BINS
---
bins:
    type: RectilinearBinMapper
    # add buffer region: strict starting state definition
    boundaries: [[0.0, 4.0, 6.0, 'inf']]
EOF
cat << EOF > STATES
---
states:
  - label: a
    coords:
      - [1.0]

  - label: b
    coords:
      - [7.0]
EOF

# create module.py file to process 1D or 2D scheme
cat << EOF > module.py
#!/usr/bin/env python

import numpy

def pull_data_1d(n_iter, iter_group):
    '''
    This function reshapes auxiliary data for each iteration and returns it.
    '''
    auxdata = iter_group['pcoord'][...]
    #data = auxdata[:,:,numpy.newaxis]
    data = auxdata[:,:,2] # pcoord min contact
    return auxdata

def pull_data_2d(n_iter, iter_group):
    '''
    This function reshapes 2 auxiliary datasets for each iteration and returns it.
    '''
    #auxdata1 = iter_group['auxdata']['${AUX_A}'][...]
    # special for pcoord
    auxdata1 = iter_group['pcoord'][...]
    auxdata2 = iter_group['pcoord'][...]
    data = numpy.dstack((auxdata1[2]))
    return data
EOF

# run w_assign to assign macrostates based off of defined BINS and STATES
w_assign -W $WEST --bins-from-file BINS --states-from-file STATES -o assign.h5 --construct-dataset module.pull_data_1d --serial &&
echo "    w_assign finished: assign.h5 file created in $SCHEME directory"

# run w_direct to then calculate multiple flux values and rate constants
# 'all' is equivalent to w_kinetics, w_kinavg, and w_stateprobs
# step size of 1 is used, and can later be window averaged based off of ACF decay lag time, or can use default
# w_direct by default will set 'correl' = True, which calculates the ACF lag time and window averages
#w_direct kinetics -a assign.h5 -W $WDIR/$WEST -o direct.h5 --evolution-mode cumulative --step-iter 1 &&

w_direct all -a assign.h5 -W $WEST -o direct.h5 --evolution-mode cumulative --step-iter 1 &&

#w_direct all -a assign.h5 -W $WDIR/$WEST -o direct.h5 --step-iter 1 &&
#w_direct all -a assign.h5 -W $WDIR/$WEST -o direct.h5 &&
echo "    w_direct finished: direct.h5 file created in $SCHEME directory"

