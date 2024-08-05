#!/bin/bash

for DIR in */ ; do

rsync -axvhP dty7@ls6.tacc.utexas.edu:/home1/09416/dty7/scratch/eg5/ls6-we/${DIR}west.h5 $DIR

done

