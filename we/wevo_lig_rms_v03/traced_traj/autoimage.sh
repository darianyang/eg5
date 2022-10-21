#!/bin/bash

ITER=$1
TRAJ=$2


cat << EOF > autoimage.cpp
parm ../common_files/1x88_dry.prmtop
trajin ${ITER}_${TRAJ}.dcd
autoimage
trajout ${ITER}_${TRAJ}_auto.dcd
go
quit
EOF

cpptraj -i autoimage.cpp

rm -v autoimage.cpp
