#!/bin/bash

ITER=180
TRAJ=7


cat << EOF > autoimage.cpp
parm ../common_files/1x88_solv.prmtop
trajin ${ITER}_${TRAJ}.ncdf
autoimage
trajout ${ITER}_${TRAJ}_auto.dcd
go
quit
EOF

cpptraj -i autoimage.cpp

rm -v autoimage.cpp
