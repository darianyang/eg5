import numpy as np
import sys

import h5py
import pickle
import os
import subprocess
import mdtraj as md
import tqdm

import logging
logging.basicConfig(filename="msm_we.log", encoding='utf-8', level=logging.DEBUG) 
log = logging.getLogger()

from msm_we import msm_we
"""
This code
1. Initializes a model with the given PDBs and file globs
2. Opens the associated H5 file, in append mode
3. Reads through each of the trajectory files, putting the coordinates into the H5 file's auxdata/coord (if they don't already exist)
"""

we_h5filename = sys.argv[1]
refPDBfile = sys.argv[2]
WEfolder = sys.argv[3]

log.info('\n\n\nSTART COLLECT COORDS: Preparing coordinates...')
log.debug(f'Doing collectCoordinates on WE file {we_h5filename}')

# open west.h5 file in append mode (read/write but create if it doesn't exist)
h5file = h5py.File(we_h5filename, "a")
# set max_iters from west.h5 (note, empty iter for last val (not needed in loop))
maxIter = h5file.attrs["west_current_iteration"]
# set pdb file based mdtraj object
pdb = md.load(refPDBfile, top=refPDBfile)
# define coord name to look for in each traj_segs/ sub-dir
traj = "seg-nowat.ncrst"

n_iter = None
for n_iter in tqdm.tqdm(range(1, maxIter)):

    # get n_segs for each iteration
    nS = h5file[f"iterations/iter_{n_iter:08d}/seg_index"][:].shape[0] 

    # initialize coord datasets
    coords = np.zeros((nS, 2, pdb.n_atoms, 3))
    dsetName = f"/iterations/iter_{int(n_iter):08d}/auxdata/coord"

    coords_exist = False
    try:
        dset = h5file.create_dataset(dsetName, np.shape(coords), 
                                     compression=4, scaleoffset=6, chunks=True)
    except (RuntimeError, ValueError) as e:
        log.debug(f"{e} : coords exist for iteration {str(n_iter)} NOT overwritten")
        coords_exist = True
        continue

    for iS in range(nS):
        # define path of current segment
        childPath = WEfolder + f"/traj_segs/{n_iter:06d}/{iS:06d}"
        # define *parent* path of current segment
        parent_iS = int(h5file[f"iterations/iter_{n_iter:08d}/seg_index"]["parent_id"][iS])
        #log.debug(f"ITER: {n_iter} | SEG: {iS} | PARENT: {parent_iS}")
        # recycled segments have a negative parent_id, use ref for these parents
        # note only sometimes negative parent_id since most of the time the 
        # recycled traj is low weight and gets merged
        if parent_iS < 0:
            coord0 = np.squeeze(pdb._xyz)
            log.info(f'Negative parent_id for ITER:{n_iter} SEG:{iS} (recycled)')
        else:
            # if not recycled, can just create path var
            parentPath = WEfolder + f"traj_segs/{n_iter-1:06d}/{parent_iS:06d}"

            try:
                coord0 = np.squeeze(md.load(f'{parentPath}/{traj}', top=pdb.topology)._xyz)
            except OSError as e:
                log.warning(f"{e} : Parent traj file doesn't exist for iteration {str(n_iter)}, loading reference structure coords")
                coord0 = np.squeeze(pdb._xyz)

        coord1 = np.squeeze(md.load(f'{childPath}/{traj}', top=pdb.topology)._xyz)

        coords[iS, 0, :, :] = coord0
        coords[iS, 1, :, :] = coord1

    if not coords_exist:
        dset[:] = coords

log.debug(f"Wrote coords for {n_iter} iterations.")
h5file.close()
