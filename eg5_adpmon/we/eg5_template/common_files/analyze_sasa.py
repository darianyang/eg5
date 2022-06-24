"""
Calculate SASA during an WESTPA simulation.
"""

import numpy as np
import mdtraj
import sys

def calc_sasa(selection):
    """
    Parameters
    ----------
    selection : str
        MDTraj style selection string for decomposed SASA region.

    Returns
    -------
    sasa : float
        SASA (A^2) of the selected region.
    """
    # set up data path from args
    wpath = sys.argv[1]

    # load the parent rst frame and current segment coords
    traj_parent = mdtraj.load(f"{wpath}/parent.ncrst", top="eg5_holo.prmtop") 
    traj_segment = mdtraj.load(f"{wpath}/seg.nc", top="eg5_holo.prmtop")
    traj = traj_parent+traj_segment

    #traj = mdtraj.load(f"{wpath}/eg5_2022.pdb")
    #traj = mdtraj.load("../bstates/bstate.ncrst", top="eg5_holo.prmtop")

    # select a portion of the entire system
    # note that some amber models (e.g. TIP4P) will have virtual sites (VS)
    # these have no dedicated radii in the MDTraj _ATOMIC_RADII library
    # https://github.com/mdtraj/mdtraj/blob/master/mdtraj/geometry/sasa.py#L39
    # https://github.com/mdtraj/mdtraj/issues/1151
    sel = traj.topology.select(selection + " and not element VS")

    # slice the trajectory to get atoms coordinates for the selection
    sel_atoms = traj.atom_slice(sel)
    
    # calculate sasa and sum to get total per frame
    sasa = mdtraj.shrake_rupley(sel_atoms, mode="atom").sum(axis=1) 
    
    # convert nm^2 to Angstroms^2 and return
    return sasa * 10**2 

receptor = calc_sasa("resid 23 to 29 or resid 71 to 77 or resid 80 \
                     or resid 102 to 117 or resid 131 or resid 231 \
                     to 232 or resid 264 to 269 or resid 273 or \
                     resid 332 or resid 334 to 337")
adp = calc_sasa("resname ADP")
#adp = calc_sasa("resid 369")

np.savetxt("rec_sasa.dat", receptor) 
np.savetxt("lig_sasa.dat", adp)
