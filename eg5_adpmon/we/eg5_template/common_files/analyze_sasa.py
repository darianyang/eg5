"""
Calculate SASA during an WESTPA simulation.
"""

import numpy as np
import mdtraj

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

    # load the parent rst frame and current segment coords
    #traj_parent = mdtraj.load("parent.ncrst", top="./eg5_holo.prmtop") 
    #traj_segment = mdtraj.load("seg.nc", top="./eg5_holo.prmtop")
    #traj = traj_parent+traj_segment

    # test
    traj = mdtraj.load("eg5_2022.pdb")    

    # select receptor and ligand (ADP)
    #    rec_sel = np.concatenate((traj.topology.select("resid 23 to 29"),
    #                                 traj.topology.select("resid 71 to 77"),
    #                                 traj.topology.select("resid 80"),
    #                                 traj.topology.select("resid 102 to 117"),
    #                                 traj.topology.select("resid 131"),
    #                                 traj.topology.select("resid 231 to 232"),
    #                                 traj.topology.select("resid 264 to 269"),
    #                                 traj.topology.select("resid 273"),
    #                                 traj.topology.select("resid 332"),
    #                                 traj.topology.select("resid 334 to 337"))) 
    #    
    #    lig_sel = traj.topology.select("resid 369")

    sel = traj.topology.select(selection)    

    # slice the trajectory to get atoms coordinates for the selection
    sel_atoms = traj.atom_slice(sel)
    
    # calculate sasa and sum to get total per frame
    sasa = mdtraj.shrake_rupley(sel_atoms, mode="atom").sum(axis=1) 
    
    # convert nm^2 to Angstroms^2 and return
    return sasa * 10**2 

receptor = calc_sasa("resid 23 to 29 or resid 71 to 77 or resid 80 or resid 102 to 117 or resid 131 or resid 231 to 232 or resid 264 to 269 or resid 273 or resid 332 or resid 334 to 337")
adp = calc_sasa("resid 369")

np.savetxt("rec_sasa.dat", receptor) 
np.savetxt("lig_sasa.dat", adp)

