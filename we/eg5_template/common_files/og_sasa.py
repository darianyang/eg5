import numpy
import mdtraj
import sys

#traj_parent = mdtraj.load("parent.ncrst", top="./eg5_holo.prmtop") # load the parent
#traj_segment = mdtraj.load("seg.nc", top="./eg5_holo.prmtop") # load the segment coords
#traj = traj_parent+traj_segment

traj = mdtraj.load("eg5_2022.pdb")

rec_sel = numpy.concatenate((traj.topology.select('resid 23 to 29'),
                             traj.topology.select('resid 71 to 77'),
                             traj.topology.select('resid 80'),
                             traj.topology.select('resid 102 to 117'),
                             traj.topology.select('resid 131'),
                             traj.topology.select('resid 231 to 232'),
                             traj.topology.select('resid 264 to 269'),
                             traj.topology.select('resid 273'),
                             traj.topology.select('resid 332'),
                             traj.topology.select('resid 334 to 337'))) # selection of receptor atoms

lig_sel = traj.topology.select('resid 369')

rec_atoms = traj.atom_slice(rec_sel) # slice the trajectory to get atoms coordinates for the selection

lig_atoms = traj.atom_slice(lig_sel) # slice the trajectory to get atoms coordinates for the selection

rec_sasa = mdtraj.shrake_rupley(rec_atoms, mode="atom").sum(axis=1) # calculate sasa and sum to get total per frame

lig_sasa = mdtraj.shrake_rupley(lig_atoms, mode="atom").sum(axis=1) # calculate sasa and sum to get total per frame

rec_sasa *= 10 # convert to Angstroms

lig_sasa *= 10 # convert to Angstroms

numpy.savetxt("rec_sasa.dat", rec_sasa) # save to text file

numpy.savetxt("lig_sasa.dat", lig_sasa) # save to text file
