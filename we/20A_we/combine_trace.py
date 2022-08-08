#!/bin/env python3
"""
Combine and analyze traj seg seg.nc files in MDAnalysis.
"""

import MDAnalysis as mda
import numpy as np
import h5py
import os

def get_parents(walker_tuple, h5_file):
    it, wlk = walker_tuple
    parent = h5_file["iterations/iter_{:08d}".format(it)]["seg_index"]["parent_id"][wlk]
    return it-1, parent

def trace_walker(walker_tuple, h5_file):
    """
    Parameters
    ----------
    walker_tuple : tuple
        A tuple of the desired iteration and walker to trace from.
        e.g. (100,1) would trace from iteration 100, walker 1.
    h5_file : str
        Path to west.h5 file from WESTPA simulation. 

    Returns
    -------
    trace : ndarray
        Array of 2 element arrays (interation, walker) tracing up to
        the target walker from iteration 1.

    """
    h5 = h5py.File(h5_file, "r")
    # Unroll the tuple into iteration/walker 
    it, wlk = walker_tuple
    # Initialize our path
    path = [(it,wlk)]
    # And trace it
    while it > 1: 
        it, wlk = get_parents((it, wlk), h5)
        path.append((it,wlk))
    trace = np.array(sorted(path, key=lambda x: x[0]))
    return trace


def concat_iter(trace, top, out_path=None):
    """
    Parameters
    ----------
    trace : ndarray
        Array of 2 element arrays (iteration, walker) tracing up to 
        the target walker from iteration 1.
    top : str
        Path string to topology/parameter file.
    out_path : str
        Optional path string for trajectory output.    
        File extension must be MDA compatible (e.g. *.ncdf).

    Returns
    -------
    trace_traj : MDA Universe object
         Single continuous trajectory of traced iteration segments.
    """
    # paths to each iteration walker in the trace
    paths = [f"traj_segs/{it:06d}/{wlk:06d}/seg.nc" for it, wlk in trace]    

    # load all segment trajectories for trace into Universe object
    #trace_traj = mda.Universe(top, paths, in_memory=True, in_memory_step=10)
    trace_traj = mda.Universe(top, paths)

    if out_path:
        #protein = trace_traj.select_atoms("protein")
        protein = trace_traj.select_atoms("all")
        with mda.Writer(out_path, protein.n_atoms) as W:
            for ts in trace_traj.trajectory:
                W.write(protein)

    return trace_traj

#WE = "2kod_we_1d_sim_v03/v03.02"
#WE = "2kod_we_1d_sim_v02"
WE = "v00"
os.chdir(WE)

trace = trace_walker((180, 7), "west.h5")
concat_iter(trace, f"common_files/1x88_solv.prmtop", out_path="traced_traj/180_7.ncdf")

# TODO: need to use cpptraj autoimage and maybe rms fit to tame the unweildly output
# add this to the trace walker function through os commands?
# or use mda equivalent?
