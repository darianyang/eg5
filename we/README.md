# WESTPA Simulations
* See 20A_we for previous WESTPA runs using the 1x88 system with a 20A box size
    * current we simulation directories are using the 16A box equilibrated 1x88

* TEMPLATE: 2D ss WE simulations of ADP unbinding from eg5 (1x88)
    * pcoords: ADP Base RMSD, ADP Phos RMSD
    * v00-v04 and thresh12/22
        * TODO: threshold for max as well at 0.1
    * 50 ps tau
* TEMPLATE2: 1D ss WE simulation of ADP unbining from eg5 (1x88)
    * pcoord: Mg2+ RMSD
    * 100ps tau
    * auto strip water but not ions
    * 8 walkers and 20 MAB bins
* ALL TEMPLATES:
    * no HDF5 framework for iteration.h5 files
    * 10 initial segments per state in w_init

* 1d_lig_rms_thresh12
    * same as template 2 but with 1D pcoord of ADP and Mg RMSD
    * also now has 1e-12 threshold for min and 0.1 for max

* 1d_lig_rms_from_106_1:
    * Took a later adp base moved segment and started new WE simulation
    * Constrained MAB to be from 8Å pcoord (ADP and Mg RMSD) to inf
    * Also put more 1D MAB bins (20 --> 25)

Note: to turn off hdf5 framework for iteration h5 files, remove the datarefs/iteration line
