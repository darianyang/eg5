# WESTPA Simulations
* TEMPLATE: 2D ss WE simulations of ADP unbinding from eg5 (1x88)
    * pcoords: ADP Base RMSD, ADP Phos RMSD
    * v00-v04 and thresh12/22
        * TODO: threshold for max as well at 0.1
    * 50 ps tau
* TEMPLATE2: 1D ss WE simulation of ADP unbining from eg5 (1x88)
    * pcoord: Mg2+ RMSD
    * 100ps tau
    * auto strip water but not ions
* 1d_lig_rms_thresh12
    * same as template 2 but with 1D pcoord of ADP and Mg RMSD
    * also now has 1e-12 threshold for min and 0.1 for max

Note: to turn off hdf5 framework for iteration h5 files, remove the datarefs/iteration line
