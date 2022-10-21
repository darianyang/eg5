## V00: using a 1D pcoord of ADP + Mg RMSD
* With MAB from 0-9Å, 8 rectilinear bins from 9-11Å, and MAB from 11-inf
* Going back to 50ps tau
* using 10 initial and later MAB bins
    * after initial run (low bin occupancy (3)), switched to 20 and 20 MAB bins
* using 5 walkers per bin
    * need more sampling, after 500 iterations, started using 10 walkers per bin

## V01: same as V00 but without the middle rectilinear bins
* Started with bins at 0-10 and 10+
* After 294 iterations, decided to switch to 7.5Å for the interface
