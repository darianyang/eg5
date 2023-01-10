#Mixed mab on bridges

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

----------------------
# multi mab on h2p

## v00: using 5 mab bins
* using mixed mab v01 as template
* changed to 4 walkers per bin
* 5 mab bins: [0, 3.5, 7.5, 11.5, 15.5, 'inf']
    * each with 3 mab bins
    * no forward directionality from 7.5 to 15.5
* recycling at >6A min contact distance 

## v01: using v00 as template
* trying again but with 3 mab bins instead
* [0, 7.5, 15.5, 'inf']
    * with middle bin no directionality
