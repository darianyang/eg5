V01:
* Switched 100ps tau to 50ps
* Added int ene to wevo input
    * Now it is 3d input with adp rms, min contact dist, int ene

V02:
* Now using 40 segments instead of 100
* Also starting with 200 kT prob min (1e-84)
* This version is a copy of ATB resampler coordinates
    * Coord 1: ADP RMSD
    * Coord 2: Receptor Int Ene where receptor is all protein and ligand is Mg and ADP
    * Coord 3 (1 bin) : min contact distance

