import numpy as np
import sys

# grab args
fin = str(sys.argv[1])
fout = str(sys.argv[2])

# using cpptraj energy output file
energies = np.loadtxt(fin)

# needs to be commonly indexable for 1d (e.g. pcoord 1 frame) or 2d (e.g. runseg n frames)
energies = np.atleast_2d(energies)

# first col is frame, then every 9 columns is total energy 
ene_complex = energies[:,9]
ene_recept = energies[:,18]
ene_ligand = energies[:,27]

# interaction energy calculation
ene_int = ene_complex - (ene_recept + ene_ligand)

# return interaction energy as a positive value
ene_int *= -1

# goes from 1 column of seg ene vals to 1 row
# makes it compatible with west_return / west.h5
ene_int = np.column_stack(ene_int)

# saveout
np.savetxt(fout, ene_int)
