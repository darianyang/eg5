import numpy
import sys

fin = str(sys.argv[1])
fout = str(sys.argv[2])

ene_cplx = numpy.loadtxt(fin, usecols=0)
ene_pro = numpy.loadtxt(fin, usecols=1)
ene_lig = numpy.loadtxt(fin, usecols=2)

ene_int = ene_cplx - (ene_pro + ene_lig)

ene_int *= -1

arr = numpy.column_stack(numpy.array([ene_int]))

numpy.savetxt(fout, arr)
