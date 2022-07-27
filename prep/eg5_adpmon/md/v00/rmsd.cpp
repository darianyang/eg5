parm 1x88_dry.prmtop
trajin 06_prod.nc 1 last 100
trajin 07_prod.nc 1 last 100
trajin 08_prod.nc 1 last 100
autoimage
reference 1x88_solv.pdb [REF]
rms prot :1-368&!@H= out 1us/rms_prot_heavy.dat ref [REF]
rms complex :1-368,370-371&!@H= out 1us/rms_complex_heavy.dat ref [REF] nofit
rms adp :370&!@H= out 1us/rms_adp_heavy.dat ref [REF] nofit
rms mon :371&!@H= out 1us/rms_mon_heavy.dat ref [REF] nofit
# rmsf
average crdset avg_crd
run
rms :1-368@CA,C,O,N ref avg_crd
atomicfluct out 1us/rmsf.dat :1-368@CA,C,O,N byres
rms protcut :18-364&!@H= out 1us/rms_protcut_heavy.dat ref [REF]
rms complexcut :18-364,370-371&!@H= out 1us/rms_complexcut_heavy.dat ref [REF] nofit
go
quit
