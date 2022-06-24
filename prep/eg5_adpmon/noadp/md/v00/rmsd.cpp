parm 1x88_apo_solv.prmtop
trajin 03_eq1.nc 1 last 100
trajin 04_eq2.nc 1 last 100
trajin 05_eq3.nc 1 last 100
trajin 25ns.nc 1 last 100
autoimage
reference 1x88_apo_solv.pdb [REF]
rms prot :1-368&!@H= out equil/rms_prot_heavy.dat ref [REF]
#rms complex :1-368,370-371&!@H= out equil/rms_complex_heavy.dat ref [REF] nofit
#rms adp :370&!@H= out equil/rms_adp_heavy.dat ref [REF] nofit
#rms mon :371&!@H= out equil/rms_mon_heavy.dat ref [REF] nofit
rms mon :369&!@H= out equil/rms_mon_heavy.dat ref [REF] nofit
# rmsf
average crdset avg_crd
run
rms :1-368@CA,C,O,N ref avg_crd
atomicfluct out equil/rmsf.dat :1-368@CA,C,O,N byres
rms protcut :18-364&!@H= out equil/rms_protcut_heavy.dat ref [REF]
#rms complexcut :18-364,370-371&!@H= out equil/rms_complexcut_heavy.dat ref [REF] nofit
go
quit
