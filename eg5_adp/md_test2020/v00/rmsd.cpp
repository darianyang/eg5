parm eg5_2020_solv.prmtop
trajin 03_eq1.nc
trajin 04_eq2.nc
trajin 05_eq3.nc
autoimage anchor :113
reference eg5_2020_solv.pdb [REF]
rms complex :1-368,370&!@H= out equil/rms_complex_heavy.dat ref [REF]
rms prot :1-368&!@H= out equil/rms_prot_heavy.dat ref [REF] nofit 
rms adp :370&!@H= out equil/rms_adp_heavy.dat ref [REF] nofit
go
quit
