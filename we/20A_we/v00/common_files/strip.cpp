parm 1x88_solv.prmtop
trajin 1x88_solv.inpcrd
strip :WAT,Na+,Cl- parmout 1x88_dry.prmtop
strip :WAT,Na+,Cl-,Mg2+,ADP,MON parmout 1x88_dry_apo.prmtop
run
quit
