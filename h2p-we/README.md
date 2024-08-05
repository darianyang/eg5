1d_lig_rms_multi_mab_v00
## v00: using 5 mab bins
* using mixed mab v01 as template
* changed to 4 walkers per bin
* 5 mab bins: [0, 3.5, 7.5, 11.5, 15.5, 'inf']
    * each with 3 mab bins
    * no forward directionality from 7.5 to 15.5
* recycling at >6A min contact distance 

1d_lig_rms_multi_mab_v01
## v01:
* trying again but with 3 mab bins instead
* [0, 7.5, 15.5, 'inf']
    * with middle bin no directionality



* wevo_nomon_v01: wevo with eg5 + adp and no monastrol
* this used wevo v03 as a template but I want to replace the basis state with the wt eg5, so no monastrol, maybe it will be easier to unbind
* replaced bstate
* added common files from atb: /jet/home/dyang4/mcb180038p/atb43/eg5/eg5_holo_test_r02b2d_ene_dist_20/common_files
* updated get_pcoord and prmtop names
* updated runseg prmtop names
* note that the true eg5 wt from atb is the same numbering as 1x88, mon is just after it (371)
* updated west.cfg

* wevo_nomon_v01:
    * replaced bstate with i1500 s3 from wevo no mon v00

# multi-mab nomon (WT Eg5)
multi-mab_nomon_v00:
    * from wevo_nomon_v00
    * 4 MAB bins, 4 pcoords 
    # ADP-Mg RMSD (A)
    - [0, 7.5, 10.5, 'inf']
    # interaction energy ADP+Mg and Eg5
    - ['-inf', 'inf']
    # min contact distance ADP and Eg5
    - [0, 6, 'inf']
    # min contact distance PO4 and Eg5
    - [0, 6, 'inf']
multi-mab_nomon_v01:
    * from v00
    * realized I didn't put a limit for int ene pcoord
    * also adjusted to make mab directionf or int ene -1 except below 0 int ene
    * 4 MAB bins, 4 pcoords 
    # ADP-Mg RMSD (A)
    - [0, 7.5, 10.5, 'inf']
    # interaction energy ADP+Mg and Eg5
    - ['-inf', 0, 'inf']
    # min contact distance ADP and Eg5
    - [0, 6, 'inf']
    # min contact distance PO4 and Eg5
    - [0, 6, 'inf']
    * after 284i, changed int ene boundary from 0 to 50, ran 24 hours 8 A100s
        * so now after <50 int ene pcoord will switch to RMSD and PO4 min contact
    * after ~324i, seems like the int ene is progressing well but a few low probability walkers
        * ran another 24 hour job on 8 A100s
        * this time with small and large weight thresholds (similar to MLS setup)
            * small = 1e-129 (300kT) and large = 0.1
        * also, switched to <10 int ene boundary, after this goes to a min PO4 distance and int ene
            * before I used min PO4 and RMSD, but probably no need for RMSD after that point

    * started a WESS run off of this run (from iteration 585)
        * wess_v00_multi-mab_nomon_v01
            * copied everything but left all iterations except for the final i585 files
        * binning off of interaction energy only, every 10 kcal from -200 to 320
            * got some empty bins with this, trying using 50 kcal from -150 to 300
            * made big bins in negative int ene region: -100 to 0 and less than -100
        * using MLS default 0.75 window size
            * changed to 1 to try to populate bins better
                * this didn't help, changed back

    * attempting to do haMSM post-analysis on this dataset

# multi-mab with monastrol
multi-mab_wmon_v00:
    * same setup as multi-mab_nomon_v01, but with monastrol bound
    * so I had to switch the bstate and the common files
        * then update these in the runseg and get_pcoord scripts
        * also added back aux data for monastrol RMSD, SASA, and energy calculations
        * note I had to change analyze_sasa py script to 1x88_solv as well

# multi-mab with DL5
multi-mab_dl5_v00:
    * same setup as other multi-mab runs but with no monastrol and with a deletion of loop5
    * so switched the bstate and common files to dl5

# post TACC
multi-mab_wmon_v01:
    * replicate of v00
    * not working too well
multi-mab_wmon_v02:
    * getting rid of the PO4 dimension
    * so now, using a 3D pcoord: RMSD, mindist, IntEne
    * trying with a larger focus on IntEne
    * also trying with the 86 direction 
    * didn't really go anywhere after a day
multi-mab_wmon_v03:
    * same as v02 set up but I'm using the original directionality settings of v01 
