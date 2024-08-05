### moved eg5 from h2p to ls6

multi_mab_dl5_v01:
    * using the multi_mab_dl5_v00 from h2p-we as a base
    * updating the bstate with the corrected dl5 equilibrated structure
    * updating env.sh and node.sh and slurm sub files to ls6 specific
    * I didn't notice before but I had the weight thresholds on previously for some reason
        * turned these off for this next run
    * had to also update the prmtop file to match the updated bstate
    * also changed walkers per bin from 4 to 3 since there are 3 gpus per node on ls6
    * at ~1900 iterations, updated the binning at the final multi-mab bin to focus more on int ene and less on the PO4 distance

multi_mab_nomon_v01_r00 and r01:
    * using multi_mab_nomon_v01 from h2p-we as template
    * just updating to be compatible with ls6, env.sh, node.sh and slurm submission
    * didn't change anything else, should have changed walkers per bin to 3 but kept it the same as before for consistency since this is a replicate

multi_mab_dl5_v02:
    * using same settings as before (from adjusted bins of v01)
    * back to simulating with 4 walkers per bin to be similar to WT

multi_mab_wmon_v01:
    * replicate of the w_mon run from h2p
