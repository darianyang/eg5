### moved eg5 from h2p to ls6

multi_mab_dl5_v01:
    * using the multi_mab_dl5_v00 from h2p-we as a base
    * updating the bstate with the corrected dl5 equilibrated structure
    * updating env.sh and node.sh and slurm sub files to ls6 specific
    * I didn't notice before but I had the weight thresholds on previously for some reason
        * turned these off for this next run
    * had to also update the prmtop file to match the updated bstate
    * also changed walkers per bin from 4 to 3 since there are 3 gpus per node on ls6
