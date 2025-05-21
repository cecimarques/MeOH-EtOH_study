Hello :)

Welcome to further supplementary information of our paper xxx.

Inside the directory **Simulation_files** you can find the microstates attained by the end of the production run at ambient conditions for all of the 11 systems studied in our work, namely pure EtOH, pure MeOH and MeOH/EtOH binary mixtures at ratios 1:9, 1:4, 3:7, 2:3, 1:1, 3:2, 7:3, 4:1, 9:1. These are given in the form of LAMMPS data files (atom_style full). As potential parameters are included in the file, files for all four force fields used in the assessment at ambient conditions are provided. (...)

Inside the directory **Post_processing_codes** you can find all python codes written in-house that were used in the post-processing of simulation files aiming to get results. These are Hbond_lifetime_ME.py, (xxx).
- The code _Hbond_lifetime_ME.py_ corresponds to the code used to compute the continuous hydrogen bond life time between hydrogens of methanol and oxygens of ethanol. There was one code for each of the four hydrogen types found in the binary methanol/ethanol mixtures, namely H(MeOH)...O(MeOH), H(MeOH)...O(EtOH), H(EtOH)...O(MeOH), H(EtOH)...O(EtOH), with the code in each of the four cases following the same setup (only minor changes/adaptations are required). The calculation follows the specifics described in the supplementary information of the paper.
