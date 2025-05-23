units		real
atom_style	atomic
boundary        p p p

# ------------------------------------------------------------
# This is a LAMMPS data file built using the first configuration of the coarsened trajectory (obtained using
# the python code from the all-atom configuration for each system). 
# Reminder: The MeOH and EtOH molecules are assigned # as "atom types" 1 and 2, respectively.
read_data       configuration.data

mass		1 32
mass		2 46

# The force field is completely made up just for the sake of being able to do the rerun. I am putting whatever:
# it does not matter for the sake of this calculation.
pair_style      lj/cut 12
pair_coeff	1 1 0.6 3.0
pair_coeff	1 2 0.6 3.0
pair_coeff	2 2 0.6 3.0

neighbor        0.5 bin 
neigh_modify    delay 0 every 1 check yes

# Useful information for the setup of these commands: the configurations in the coarsened trajectory are 
# taken as if they were consecutively 1 timestep away from one another when creating the trajectory.
compute         myRDF all rdf 1200 1 1 1 2 2 2
fix             myRDF1 all ave/time 1 1999 1999  c_myRDF[*] file rdf.dat mode vector
rerun		com_traj.out first 0 last 1999 every 1 dump x y z 
