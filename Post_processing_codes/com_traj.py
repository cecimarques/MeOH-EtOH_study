# -*- coding: utf-8 -*-
"""
@author: Cecilia Alvares
"""
# ---------------------------------------------------------------------------
# This code is going to read the com.out file output by running the script
# lammps_rdf-com1.md in LAMMPS and create a coarsened trajectory out of it.

import numpy as np

# ---------------------------------------------------------------------------
# A total of 2000 configurations, saved in the production run, will be used
# here.
number_of_configurations = 2000

# Input here the number of molecules composing the simulation domain for the 
# given system (MeOH and EtOH included).
lines = 729

com_traj = np.zeros((1,4))
ofi2 = open("com.out", 'r')
ofi2.readline()
ofi2.readline()
ofi2.readline()
for it_0 in range (0, number_of_configurations):
    ofi2.readline()
    temporary_array = []
    for it_1 in range (0, lines):
        dump2 = ofi2.readline()
        for e,it_2 in zip(dump2.split(' '), range(4)):
            temporary_array.append(float(e))
    temporary_array = np.array(temporary_array,float)
    temporary_array = temporary_array.reshape(lines,4)
    com_traj = np.concatenate((com_traj, temporary_array), axis = 0)

com_traj = np.delete(com_traj, 0 , 0)

# -----------------------------------------------------------------------
# There will be other auxiliary info that I will need to create the traj.
# This is (1) the dimensions of the simulation box in each configuration 
# read above and (2) a definition of "atom type" that will identify whe-
# ther my molecule is a MeOH or EtOH molecule.

# The first item will be accomplished by simply reading the section of the
# log.lammps file that I've output in the rerun and that contains the info
# on the values of xlo, xhi, ylo, yhi, zlo, zhi for each of the configura-
# tions. 
# DISCLAIMER: I dont know why in the past I didnt decide to simply get
# this from the all atom trajectory I am reading: much easier.
box_size = []
ofi1 = open("log.lammps", 'r')
for it_1 in range(0, number_of_configurations):
    dump = ofi1.readline()
    #dump = dump[0:32]
    for e,it_2 in zip(dump.split('\t'), range(7)):
        box_size.append(float(e))
box_size = np.array(box_size,float)
box_size = box_size.reshape(number_of_configurations,7)

# The second item will be achieved using the pos section of the first ever 
# configuration used to do the "NPT simulation". I will simply count how 
# many atoms composed a molecule having the same molecule ID and then figure 
# out whether it molecule is a MeOH (type 1) or EtOH (type 2). Note that se-
# veral other criteria could be chosen.

# Input here the number of atoms in the simulation domain.
NA = 6258

A = []
ofi1 = open("positions", 'r')
for it_1 in range(0, NA):
    dump = ofi1.readline()
    #dump = dump[0:32]
    for e,it_2 in zip(dump.split('\t'), range(7)):
        A.append(float(e))
A = np.array(A,float)
A = A.reshape(NA,7)

# This array here will take in the first column a molecule ID and in the
# second, the type of the molecule as defined above.
molecule_type = np.zeros((lines,2))

# This is the threshold for how many atoms in the count of the loop below
# I need to have to already make sure it's an EtOH and therefore be able
# to break the loop with certainty.
threshold = 7

# Creating the molecule ID column
for it_1 in range (0, lines):
    molecule_type[it_1,0] = it_1 + 1

# I will put the second column as if all were EtOH molecules and 
# then simply let the loop below change it accordingly for when
# it's MeOH.
molecule_type[:,1] = 2

for it_1 in range (0, lines):
    given_molecule = molecule_type[it_1,0]
    count_atoms = 0
    for it_2 in range (0, NA):
        if A[it_2,1] == given_molecule:
            count_atoms = count_atoms + 1
        if count_atoms == threshold:
            break
    if count_atoms <= (threshold - 1):
        molecule_type[it_1,1] = 1
    
#------------------------------------------------------
# Finally, let's now export the trajectory for the com of the molecules.
# Note that I am not going to wrap things back in the box (the coordina-
# tes of the com are may be unwrapped). This should be okay: LAMMPS will 
# understhand it and map it back to the simulation domain when reading
# this trajectory in a rerun command.
# One thing that is important to note here is that the com/chunk command
# in LAMMPS will assign a chunk ID equal to the molecule ID within the way 
# the script lammps_rdf-com1.in was specifically written. This is key for
# the coherence of everything that was done here.
ofi = open("com_traj.out", 'w')   
for it_1 in range(0, number_of_configurations):
    # ----------------------- HEADER ---------------------
    ofi.write('ITEM: TIMESTEP')
    ofi.write('\n')
    # We dont care about the specifics of this.
    ofi.write(str(int(it_1)))
    ofi.write('\n')
    ofi.write('ITEM: NUMBER OF ATOMS')
    ofi.write('\n')
    ofi.write(str(int(lines)))
    ofi.write('\n') 
    ofi.write('ITEM: BOX BOUNDS pp pp pp') 
    ofi.write('\n') 
    ofi.write(str(box_size[it_1,0]))
    ofi.write(' ') 
    ofi.write(str(box_size[it_1,1]))
    ofi.write('\n') 
    ofi.write(str(box_size[it_1,2]))
    ofi.write(' ') 
    ofi.write(str(box_size[it_1,3]))
    ofi.write('\n') 
    ofi.write(str(box_size[it_1,4]))
    ofi.write(' ') 
    ofi.write(str(box_size[it_1,5]))
    ofi.write('\n') 
    ofi.write('ITEM: ATOMS id type xu yu zu')
    ofi.write('\n') 
    # -------------------------------------------------
    for it_2 in range(0, lines):
        molecule_ID = com_traj[int(it_1*(lines) + it_2), 0]
        ofi.write(str(int(molecule_ID)))
        ofi.write(' ')
        ofi.write(str(int(molecule_type[int(molecule_ID - 1), 1])))
        for it_3 in range (1, 4):
            ofi.write(' ')
            ofi.write(str("{0:.6f}".format(com_traj[int(it_1*(lines) + it_2), it_3])))
        ofi.write('\n')
