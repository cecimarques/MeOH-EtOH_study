# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 13:57:38 2024

@author: u2371853
"""

# This code computes the function C(t) as defined by the paper by Richard and
# Paola (https://doi.org/10.1063/1.4922445, equation 9).
# The choice here is to assess the *continuous* H bond life time.

# ---------------------------------------------------------------------------
import numpy as np

# Input here the number of configurations used and the time they span in ps.
number_of_configurations = 2000
production_length = 1000 # ps

# Input here the number of atoms in your simulation domain
NA = 5982

# These values will receive the atom types of the O and H atoms of the MeOH
# and EtOH  molecules in the given system. Note that "atom type" is a definition
# within how simulations are carried out in LAMMPS and that these atom types sh-
# ould not be shared by any other atom for this code to work properly.
O_MeOH = 4
H_MeOH = 2
H_EtOH = 5
O_EtOH = 10

traj = []
box_size = []

# This bit reads the trajectory file output by LAMMPS during the production run.
# For each configuration, LAMMPS atom ID, LAMMPS atom type and x, y, z coordina-
# tes (unwrapped) were saved. The dimensions of the simulation domain, which he-
# re is cubic, at each configuration is also stored in the trajectory.
ofi = open("dump.pos", 'r')
for it_1 in range(0, number_of_configurations):
    ofi.readline()
    ofi.readline()
    ofi.readline()
    ofi.readline()
    ofi.readline()
    # -------------------------------------------------
    # Reading the dimensions of the simulation domain. There are 3 lines, each
    # containing the starting and ending coordinate of the simulation domain in
    # the given direction (i.e. x, y and z).
    for it_2 in range (0,3):
        dump = ofi.readline()
        for e,it_3 in zip(dump.split(' '), range(2)):
            box_size.append(float(e))
    # -------------------------------------------------
    dump = ofi.readline()
    # Reading information about the atoms in the simulation domain. Information
    # for each atom lies in a line and there are 5 columns containing, sequen-
    # tially, the LAMMPS atom ID, LAMMPS atom type and x, y, z coordinates (un-
    # wrapped).
    # Since I only need information about the hydrogens and oxygens to calculate
    # what this code is meant to calculate, I will not store information about
    # other atoms.
    for it_2 in range(0, NA):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_3 in zip(dump.split(' '), range(5)):
            if it_3 == 1:
                if (e != str(O_MeOH)) & (e != str(H_MeOH)) & (e != str(O_EtOH)) & (e != str(H_EtOH)):
                    if len(traj) != 0:
                        traj.pop()
                    break
            traj.append(float(e))
    # --------------------------------------------------
# Since I only stored O-Hs of MeOH and EtOH in the traj array,
# I need to redefine the value of NA, which now no longer is
# the number of atoms in the simulation domain but rather of
# O_MeOH, O_EtOH, H_MeOH and H_EtOH in the simulation domain. 
# The expression above for calculating it is solid: the divi-
# sion by 5 justifies itself on the fact that traj is still a 
# list (1D "python object") and each line, which contain data 
# for 1 atom, has 5 "occurances" in this list.
NA = int(len(traj)/(5*number_of_configurations))

traj = np.array(traj,float)
traj = traj.reshape(int(number_of_configurations*NA),5)
box_size = np.array(box_size,float)
box_size = box_size.reshape(int(number_of_configurations*3),2)

# ---------------------------------------------------------------------
# I will now read the Bonds section of the LAMMPS data file used in the si-
# mulation for this given system. It will be useful later when I need to as-
# ses the 3-atom angle value to determine whether or not a given O...H pair
# is hydrogen bonded.

# This is the number of bonds between all atoms composing the simulation 
# domain.
number_of_bonds = 5253

bonds = []
ofi1 = open("bonds", 'r')
for it_1 in range(0, number_of_bonds):
    dump = ofi1.readline()
    #dump = dump[0:32]
    for e,it_2 in zip(dump.split('\t'), range(4)):
        bonds.append(float(e))
bonds = np.array(bonds,float)
bonds = bonds.reshape(number_of_bonds,4)

# In order to work only with useful information (i.e. information necessa-
# ry for the code), I am going to delete all bonds that are not O-H in MeOH
# and EtOH molecules.
# Input below the bond types for the previously mentioned bonds (definition
# within LAMMPS).
bond_MeOH = 3
bond_EtOH = 8
for it_1 in reversed(range(0, number_of_bonds)):
    if (bonds[it_1,1] != bond_MeOH) & (bonds[it_1,1] != bond_EtOH):
        bonds = np.delete(bonds, it_1, 0)

# Let's now redefine the number of bonds after deleting the ir-
# relevant types.
number_of_bonds = len(bonds)

# ----------------------------------------------------------------
# For the methodology I am setting up to work, I will have to also order
# the traj in ascending order with respect to ID (basically, for my cal-
# culation of hydrogen bond lifetime, I am relying a 1D array for each con-
# figuration that will tell me if in a given configuration a given pair
# ij is H bonded or not. Having all configurations in the trajectory be
# ordered in ascending order with atom ID will ensure that a pair ij is
# found in the same iteration of loops to-be-further-on-defined an thus
# have each line of the previously mentioned array correspond to a same
# pair throughout all configurations considered in the H bond lifetime
# calculation.
temporary_array = np.zeros((1,5))
for it_1 in range (0, number_of_configurations):
    for it_2 in range(0, NA):
        for it_3 in range(it_2+1, NA):
            if traj[int(it_1*NA + it_2),0] > traj[int(it_1*NA + it_3),0]:
                temporary_array[0,:] = traj[int(it_1*NA + it_2),:]
                traj[int(it_1*NA + it_2),:] = traj[int(it_1*NA + it_3),:]
                traj[int(it_1*NA + it_3),:] = temporary_array[0,:]

# For the sake of organization, let's now map everything back to inside
# of the box (the trajectories were, afterall, saved in an unwrapped fa-
# shion).
for it_1 in range(0, int(number_of_configurations*NA)):
    xlo = box_size[int(int(it_1/NA)*3),0]
    xhi = box_size[int(int(it_1/NA)*3),1]
    lx = xhi - xlo
    ylo = box_size[int(int(it_1/NA)*3)+1,0]
    yhi = box_size[int(int(it_1/NA)*3)+1,1]
    ly = yhi - ylo
    zlo = box_size[int(int(it_1/NA)*3)+2,0]
    zhi = box_size[int(int(it_1/NA)*3)+2,1]
    lz = zhi - zlo
    k = 1 
    # k will only be 0 when I did not need to enter any of the if con-
    # ditions within the loop above, meaning that the atom is already
    # within the simulation domain.
    while k == 1:
        k = 0
        # ------------------------------
        # Checking for the x coordinate
        if traj[it_1,2] < xlo:
            traj[it_1,2] = traj[it_1,2] + lx
            k = 1
        if traj[it_1,2] > xhi:
            traj[it_1,2] = traj[it_1,2] - lx
            k = 1
        # ------------------------------
        # Checking for the y coordinate
        if traj[it_1,3] < ylo:
            traj[it_1,3] = traj[it_1,3] + ly
            k = 1
        if traj[it_1,3] > yhi:
            traj[it_1,3] = traj[it_1,3] - ly
            k = 1
        # ------------------------------
        # Checking for the z coordinate
        if traj[it_1,4] < zlo:
            traj[it_1,4] = traj[it_1,4] + lz
            k = 1
        if traj[it_1,4] > zhi:
            traj[it_1,4] = traj[it_1,4] - lz
            k = 1
        
# ----------------------------------------------------------------
import math as math

# I will also take the chance to check how many MeOH and EtOH molecules
# there are: these values will be useful for later. Since the number of
# molecules is constant throughout the simulation, I will simply check
# the first configuration. I will use the number of oxygens in the mole-
# cules to make the molecule count.
N_MeOH = 0
N_EtOH = 0
for it_1 in range(0,NA):
    if traj[it_1,1] == O_MeOH:
        N_MeOH = N_MeOH + 1
    if traj[it_1,1] == O_EtOH:
        N_EtOH = N_EtOH + 1

# -------------------------------------------------------------------------
# For reasons that will become understandable later, I will create a 1D 
# array for the given type of H bond I care about here, which will have 
# number of lines = number of corresponding ij pairs and 1 column. This 
# array will allow me to evaluate if a H-bond that existed in the previ-
# ous configuration broke apart in the next configuration or not.
# There is one python code for each possible H-bond that exists in the
# given system. In the nomenclature-scheme below, the first molecule al-
# ways refer to the H acceptor and the second, naturally, to the H donor.
MeOH_EtOH = np.zeros((int(N_MeOH*(N_EtOH)),1))
# Just "for the sake of it" these would be how the other arrays would lo-
# ok like if you were to assess other hydrogen bonds:
# MeOH_EtOH = np.zeros((int(N_MeOH*N_EtOH),1))
# EtOH_MeOH = np.zeros((int(N_EtOH*N_MeOH),1))
# EtOH_EtOH = np.zeros((int(N_EtOH*(N_EtOH-1)),1))

# Input in "sample_length" the time (in ps) the C(t) function will span.
# Input in "bin_size" the time space (in ps) between two configurations.
# This variable is called "bin size" because it marks values of times for
# which there will be values in the C(t) function.
# The number of bins in the function (i.e. number of values of t for which
# there will be data) follow from the previously made choices. 
# Note that the code works with equally spaced configurations.
sample_length = 100 # ps 
bin_size = 0.500 # ps
number_of_bins = int(sample_length/bin_size)  

# Definining the function C(t) for *one* subtrajectory
C_ME = np.zeros((number_of_bins,2))
for it_1 in range (0, number_of_bins):
    C_ME[it_1,0] = bin_size*it_1

# I chose (arbitrarily) that I will have 100 samples to build the Cx(t)
# function. Since all of them must praise for a simulation 100 ps long
# as defined earlier and my trajectory lasts for a total of 1000 ps, I 
# will need to create a new subjtractory (built based on a new time re-
# ference) every amount of time sample_step in case I want these samples 
# to be equally spaced in the sequence *AND* have the last subtrajectory
# last up the last configuration saved in the trajectory.
samples = 100  
sample_step = (production_length - sample_length)/samples  # ps
# Given that each configuration is 0.5 ps (which is 500 fs) apart from 
# one another, this amount of time corresponds to the following delta
# in terms of configuration
sample_step = sample_step/0.500  # configurations

# Definining the function C(t) for all subtrajectories averaged. This is
# end goal!!!!
average_CME = np.zeros((len(C_ME),2))
average_CME[:,0] = C_ME[:,0]

for a in range(0, samples):
    # For every subjtraectory I will have to reset the function Cx(t).
    C_ME[:,1] = 0
    # I also need to reset the array that counts if a given O...H pair
    # of MeOH molecules are H bonded or not.
    # Note that this array is ONLY reset here, when starting the asses-
    # ment of a new subtrajectory, AND NOT by the end of reading a con-
    # figuration. As it will be evident later, this allows keeping ties 
    # with pairs that were hydrogen bonded in the configuration previo-
    # us to the "current" one assessed in the loop it_1 below.
    MeOH_EtOH[:,0] = 0
    # ----------------------------------------------------------------
    # These here mark the starting and ending point of this subtrajec-
    # tory.
    # Since the configurations are 500 fs spaced and I want each sample
    # to go over a length of 0.1 ns, I should read up to 200 configurati-
    # ons after tstart.
    tstart = int(a*sample_step)
    tend = int(tstart + 200)
    # ----------------------------------------------------------------
    for it_1 in range(tstart, tend):
        A = np.zeros((NA,5))
        # Here I am adding the position of all atoms (which in this case
        # are only Os and Hs since I deleted the rest) of this given con-
        # figuration in the A array.
        A[:,:] = traj[int(it_1*NA):int((it_1+1)*NA),:]
        
        # I will also need later the size of the simulation domain.
        xlo = box_size[int(it_1*3),0]
        xhi = box_size[int(it_1*3),1]
        lx = xhi - xlo
        ylo = box_size[int(it_1*3)+1,0]
        yhi = box_size[int(it_1*3)+1,1]
        ly = yhi - ylo
        zlo = box_size[int(it_1*3)+2,0]
        zhi = box_size[int(it_1*3)+2,1]
        lz =  zhi - zlo
        
        # So, here, I will be checking for a possible H(MeOH)...O(MeOH)
        # hydrogen bond.
        # I will consider all possible pairs and identify whether 
        # or not a H bond exists by using the angle O-H...O and the 
        # distance H....O with the thresholds mentioned in the manuscript.
        # -------------------------------------------------------
        line_ME = 0
        for it_2 in range (0, NA):
            if A[it_2,1] == H_MeOH:
                ID_j = A[it_2,0]
                pos_j = np.array([[A[it_2,2], A[it_2,3], A[it_2,4]]])
                
                # -------------------------------------------------------
                # Let's find the O chemically bonded to this hydrogen.
                for it_3 in range(0, len(bonds)):
                    if bonds[it_3,2] == ID_j:
                        ID_i = bonds[it_3,3]
                        break
                    if bonds[it_3,3] == ID_j:
                        ID_i = bonds[it_3,2]
                        break
                # Now I gotta find where this given -O is in the A array
                # and get its position.
                for it_3 in range (0,NA):
                    if A[it_3,0] == ID_i:
                        pos_i = np.array([[A[it_3,2], A[it_3,3], A[it_3,4]]])
                        break
                # I will put atomi sitting the closest to atomj, in case it
                # isnt (due to the periodic boundary conditions).
                dist_x = pos_j[0,0] - pos_i[0,0]
                dist_y = pos_j[0,1] - pos_i[0,1]
                dist_z = pos_j[0,2] - pos_i[0,2]

                if (abs(dist_x) > lx/2) & (dist_x > 0):
                    pos_i[0,0] = pos_i[0,0] + lx
                if (abs(dist_x) > lx/2) & (dist_x < 0):
                    pos_i[0,0] = pos_i[0,0] - lx
                # --------------------------------------
                if (abs(dist_y) > ly/2) & (dist_y > 0):
                    pos_i[0,1] = pos_i[0,1] + ly
                if (abs(dist_y) > ly/2) & (dist_y < 0):
                    pos_i[0,1] = pos_i[0,1] - ly
                # --------------------------------------
                if (abs(dist_z) > lz/2) & (dist_z > 0):
                    pos_i[0,2] = pos_i[0,2] + lz
                if (abs(dist_z) > lz/2) & (dist_z < 0):
                    pos_i[0,2] = pos_i[0,2] - lz
                    
                # -------------------------------------------------------
                # Now I simply need to loop over all the O(EtOH) toand see 
                # if any of them is hydrogen-bonded to the given H(MeOH).
                for it_3 in range (0, NA):
                    if (A[it_3,1] == O_EtOH):
                        ID_k = A[it_3,0]
                        pos_k = np.array([[A[it_3,2], A[it_3,3], A[it_3,4]]])
                        
                        # Let's now evaluate whether or not I need to trans-
                        # late atom k with respect to j before checking for
                        # distances and angles.
                        dist_x = pos_j[0,0] - pos_k[0,0]
                        dist_y = pos_j[0,1] - pos_k[0,1]
                        dist_z = pos_j[0,2] - pos_k[0,2]
                        
                        if (abs(dist_x) > lx/2) & (dist_x > 0):
                            pos_k[0,0] = pos_k[0,0] + lx
                        if (abs(dist_x) > lx/2) & (dist_x < 0):
                            pos_k[0,0] = pos_k[0,0] - lx
                        # --------------------------------------
                        if (abs(dist_y) > ly/2) & (dist_y > 0):
                            pos_k[0,1] = pos_k[0,1] + ly
                        if (abs(dist_y) > ly/2) & (dist_y < 0):
                            pos_k[0,1] = pos_k[0,1] - ly
                        # --------------------------------------
                        if (abs(dist_z) > lz/2) & (dist_z > 0):
                            pos_k[0,2] = pos_k[0,2] + lz
                        if (abs(dist_z) > lz/2) & (dist_z < 0):
                            pos_k[0,2] = pos_k[0,2] - lz
                            
                        # Good, now let's check if there is a hydrogen
                        # bond O(MeOH)-H(MeOH)...O(MEtOH) between these a-
                        # toms.
                        vector_u = np.array([[pos_j[0,0] - pos_i[0,0], pos_j[0,1] - pos_i[0,1], pos_j[0,2] - pos_i[0,2]]])
                        vector_v = np.array([[pos_j[0,0] - pos_k[0,0], pos_j[0,1] - pos_k[0,1], pos_j[0,2] - pos_k[0,2]]])
                        mod_u = (vector_u[0,0]**2 + vector_u[0,1]**2 + vector_u[0,2]**2)**(1/2)
                        mod_v = (vector_v[0,0]**2 + vector_v[0,1]**2 + vector_v[0,2]**2)**(1/2)
                        
                        distance = mod_v
                        cos_value = (vector_u[0,0]*vector_v[0,0] + vector_u[0,1]*vector_v[0,1] + vector_u[0,2]*vector_v[0,2])/(mod_u*mod_v)
                        # --------------------------------------------
                        # Let's correct for possible float errors just in
                        # case
                        if cos_value > 1:
                            cos_value = 1
                        if cos_value < -1:
                            cos_value = -1
                        # --------------------------------------------
                        angle_value = math.acos(cos_value)*180/math.pi
                        # If I meet the condition below, I will declare
                        # that the hydrogen bond exists.
                        if (130 < angle_value <= 180) & (distance < 3):
                            # If this is a hydrogen bond and I am in the
                            # first configuration, I will declare the pa-
                            # ir is H bonded: If it is not, I dont need to
                            # do anything as the slot already has a 0 by
                            # default.
                            if it_1 == tstart: 
                                MeOH_EtOH[line_ME,0] = 1
                        # According to what I read, I will enter here
                        # everytime I dont enter the first if (i.e., all
                        # other cases).
                        else:
                            # If I am not in the first configuration, I
                            # will only need to alter the value of the
                            # array on this given line IF there was a H
                            # bond between them previously AND this bond
                            # just broke, which is when the condition of
                            # the else statement + of the if statement be-
                            # low is met (and is the only case in which it
                            # is met).
                            # If the bond did not break, the value at the
                            # slot will still be 1 as I am not changing
                            # it compared to the configuration before. 
                            # If a H bond formed, I do not care since I
                            # am calculating the continuous H bond life-
                            # time , so I should leave it as 0 as it al-
                            # ready should be from the previous iteration.
                            # Finally, the last "possible" possibility is
                            # if there was not a H bond before. In this ca-
                            # se the slot is already is 0 from the configu-
                            # ration before (or since the very definition
                            # of the MeOH_EtOH when we started this sam-
                            # ple in case there was never a H bond).
                            if (it_1 != tstart) & (MeOH_EtOH[line_ME,0] == 1):
                                MeOH_EtOH[line_ME,0] = 0
                        
                        line_ME = line_ME + 1
                        
        # ----------------------------------------------
        # Once I reach this point, I have finished the loop of it_2, which go-
        # es over a given configurations of the given subjtrajectory.
        # If I happen to be in the initial configuration of this subtrajectory,
        # I will create the denominator of the function Cx(t) which is based on  
        # the MeOH_EtOH array that I've just built (see equation 9 of Gowers and
        # Paola's paper or the equation shown in our manuscript).
        if it_1 == tstart:
            denominator = 0
            for it_2 in range (0, len(MeOH_EtOH)):
                denominator = denominator + MeOH_EtOH[it_2,0]
            # Furthermore, if I am in the initial configuration of this batch I 
            # also need to save the looks of this array MeOH_EtOH to later build 
            # my numerator of the C(t) function. I will save it under the name 
            # "factor".
            factor = np.zeros((len(MeOH_EtOH),1))
            factor[:,0] = MeOH_EtOH[:,0]
            # This array will now have 1 line and len(MeOH_EtOH) columns with 
            # the info on the columns ordered in the same order as it was in 
            # the lines.
            factor = np.transpose(factor)
            
        # ----------------------------------------------
        # Once I reach this point here, I can build the value of C(t) for the 
        # given time instant within this subjtraectory.
        MeOH_EtOHt = np.zeros((len(MeOH_EtOH),1))
        MeOH_EtOHt[:,0] = MeOH_EtOH[:,0]
        numerator = np.matmul(factor, MeOH_EtOHt)    
        corresponding_bin = int(it_1 - tstart)
        C_ME[corresponding_bin,1] = float(np.sum(numerator))/denominator
    
    # --------------------------------------------------------------
    # WHen I get here, I will have finished building a function C(t) corres-
    # ponding to a given subjtraectory.
    # This will store (a word that here means "sum") the function C_ME for all 
    # subtrajectories in the variable average_CME. Later I will be averaging
    # this
    average_CME[:,1] = average_CME[:,1] + C_ME[:,1]

average_CME[:,1] = average_CME[:,1]/samples

# Output the function we want! :)
ofi = open("Cx.dat", 'w')   
for it_1 in range(0, len(average_CME)):
    ofi.write(str("{0:.5f}".format(average_CME[it_1,0])))
    ofi.write(' ')
    ofi.write(str("{0:.5f}".format(average_CME[it_1,1])))
    ofi.write('\n')
    # --------------------------------------------------------
