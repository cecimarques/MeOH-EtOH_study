# -*- coding: utf-8 -*-
"""
@author: Cecilia Alvares
"""

# These code aims to create a frequency count of the the different possible
# "hydrogen bond scenarios" that happen for a given oxygen of methanol and
# ethanol molecules (O_MeOH or O_EtOH) in a given system. The different pos-
# sible "hydrogen bond scenarios"  are having an O_MeOH or O_EtOH hydrogen
# bonded to (a) no H, (b) 1 H_MeOH, (c) 1 H_EtOH, (d) 2 H_MeOH, (d) 2 H_EtOH, 
# (e) 1 H_MeOH and 1 H_EtOH, (f) none of the previous alternatives. Option 
# (f) is expected to have either 0 or very few counts stemming from acciden-
# tal events.
# There is one code for O(MeOH) and another code for O(EtOH), with this one
# being for O(EtOH).

# ----------------------------------------------------------------------
import numpy as np
import math as math

# Input here the number of configurations used.
number_of_configurations = 2000

# # Input here the number of atoms and number of bonds in your simulation do-
# main
NA = 6258
number_of_bonds = 5529

# These values will receive the atom types of the O and H atoms of the MeOH
# and EtOH  molecules in the given system. Note that "atom type" is a definition
# within how simulations are carried out in LAMMPS and that these atom types sh-
# ould not be shared by any other atom for this code to work properly.
H_MeOH = 2
O_MeOH = 4
O_EtOH = 10
H_EtOH = 5

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
                if (e != str(H_EtOH)) & (e != str(O_EtOH)) & (e != str(O_MeOH)) & (e != str(H_MeOH)):
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
    if (bonds[it_1,1] != bond_EtOH) & (bonds[it_1,1] != bond_MeOH):
        bonds = np.delete(bonds, it_1, 0)

# Let's now redefine the number of bonds after deleting the ir-
# relevant types.
number_of_bonds = len(bonds)

# ------------------------------------------------------------------------
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
# Let's make some definitions for the hydrogen bond coordination number histo-
# grams I am trying to create. 

# This array is symbolic. It has 7 indexes and each of them represent the
# (a)-(g) "hydrogen bond scenarios" that were initially mentioned (and are 
# also presented below to refresh your memory).
# Scenarios: (a) none, (b) 1 H_MeOH, (c) 1 H_EtOH, (d) 2 H_MeOH, (e) 2 H_EtOH, 
# (f) 1 H_MeOH and 1 H_EtOH, (g) none of the previous alternatives.
x_axis = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]])
# This will count how many Os falling within each of the previous cases there 
# are.
O_coord = np.zeros((1, 7))

# These are the angle & bond criteria for deciding whether or not an O...H
# pair is hydrogen bonded.
# You may change accordingly if you wish.
angle_min = 130
dist_max = 3

for it_1 in range(0, number_of_configurations):
    for it_2 in range(0, NA):
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
        
        # ----------------------------------------------------------
        # If I have found an atom O_EtOH when running the loop of it_2, I
        # will then search its coordination number for possible O(EtOH)...
        # H(MeOH) and O(EtOH)...H(EtOH) hydrogen bonds.
        # I will be referring to this O_EtOH as "atom i" in the variables
        # and further comments, while the H_EtOH and H_MeOH are referred
        # to as "atom j".
        if A[it_2,1] == O_EtOH:
            ID_i = A[it_2,0]
            pos_i = np.array([[A[it_2,2], A[it_2,3], A[it_2,4]]])
            # This variable will receive the number of H_MeOH and H_EtOH I 
            # find hydrogen bonded to the given O_EtOH at the given confi-
            # guration.
            counter_Hm = 0
            counter_He = 0
            for it_3 in range (0, NA):
                if A[it_3,1] == H_EtOH:
                    ID_j = A[it_3,0]
                    pos_j = np.array([[A[it_3,2], A[it_3,3], A[it_3,4]]])
                    # -----------------------------------------------------
                    # I will put atomi sitting the closest to atomj, in case it
                    # isnt and recompute the distances later.
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
                    
                    dist_total = ((pos_j[0,0] - pos_i[0,0])**2 + (pos_j[0,1] - pos_i[0,1])**2 + (pos_j[0,2] - pos_i[0,2])**2)**(1/2)
                    # I only need to go further if this pair is hydrogen bonded.
                    # Note that I will also meet this condition for the H_EtOH
                    # who is chemically bonded to this given O_EtOH. I will
                    # take acre of this later.
                    if dist_total >= dist_max:
                        continue
                    # ----------------------------------------------------
                    # If this is a potential H bond as the condition above 
                    # is NOT met, I will then check the O(EtOH)...H(EtOH)-O
                    # angle to confirm whether or not this is a hydrogen
                    # bond and thus if I indeed should be counting this H
                    # in the coordination number of this given O(EtOH) in
                    # this configuration.
                    # I will be refering to this oxygen bonded to atom H(Et-
                    # OH) as atom k.
                    for it_4 in range (0, number_of_bonds):
                        if bonds[it_4,2] == ID_j:
                            ID_k = bonds[it_4,3]
                            break
                        if bonds[it_4,3] == ID_j:
                            ID_k = bonds[it_4,2]
                            break
                    # ----------------------------------------------------
                    # Now, if this is the H_EtOH chemically bonded to the gi-
                    # ven O_EtOH, it is not the one I am looking for and thus
                    # I will stop this iteration and go to the next iteration
                    # in loop of it_3.
                    if ID_k == ID_i:
                        continue
                    # ----------------------------------------------------
                    for it_4 in range (0, NA):
                        if A[it_4,0] == ID_k:
                            pos_k = np.array([[A[it_4,2], A[it_4,3], A[it_4,4]]])
                            break
        
                    # Let's put j and k in the minimum distance in case they
                    # are not.
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
                    # -----------------------------------------------------
                    
                    vector_u = np.array([[pos_j[0,0] - pos_i[0,0], pos_j[0,1] - pos_i[0,1], pos_j[0,2] - pos_i[0,2]]])
                    vector_v = np.array([[pos_j[0,0] - pos_k[0,0], pos_j[0,1] - pos_k[0,1], pos_j[0,2] - pos_k[0,2]]])
                    mod_u = (vector_u[0,0]**2 + vector_u[0,1]**2 + vector_u[0,2]**2)**(1/2)
                    mod_v = (vector_v[0,0]**2 + vector_v[0,1]**2 + vector_v[0,2]**2)**(1/2)
                    
                    cos_value = (vector_u[0,0]*vector_v[0,0] + vector_u[0,1]*vector_v[0,1] + vector_u[0,2]*vector_v[0,2])/(mod_u*mod_v)
                    # --------------------------------------------
                    # Let's correct for possible float errors just in
                    # case.
                    # PS: I am going to set them smaller than 0.9999 simply
                    # for the sake of later being able to accumulate them
                    # in the last biny value (this can make more sense on-
                    # ce you get the bottom of the code).
                    if cos_value > 1:
                        cos_value = 0.9999
                    if cos_value < -1:
                        cos_value = -0.9999
                    # --------------------------------------------
                    angle_value = math.acos(cos_value)*180/math.pi
                    # Note that I am no longer checking the criterion for
                    # the O_EtOH....H_EtOH distance, since this has alre-
                    # ady been tacked in the if condition above.
                    if (angle_min < angle_value <= 180):
                        counter_He = counter_He + 1
                
                # Now, I may also have a O(EtOH)....H(MeOH) bond with this
                # given O(EtOH). So let's check for it. This calculation
                # follows the same frame as for H_EtOH above.
                if A[it_3,1] == H_MeOH:
                    ID_j = A[it_3,0]
                    pos_j = np.array([[A[it_3,2], A[it_3,3], A[it_3,4]]])
                    # -----------------------------------------------------
                    # I will put atomi sitting the closest to atomj, in case it
                    # isnt and recompute the distances later.
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
                    
                    dist_total = ((pos_j[0,0] - pos_i[0,0])**2 + (pos_j[0,1] - pos_i[0,1])**2 + (pos_j[0,2] - pos_i[0,2])**2)**(1/2)
                    if dist_total >= dist_max:
                        continue
                    # ----------------------------------------------------
                    for it_4 in range (0, number_of_bonds):
                        if bonds[it_4,2] == ID_j:
                            ID_k = bonds[it_4,3]
                            break
                        if bonds[it_4,3] == ID_j:
                            ID_k = bonds[it_4,2]
                            break

                    for it_4 in range (0, NA):
                        if A[it_4,0] == ID_k:
                            pos_k = np.array([[A[it_4,2], A[it_4,3], A[it_4,4]]])
                            break
        
                    # Let's put j and k in the minimum distance in case they
                    # are not.
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
                    # -----------------------------------------------------
                    
                    vector_u = np.array([[pos_j[0,0] - pos_i[0,0], pos_j[0,1] - pos_i[0,1], pos_j[0,2] - pos_i[0,2]]])
                    vector_v = np.array([[pos_j[0,0] - pos_k[0,0], pos_j[0,1] - pos_k[0,1], pos_j[0,2] - pos_k[0,2]]])
                    mod_u = (vector_u[0,0]**2 + vector_u[0,1]**2 + vector_u[0,2]**2)**(1/2)
                    mod_v = (vector_v[0,0]**2 + vector_v[0,1]**2 + vector_v[0,2]**2)**(1/2)
                    
                    cos_value = (vector_u[0,0]*vector_v[0,0] + vector_u[0,1]*vector_v[0,1] + vector_u[0,2]*vector_v[0,2])/(mod_u*mod_v)
                    # --------------------------------------------
                    if cos_value > 1:
                        cos_value = 0.9999
                    if cos_value < -1:
                        cos_value = -0.9999
                    # --------------------------------------------
                    angle_value = math.acos(cos_value)*180/math.pi
                    if (angle_min < angle_value <= 180):
                        counter_Hm = counter_Hm + 1
  
            # Once I get here I will have finished iterating over the loop
            # of it_3 and thus I will have counted how many Hs are hydro-
            # gen bonded to the given O_EtOH in this given configuration. 
            # It is then time to save this information accordingly.
            # First, lets discover which "hydrogen bond scenario" this is.
            # Reminder: (a) none, (b) 1 H_MeOH, (c) 1 H_EtOH, (d) 2 H_MeOH, 
            # (e) 2 H_EtOH, (f) 1 H_MeOH and 1 H_EtOH, (g) none of the pre-
            # vious alternatives.
            # PS: the variable "flag" will help me realize whether or not a-
            # ny of the six more conventional cases happened or not so that
            # I can count the anomalous event (whatever the specifics of it).
            flag = 0
            if (counter_He == 0) & (counter_Hm == 0):
                O_coord[0,0] = O_coord[0,0] + 1
                flag = 1
            # ------------------------------------------------------------
            if (counter_He == 0) & (counter_Hm == 1):
                O_coord[0,1] = O_coord[0,1] + 1
                flag = 1
            # ------------------------------------------------------------
            if (counter_He == 1) & (counter_Hm == 0):
                O_coord[0,2] = O_coord[0,2] + 1
                flag = 1
            # ------------------------------------------------------------
            if (counter_He == 0) & (counter_Hm == 2):
                O_coord[0,3] = O_coord[0,3] + 1
                flag = 1
            # ------------------------------------------------------------
            if (counter_He == 2) & (counter_Hm == 0):
                O_coord[0,4] = O_coord[0,4] + 1
                flag = 1
            # ------------------------------------------------------------
            if (counter_He == 1) & (counter_Hm == 1):
                O_coord[0,5] = O_coord[0,5] + 1
                flag = 1
            # ------------------------------------------------------------
            if flag == 0:
                O_coord[0,6] = O_coord[0,6] + 1
            # ------------------------------------------------------------

# -------------------------------------------------------------------------
# Creating the output file
output = np.concatenate((x_axis, O_coord), axis = 0)
output = np.transpose(output)

# Ideally, a normalization should be made here in order to be able to com-
# pare the different systems, since these contemplate a different number
# of O_MeOH and O_EtOH in the simulation domain. The normalization to be ma-
# de here is a division over (number_of_configurations)*O_EtOH: this is a-
# fterall the total amount of counts in the array "output". Let's find 
# out how many O_EtOH we have (note that I could just have input this, but
# I decided to calculate here for practicity). 
# For the sake of that, I am going to find how many O_EtOH I have in the
# last configuration I have previously ran in the loop of it_1.
NO = 0
for it_1 in range(0, NA):
    if A[it_1,1] == O_EtOH:
        NO = NO + 1
        
output[:,1] = output[:,1]/(number_of_configurations*NO)
# ---------------------------------------------------------------------

ofi = open("output.dat", 'w')   
for it_1 in range(0, len(output)):
    ofi.write(str("{0:.2f}".format(output[it_1,0])))
    ofi.write(' ')
    ofi.write(str("{0:.2f}".format(output[it_1,1])))
    ofi.write('\n')
