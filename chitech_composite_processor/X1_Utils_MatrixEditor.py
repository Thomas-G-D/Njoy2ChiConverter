# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 10:00:10 2022

@author: tdeguire
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 12:52:49 2022

@author: tdeguire
"""

"""Combined Chi-tech data from multiple files to form composite cross sections"""

import numpy as np
import matplotlib.pyplot as plt
#def EditTransferMatrix(tape, data):
    
tape = 'tape41'    # input file
counter = 0 
counter2 = 0
test_mat = []
test_mat.append(np.zeros([172,172]))
test_mat.append(np.zeros([172,172]))
test_mat.append(np.zeros([172,172]))
test_mat.append(np.zeros([172,172]))
test_mat.append(np.zeros([172,172]))
test_mat.append(np.zeros([172,172]))
test_mat.append(np.zeros([172,172]))
test_mat.append(np.zeros([172,172]))

# %% this loop is to read elastic scattering matrix
for i in range (1):
    cf = open('tape41', 'r') #tape[i]
    structure = []  # initialize variable to determine matxsr structure
    nelas = []      # initialize variable to contain elastic scattering information
    for line in cf:
        # Get the group structure
        if line.find('8d     nelas') != -1:           # find where elastic scattering information starts
            structure_line = cf.readline().split()    # read the lines for the structure
            structure.append(structure_line)  
            while (structure_line[0] != '9d'):              # while still reading structure
                    structure_line = cf.readline().split()  # read structure line
                    structure.append(structure_line)        # add to structure variable
            nelas.append(structure_line)                    # append elastic scattering variable
            structure[0] = structure[0][2:] # remove random values
            structure = structure[:-1]      # remove last row (code reads one line too far)
            transfer_groups = [item for sublist in structure for item in sublist] # determine transfer group structure
            starting_groups = transfer_groups[-172:]    # determine which group is the starting group for each set of xs
            transfer_groups = transfer_groups[:-172]    # remove those lines from structure varible
            #transfer_groups = transfer_groups[:-172]   # do not need the last 172 variables: this is just energy groups
            
            nelas_line = cf.readline().split()  # continue reading transfer matrix
            nelas.append(nelas_line)  
            
            ### note this while loop may need appending; see what happens if only transfer matrix
            while (nelas_line[0] != '6d'):         # while still reading the transfer matrix
                nelas_line = cf.readline().split() # read line
                if len(nelas_line) < 6: # this if loop handles the case when negative values are present
                                        # if there are negative values, there are no spaces between numbers
                    new_line = ' '.join([str(item) for item in nelas_line]) # read line into a string variable
                    for char in range(1, len(new_line)):                    # loop over each character
                        # if the character is negative and not part of scientific notation, add a space
                        if new_line[char] == '-' and new_line[char-1] != 'E' and new_line[char-1] != ' ':
                            new_line = new_line[:char] +' '+new_line[char:]
                            char += 2 ### this line may be unnecessary
                    nelas.append(new_line.split())
                else: # if the line has the expected number of outputs, proceed normally
                    nelas.append(nelas_line)
            
            nelas = nelas[:-1]  # remove last line from list it reads an extra line
            nelas_list = [item for sublist in nelas for item in sublist] # put all values into a single list variable
            occurences = nelas_list.count('9d') # count how many times 9d occurs in the list variable
            # this loop removes all instances of 9d
            for c in range(occurences):
                nelas_list.remove('9d')
            grps_done = 0    # variable to index groups done, this is needed to follow the P band format
            nelas_list = [float(xs) for xs in nelas_list] # convert strings to floats
            transfer_groups = [int(grp) for grp in transfer_groups]     # convert str to int
            starting_groups = [int(grp) for grp in starting_groups]     # convert str to int
            
            # this loop builds the transfer matrix
            for dep_grps in range(len(transfer_groups)):                        # departing group
                for p in range(0,8):                                            # legendre coefficient number
                    for arr_grps in reversed(range(transfer_groups[dep_grps])): # arrival group
                        # the grp variable indexes the row in the transfer matrix
                        # the tape file is formatted to read on the order of the arrival group
                        # by descending transfer group, by increasing P
                        # e.g. group 172 -> 170 P_0, group 171 -> 170 P_0, 170 -> 170, P_0, 172->170 P_1 ...
                        # start index at the departing grooup
                        # arrival groups is based on the number of groups which scatter to that group, so add it
                        # subtract the total number of groups which scatter to that group
                        # add one to account for the fact that transfer_groups >= 1, so need to avoid negatives
                        grp = dep_grps+arr_grps-transfer_groups[dep_grps]+1
                        # here we determine the value to start at
                        # grps done indicates how many xs have already been indexed,
                        # multiply it by 8 to account for P_0 through P_7
                        # do p*transfer_groups to make sure the right P band is being selected
                        # add the total number of groups which scatter into the group of interest
                        # subtract the index of the arrival group value to look at each group in the for loop
                        # subtract 1 to account for indexing starting at 0
                        grp = dep_grps+arr_grps-transfer_groups[dep_grps]+1

                        row = starting_groups[dep_grps]-1
                        if abs(nelas_list[grps_done*8+p*transfer_groups[dep_grps]+arr_grps]) > 1e-18:
                            test_mat[p][row-arr_grps][dep_grps] = nelas_list[grps_done*8+p*transfer_groups[dep_grps]+arr_grps]
                        #test_mat[p][grp][dep_grps] = nelas_list[grps_done*8+p*transfer_groups[dep_grps]+transfer_groups[dep_grps]-arr_grps-1]

                grps_done += transfer_groups[dep_grps] # track how many groups have already been done
                
# %% this loop handles free gas scattering
        if line.find('8d     free') != -1: # find where free gas scattering starts
            structure = []  # initialize structure variable
            free_gas  = []  # initialize free gas variable
            structure_line = cf.readline().split() # read in structure of free gas matrix
            structure.append(structure_line)
            while (structure_line[0] != '9d'):              # while reading structure
                    structure_line = cf.readline().split()  # split each line
                    structure.append(structure_line)        # add it to the structure variable
            free_gas.append(structure_line)                 # add the last line of structure to the free gas variable
            structure[0] = structure[0][2:]                 # ignore the first two values
            structure = structure[:-1]                      # ignore the last row
            transfer_groups = [item for sublist in structure for item in sublist] # create transfer groups variable
            starting_groups = transfer_groups[-172:]    # determine which group is the starting group for each set of xs
            transfer_groups = transfer_groups[:-172]    # remove those lines from structure varible
            
            fg_line = cf.readline().split()             # read the lines containing free gas scattering info  
            free_gas.append(fg_line)                    # append these lines to free gas variable
            while (len(fg_line) != 0):# or (('6d') not in fg_line): ### need this to work with the or statement
                fg_line = cf.readline().split()  # keep reading lines

                if len(fg_line) < 6:    # this if loop handles the case when negative values are present
                                        # if there are negative values, there are no spaces between numbers
                    new_line = ' '.join([str(item) for item in fg_line]) # read line into a string variable
                    for char in range(1, len(new_line)):                 # loop over each character
                        # if the character is negative and not part of scientific notation, add a space
                        if new_line[char] == '-' and new_line[char-1] != 'E' and new_line[char-1] != ' ':
                            new_line = new_line[:char] +' '+new_line[char:]
                            char += 2 ### this line may be unnecessary
                    free_gas.append(new_line.split())
                else: # if the line has the expected number of outputs, proceed normally
                    free_gas.append(fg_line)
            
            free_gas = free_gas[:-1] # eliminate the last entry, while loop increments one past desired position
            nelas_list = [item for sublist in free_gas for item in sublist] # add all values into one list
            occurences = nelas_list.count('9d') # count all occurences of the string 9d
            for c in range(occurences):         # remove each occurence 
                nelas_list.remove('9d')

            nelas_list      = [float(xs) for xs in nelas_list]          # convert str to float
            transfer_groups = [int(grp) for grp in transfer_groups]     # convert str to int
            starting_groups = [int(grp) for grp in starting_groups]     # convert str to int
            position = 0    # count how many groups have already been done
            
            # this loop adds free gas scattering to the transfer matrix
            for dep_grps in range(len(transfer_groups)):
                for p in range(0,8):
                    for arr_grps in reversed(range(transfer_groups[dep_grps])):
                        # the grp variable indexes the row in the transfer matrix
                        # the tape file is formatted to read on the order of the arrival group
                        # by descending transfer group, by increasing P
                        # e.g. group 172 -> 170 P_0, group 171 -> 170 P_0, 170 -> 170, P_0, 172->170 P_1 ...
                        # start index at the departing grooup
                        # arrival groups is based on the number of groups which scatter to that group, so add it
                        # subtract the total number of groups which scatter to that group
                        # add one to account for the fact that transfer_groups >= 1, so need to avoid negatives
                        grp = dep_grps+arr_grps-transfer_groups[dep_grps]+1

                        row = starting_groups[dep_grps]-1
                        if abs(nelas_list[position*8+p*transfer_groups[dep_grps]+arr_grps]) > 1e-18:
                            test_mat[p][row-arr_grps][dep_grps] = nelas_list[position*8+p*transfer_groups[dep_grps]+arr_grps]


                position += transfer_groups[dep_grps]



M_sig_gp_to_g = np.zeros([172,172])
M_sig_gp_to_g[:,:] = test_mat[0][:,:]
  
#======================================= Solve for psi
S = M_sig_gp_to_g.transpose()
#T = np.diag(v_sig_t)
plt.matshow(np.log(S))
plt.title('Matxsr elastic scattering')
