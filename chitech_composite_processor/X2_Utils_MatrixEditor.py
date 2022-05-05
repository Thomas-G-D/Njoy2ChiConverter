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

def EditMatrix(file, key_phrase, scat_type, test_mat):
    cf = open(file, 'r') #tape[i]
    structure = []  # initialize variable to determine matxsr structure
    xs_matrix = []      # initialize variable to contain elastic scattering information
    for line in cf:
        # Get the group structure
        if line.find(key_phrase) != -1:           # find where elastic scattering information starts
            structure_line = cf.readline().split()    # read the lines for the structure
            structure.append(structure_line)  
            while (structure_line[0] != '9d'):              # while still reading structure
                    structure_line = cf.readline().split()  # read structure line
                    structure.append(structure_line)        # add to structure variable
            xs_matrix.append(structure_line)                    # append elastic scattering variable
            structure[0] = structure[0][2:] # remove random values
            structure = structure[:-1]      # remove last row (code reads one line too far)
            transfer_groups = [item for sublist in structure for item in sublist] # determine transfer group structure
            starting_groups = transfer_groups[-172:]    # determine which group is the starting group for each set of xs
            transfer_groups = transfer_groups[:-172]    # remove those lines from structure varible
            #transfer_groups = transfer_groups[:-172]   # do not need the last 172 variables: this is just energy groups
            
            xs_line = cf.readline().split()  # continue reading transfer matrix
            xs_matrix.append(xs_line)  
            
            ### note this while loop may need appending; see what happens if only transfer matrix
            if scat_type == 'elastic':
                stop_line = xs_line[0]
                stop_cond = '6d'
            elif scat_type == 'free gas':
                stop_line = len(xs_line)
                stop_cond = 0
            while (stop_line != stop_cond):         # while still reading the transfer matrix
                xs_line = cf.readline().split() # read line
                if len(xs_line) < 6: # this if loop handles the case when negative values are present
                                        # if there are negative values, there are no spaces between numbers
                    new_line = ' '.join([str(item) for item in xs_line]) # read line into a string variable
                    for char in range(1, len(new_line)):                    # loop over each character
                        # if the character is negative and not part of scientific notation, add a space
                        if new_line[char] == '-' and new_line[char-1] != 'E' and new_line[char-1] != ' ':
                            new_line = new_line[:char] +' '+new_line[char:]
                            char += 2 ### this line may be unnecessary
                    xs_matrix.append(new_line.split())
                else: # if the line has the expected number of outputs, proceed normally
                    xs_matrix.append(xs_line)
                if scat_type == 'elastic':
                    stop_line = xs_line[0]
                    #stop_cond = '6d'
                elif scat_type == 'free gas':
                    stop_line = len(xs_line)
                    #stop_cond = 0
            
            xs_matrix = xs_matrix[:-1]  # remove last line from list it reads an extra line
            xs_list = [item for sublist in xs_matrix for item in sublist] # put all values into a single list variable
            occurences = xs_list.count('9d') # count how many times 9d occurs in the list variable
            # this loop removes all instances of 9d
            for c in range(occurences):
                xs_list.remove('9d')
            grps_done = 0    # variable to index groups done, this is needed to follow the P band format
            xs_list = [float(xs) for xs in xs_list] # convert strings to floats
            transfer_groups = [int(grp) for grp in transfer_groups]     # convert str to int
            starting_groups = [int(grp) for grp in starting_groups]     # convert str to int
            
            # this loop builds the transfer matrix
            for arr_grp in range(len(transfer_groups)):                        # departing group
                for p in range(0,8):                                            # legendre coefficient number
                    for dep_grp in reversed(range(transfer_groups[arr_grp])): # arrival group
                        # the grp variable indexes the row in the transfer matrix
                        # the tape file is formatted to read on the order of the arrival group
                        # by descending transfer group, by increasing P
                        # e.g. group 172 -> 170 P_0, group 171 -> 170 P_0, 170 -> 170, P_0, 172->170 P_1 ...
                        # start index at the departing grooup
                        # arrival groups is based on the number of groups which scatter to that group, so add it
                        # subtract the total number of groups which scatter to that group
                        # add one to account for the fact that transfer_groups >= 1, so need to avoid negatives
                        #grp = dep_grps+arr_grps-transfer_groups[dep_grps]+1
                        # here we determine the value to start at
                        # grps done indicates how many xs have already been indexed,
                        # multiply it by 8 to account for P_0 through P_7
                        # do p*transfer_groups to make sure the right P band is being selected
                        # add the total number of groups which scatter into the group of interest
                        # subtract the index of the arrival group value to look at each group in the for loop
                        # subtract 1 to account for indexing starting at 0
                        #grp = dep_grps+arr_grps-transfer_groups[dep_grps]+1

                        row = starting_groups[arr_grp]-1
                        if abs(xs_list[grps_done*8+p*transfer_groups[arr_grp]+dep_grp]) > 1e-18:
                            test_mat[p][row-dep_grp][arr_grp] = xs_list[grps_done*8+p*transfer_groups[arr_grp]+dep_grp]
                        #test_mat[p][grp][dep_grps] = nelas_list[grps_done*8+p*transfer_groups[dep_grps]+transfer_groups[dep_grps]-arr_grps-1]

                grps_done += transfer_groups[arr_grp] # track how many groups have already been done
                
    return test_mat
test_mat = EditMatrix('tape41', '8d     nelas', 'elastic', test_mat)
test_mat = EditMatrix('tape41', '8d     free', 'free gas', test_mat)




# M_sig_gp_to_g = np.zeros([172,172])
# M_sig_gp_to_g[:,:] = test_mat[0][:,:]
  
# #======================================= Solve for psi
# S = M_sig_gp_to_g.transpose()
# #T = np.diag(v_sig_t)
# plt.matshow(np.log(S))
# plt.title('Matxsr elastic scattering')
