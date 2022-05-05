# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 12:26:25 2022

@author: tdeguire
"""

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
    
# tape        = 'tape41'             # input file
# matxsr_data = {}                   # initialize dictionary to hold xs data
# #matxsr_data['transfer_mat'] = []   # initialize the transfer matrix entry
# test_mat = np.zeros([172,172])
# #n_groups = 172

def EditMatrix(file):
    # want this not hard coded
    matxsr_data = {}                   # initialize dictionary to hold xs data
    n_groups = 172
    transfer_mat       = [] # initialize variable to store xs info
    transfer_mat.append(np.zeros([n_groups,n_groups]))
    transfer_mat.append(np.zeros([n_groups,n_groups]))
    transfer_mat.append(np.zeros([n_groups,n_groups]))
    transfer_mat.append(np.zeros([n_groups,n_groups]))
    transfer_mat.append(np.zeros([n_groups,n_groups]))
    transfer_mat.append(np.zeros([n_groups,n_groups]))
    transfer_mat.append(np.zeros([n_groups,n_groups]))
    transfer_mat.append(np.zeros([n_groups,n_groups]))
    cf = open('tape41', 'r') #tape[i]
    structure = []   # initialize list to determine matxsr structure
    xs_matrix = []   # initialize list to contain elastic scattering information
    last_line = []   # initialize a list for containg lines at end of a specific xs
    last_line.append('start') # need this variable initialized to prevent error
    
    for line in cf:
        structure = []  # initialize variable to determine matxsr structure
        xs_matrix = []  # initialize variable to contain elastic scattering information
        
        #%% first, we read out the energy group structure
        if line.find('4d') != -1:            # find where energy groups start
            structure_line = line.split()    # split the line being read
            structure.append(structure_line) # add it to the structure
            
            # this loop reads until 5d appears, which occurs at the end of the energy loops
            while (structure_line[0] != '5d'):          # while still reading structure 
                structure_line = cf.readline().split()  # read structure line
                structure.append(structure_line)        # add to structure variable
            
            # now we need to format what we just read
            structure    = structure[:-1]   # remove the last line that contains 5d, one line too far read by while loop
            structure[0] = structure[0][1:] # remove the 4d string
            energy_bins  = [float(bins) for sublist in structure for bins in sublist] # determine transfer group str
            matxsr_data['energy_bins'] = energy_bins  # save energy bins to dictionary
            n_groups = len(energy_bins)-1
        #%% next assemble the xs values
        if line.find('6d') != -1 or last_line[0] == '6d':
            
            # this if loop initializes how the xs data will be reading depending
            # on how this portion of the code was entered
            if last_line[0] == '6d':
                structure_line = last_line             # take the last line read from the previous loop
                structure.append(structure_line)       # append it to the structure variable
                structure.append(line.split())         # then add the current line being read
            else:
                structure.append(line.split())         # add the line being read currently
                structure_line = cf.readline().split() # read the next line
                structure.append(structure_line)       # add it to the matrix
            
            # this loop reads until 7d appears, which occurs at the end of xs group structure
            while (structure_line[0] != '7d'):          # while still reading structure 
                structure_line = cf.readline().split()  # read structure line
                structure.append(structure_line)        # add to structure variable
            
            # now we need to format what we just read
            last_line      = structure[-1]       # save the last line; it is needed to the xs's
            structure      = structure[:-1]      # remove the last line that contains 7d
            structure[0]   = structure[0][1:]    # remove the 6d
            struc_info     = [info for sublist in structure for info in sublist] # read structure into one list
            num_xs         = int(len(struc_info)/3)            # determine number of xs contained in structure
            info_name      = struc_info[0:num_xs]              # determine the names of the xs
            info_grp_start = struc_info[num_xs:2*num_xs]       # determine the first energy group for the xs
            info_grp_start = [int(i) for i in info_grp_start]  # convert str to int
            info_grp_end   = struc_info[2*num_xs:3*num_xs]     # determine how many grps have values for that xs
            info_grp_end   = [int(i) for i in info_grp_end]    # cpmvert str to int
            
            # reset the structure variable; it will now be used to read the xs directly
            structure = []
            structure.append(last_line) # add the saved last line; first line of xs's
            
            # read until 8d which indicates a new xs matrix
            while (structure_line[0] != '8d'):          # while still reading structure 
                structure_line = cf.readline().split()  # read structure line
                structure.append(structure_line)        # add to structure variable
            
            # now we need to format what we just read
            last_line = structure[-1]       # save the last line; it is needed for the next xs
            structure = structure[:-1]      # remove the last line that contains 9d
            structure[0] = structure[0][1:] # remove the 8d
            xs_info = [info for sublist in structure for info in sublist] # read all xs's into one list
            
            # this loop assembles lists for each xs and adds them to the dictionary
            arr_position = 0 # initialize position in the xs array; needed to move from one xs to the next
            for single_xs in range(num_xs):  # for each xs
                xs_array = np.zeros(n_groups)     # initialize the array to be of length 172 for 172 energy groups
                for grp in range(info_grp_start[single_xs], info_grp_end[single_xs]+1): # for each group
                    index_pos          = grp+arr_position-info_grp_start[single_xs]  # index position in the array
                    grp_nmbr           = grp-1                                       # grp number
                    xs_array[grp_nmbr] = xs_info[index_pos]                          # add value to xs array
                arr_position += grp                            # increment starting position
                xs_array = [float(i) for i in xs_array]        # convert values from strings to floats
                matxsr_data[info_name[single_xs]] = xs_array   # add xs to matxsr dictionary
            continue
            
        #%% starting reading out transfer matrices
        if line.find('8d') != -1 or last_line[0] == '8d':
            structure = []  # initialize variable to determine matxsr structure
            xs_matrix = []  # initialize variable to contain elastic scattering information

            if last_line[0] == '8d':
                structure_line = last_line             # take the last line read from the previous loop
                structure.append(structure_line)       # append it to the structure variable
                structure.append(line.split())         # then add the current line being read
            else:
                structure.append(line.split())          ### this may be a mistake
                structure_line = cf.readline().split()  # append it to the structure variable
                structure.append(structure_line)        # add it to the matrix     
            
            info_name = structure[0][1]                 # xs type
            
            # this loop extracts the structure for the xs matrix
            while structure_line[0] != '9d':                # while still reading structure
                    structure_line = cf.readline().split()  # read structure line
                    structure.append(structure_line)        # add to structure variable
            
            last_line = structure[-1]      # save the last line, it is needed for building the xs matrix
            structure = structure[1:]       # remove the first line, it is unnecessary
            structure[0] = structure[0][2:] # remove unneeded values
            structure = structure[:-1]      # remove last row (code reads one line too far)
            struc_info = [item for sublist in structure for item in sublist] # determine transfer group structure
            info_grp_start = struc_info[-n_groups:] # determine which group is the starting group for each set of xs
            info_grp_end = struc_info[:n_groups]    # remove those lines from structure varible

            # this loop corrects the formatting of any negative values
            # if there are negative values, there are no spaces between numbers
            xs_info = last_line  # formerly structure_line
            if len(xs_info) < 6:
                    new_line = ' '.join([str(item) for item in xs_info])    # read line into a string variable
                    for char in range(1, len(new_line)):                    # loop over each character
                        # if the character is negative and not part of scientific notation, add a space
                        if new_line[char] == '-' and new_line[char-1] != 'E' and new_line[char-1] != ' ':
                            new_line = new_line[:char] +' '+new_line[char:]
                    xs_matrix.append(new_line.split())       # append xs matrix
            else:
                xs_matrix.append(xs_info)                    # append xs matrix
                
            # this loop reads all of the legendre coefficients into a variable
            # and stops at the next xs matrix
            while  len(xs_info) != 0 and xs_info[0] != '8d' and xs_info[0] != '6d':         # while still reading the transfer matrix
                xs_info = cf.readline().split() # read line
                if len(xs_info) < 6:    # this if loop handles the case when negative values are present
                                        # if there are negative values, there are no spaces between numbers
                    new_line = ' '.join([str(item) for item in xs_info]) # read line into a string variable
                    for char in range(1, len(new_line)):                    # loop over each character
                        # if the character is negative and not part of scientific notation, add a space
                        if new_line[char] == '-' and new_line[char-1] != 'E' and new_line[char-1] != ' ':
                            new_line = new_line[:char] +' '+new_line[char:]
                    xs_matrix.append(new_line.split())
                else: # if the line has the expected number of outputs, proceed normally
                    xs_matrix.append(xs_info)


            last_line = xs_matrix[-1]   # save the last line, this may be needed for another xs
            xs_matrix = xs_matrix[:-1]  # remove last line from list it reads an extra line
            xs_list = [item for sublist in xs_matrix for item in sublist] # put all values into a single list variable
            
            occurences = xs_list.count('9d') # count how many times 9d occurs in the list variable
            # this loop removes all instances of 9d
            for c in range(occurences):
                xs_list.remove('9d')
            
            
            xs_list        = [float(xs) for xs in xs_list]         # convert strings to floats
            info_grp_end   = [int(grp) for grp in info_grp_end]    # convert str to int
            info_grp_start = [int(grp) for grp in info_grp_start]  # convert str to int
           
            # feel like this can be done better
            xs_array       = [] # initialize variable to store xs info
            xs_array.append(np.zeros([n_groups,n_groups]))
            xs_array.append(np.zeros([n_groups,n_groups]))
            xs_array.append(np.zeros([n_groups,n_groups]))
            xs_array.append(np.zeros([n_groups,n_groups]))
            xs_array.append(np.zeros([n_groups,n_groups]))
            xs_array.append(np.zeros([n_groups,n_groups]))
            xs_array.append(np.zeros([n_groups,n_groups]))
            xs_array.append(np.zeros([n_groups,n_groups]))

            grps_done = 0    # variable to index groups done, this is needed to follow the P band format
            # this loop builds the transfer matrix
            for arr_grp in range(len(info_grp_end)):                        # departing group
                for p in range(0,8):                                        # legendre coefficient number
                    for dep_grp in reversed(range(info_grp_end[arr_grp])): # arrival group
                        # the grp variable indexes the row in the transfer matrix
                        # the tape file is formatted to read on the order of the arrival group
                        # by descending transfer group, by increasing P
                        # e.g. group 172 -> 170 P_0, group 171 -> 170 P_0, 170 -> 170, P_0, 172->170 P_1 ...
                        # start index at the departing grooup
                        # arrival groups is based on the number of groups which scatter to that group, so add it
                        # subtract the total number of groups which scatter to that group
                        # add one to account for the fact that info_grp_end >= 1, so need to avoid negatives
                        #grp = dep_grps+arr_grps-info_grp_end[dep_grps]+1
                        # here we determine the value to start at
                        # grps done indicates how many xs have already been indexed,
                        # multiply it by 8 to account for P_0 through P_7
                        # do p*info_grp_end to make sure the right P band is being selected
                        # add the total number of groups which scatter into the group of interest
                        # subtract the index of the arrival group value to look at each group in the for loop
                        # subtract 1 to account for indexing starting at 0
                        #grp = dep_grps+arr_grps-info_grp_end[dep_grps]+1

                        row = info_grp_start[arr_grp]-1
                        if abs(xs_list[grps_done*8+p*info_grp_end[arr_grp]+dep_grp]) > 1e-18:
                            xs_array[p][row-dep_grp][arr_grp] = xs_list[grps_done*8+p*info_grp_end[arr_grp]+dep_grp]

                grps_done += info_grp_end[arr_grp] # track how many groups have already been done
            matxsr_data[info_name+'_mat'] = xs_array
            
            #matxsr_data['transfer_mat'] += xs_array
            if info_name =='nx':
                trial = 1
            numbers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            numbers = str(numbers)
            loop_name = info_name + '_mat'
            if loop_name != 'free_mat':# and loop_name[2] not in numbers or loop_name == 'n01_mat':
                #print(loop_name)
                for p in range(8):
                    transfer_mat[p] += xs_array[p]
    matxsr_data['transfer_mat'] = transfer_mat
    #%% apply free gas treatment
    for p in range(8):
        for dep_grp in range(n_groups):
            for arr_grp in range(n_groups):
              if abs(matxsr_data['free_mat'][p][dep_grp][arr_grp]) > 1e-18:
                matxsr_data['transfer_mat'][p][dep_grp][arr_grp] = matxsr_data['free_mat'][p][dep_grp][arr_grp]
    return matxsr_data['transfer_mat']
        
        
        
        
        
        
        
        
        
     
    cf.close()
    return matxsr_data
# matxsr_data = EditMatrix('tape41')

# #matxsr_data['energy_bins'].append(0)
# matxsr_data['energy_bins'].reverse()
# n_bndrys, n_vals = [], []

# for grp in range(172):
#     lo_bound  = matxsr_data['energy_bins'][grp]*1.0e-6
#     hi_bound  = matxsr_data['energy_bins'][grp+1]*1.0e-6
#     bin_width = hi_bound-lo_bound
#     spectrum  = matxsr_data['nwt0'][grp]/ bin_width
#     n_bndrys += [lo_bound, hi_bound]
#     n_vals   +=  [spectrum, spectrum]
  
# n_bndrys = np.array(n_bndrys)
# n_vals = np.array(n_vals)

# # plt.figure()
# # plt.semilogy(n_bndrys, n_vals)
# # plt.figure()
# # plt.loglog(n_bndrys, n_vals)
# # need to do phi = E/de
# # and then need a way to plot bins


# M_sig_gp_to_g = np.zeros([172,172])
# M_sig_gp_to_g = matxsr_data['transfer_mat'][0]
  
# # #======================================= Solve for psi
# S = M_sig_gp_to_g.transpose()
# #T = np.diag(v_sig_t)
# plt.matshow(np.log(S))
# plt.title('Matxsr elastic scattering')
