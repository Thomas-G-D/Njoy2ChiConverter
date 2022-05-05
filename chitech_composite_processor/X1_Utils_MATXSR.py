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
    
tape        = 'tape41'             # input file
matxsr_data = {}                   # initialize dictionary to hold xs data
matxsr_data['transfer_mat'] = []   # initialize the transfer matrix entry


def EditMatrix(file):
    cf = open('tape41', 'r') #tape[i]
    structure = []   # initialize variable to determine matxsr structure
    xs_matrix = []   # initialize variable to contain elastic scattering information
   # xs        = []   # initialize variable to 
   # groups    = []   # 
    last_line = []
    last_line.append('start')
    for line in cf:
        structure = []  # initialize variable to determine matxsr structure
        xs_matrix = []      # initialize variable to contain elastic scattering information
      #  xs = []
      #  groups = []
        #%% first, we read out the energy group structure
        if line.find('4d') != -1:           # find where energy groups start
            #structure_line = cf.readline().split()    # read the lines for the structure
            structure_line = line.split()
            structure.append(structure_line)
            # this loop reads until 5d appears, which occurs at the end of the energy loops
            while (structure_line[0] != '5d'):              # while still reading structure 
                structure_line = cf.readline().split()  # read structure line
                structure.append(structure_line)        # add to structure variable
            # now we need to format what we just read
            structure = structure[:-1]  # remove the last line that contains 5d
            structure[0] = structure[0][1:] # remove the 4d
            energy_bins = [float(bins) for sublist in structure for bins in sublist] # determine transfer group str
            matxsr_data['energy_bins'] = energy_bins
        #%% next assemble the cumulative xs values
        if line.find('6d') != -1 or last_line[0] == '6d':
            if last_line[0] == '6d':
                structure_line = last_line
                structure.append(structure_line)
                structure.append(line.split()) ### I don't like that I need this line here and not below
            else:
                structure.append(line.split())
                structure_line = cf.readline().split()    # read the lines for the structure
                structure.append(structure_line)
            # this loop reads until 5d appears, which occurs at the end of the energy loops
            while (structure_line[0] != '7d'):              # while still reading structure 
                structure_line = cf.readline().split()  # read structure line
                structure.append(structure_line)        # add to structure variable
            # now we need to format what we just read
            last_line = structure[-1]
            structure = structure[:-1]      # remove the last line that contains 7d
            structure[0] = structure[0][1:] # remove the 6d
            struc_info = [info for sublist in structure for info in sublist] # determine transfer group str
            num_xs = int(len(struc_info)/3)
            info_name = struc_info[0:num_xs]   
            info_grp_start = struc_info[num_xs:2*num_xs]
            info_grp_start = [int(i) for i in info_grp_start]
            info_grp_end = struc_info[2*num_xs:3*num_xs]
            info_grp_end = [int(i) for i in info_grp_end]
            structure = []
            structure.append(last_line)
            while (structure_line[0] != '8d'):              # while still reading structure 
                structure_line = cf.readline().split()  # read structure line
                structure.append(structure_line)        # add to structure variable
            # now we need to format what we just read
            last_line = structure[-1]
            structure = structure[:-1]      # remove the last line that contains 7d
            structure[0] = structure[0][1:]
            xs_info = [info for sublist in structure for info in sublist] # determine transfer group str
            arr_position = 0
            for single_xs in range(num_xs):
                xs_array = np.zeros(172)
                for grp in range(info_grp_start[single_xs], info_grp_end[single_xs]+1):#range(info_num_grps[single_xs]):
                    
                    index_pos = grp+arr_position-info_grp_start[single_xs]
                    # if index_pos == 170:
                    #     test=1
                    grp_nmbr  = grp-1 #+ info_grp_start[single_xs]-1
                    # print(index_pos)
                    # print(grp_nmbr)
                    # print(grp)
                    xs_array[grp_nmbr] = xs_info[index_pos]
                arr_position += grp
                xs_array = [float(i) for i in xs_array]
                matxsr_data[info_name[single_xs]] = xs_array
            continue  ### recent add
            
        #%% starting reading out transfer matrices
        if line.find('8d') != -1 or last_line[0] == '8d':
            structure = []  # initialize variable to determine matxsr structure
            xs_matrix = []      # initialize variable to contain elastic scattering information
           # xs = []
           # groups = []
            if last_line[0] == '8d':
                structure_line = last_line
                structure.append(structure_line)
                structure.append(line.split()) ### I don't like that I need this line here and not below
            else: ### check if I need an append line in here
                structure.append(line.split()) ### this may be a mistake
                structure_line = cf.readline().split()    # read the lines for the structure
                structure.append(structure_line)  
            info_name = structure[0][1]
            if info_name == 'n38':
                test=1 ### missing a structure line
            while structure_line[0] != '9d':              # while still reading structure
                    structure_line = cf.readline().split()  # read structure line
                    structure.append(structure_line)        # add to structure variable
            
            last_line = structure[:-1]
            structure = structure[1:] ### new line
            structure[0] = structure[0][2:] # remove random values
            structure = structure[:-1]      # remove last row (code reads one line too far)
            struc_info = [item for sublist in structure for item in sublist] # determine transfer group structure
            info_grp_start = struc_info[-172:]    # determine which group is the starting group for each set of xs
            info_grp_end = struc_info[:172]    # remove those lines from structure varible
            #info_grp_end = info_grp_end[2:]   # do not need the last 172 variables: this is just energy groups
            
            #xs_info = cf.readline().split()  # continue reading transfer matrix
            #xs_matrix.append(xs_info)  
            
            xs_info = structure_line
            if len(xs_info) < 6: # this if loop handles the case when negative values are present
                                        # if there are negative values, there are no spaces between numbers
                    new_line = ' '.join([str(item) for item in xs_info]) # read line into a string variable
                    for char in range(1, len(new_line)):                    # loop over each character
                        # if the character is negative and not part of scientific notation, add a space
                        if new_line[char] == '-' and new_line[char-1] != 'E' and new_line[char-1] != ' ':
                            new_line = new_line[:char] +' '+new_line[char:]
                            char += 2 ### this line may be unnecessary
                    xs_matrix.append(new_line.split())
            else:
                xs_matrix.append(xs_info)                    # append elastic scattering variable
            ### need to get multiple while loops to stop
            while  len(xs_info) != 0 and xs_info[0] != '8d' and xs_info[0] != '6d':         # while still reading the transfer matrix
                xs_info = cf.readline().split() # read line
                if len(xs_info) < 6: # this if loop handles the case when negative values are present
                                        # if there are negative values, there are no spaces between numbers
                    new_line = ' '.join([str(item) for item in xs_info]) # read line into a string variable
                    for char in range(1, len(new_line)):                    # loop over each character
                        # if the character is negative and not part of scientific notation, add a space
                        if new_line[char] == '-' and new_line[char-1] != 'E' and new_line[char-1] != ' ':
                            new_line = new_line[:char] +' '+new_line[char:]
                            char += 2 ### this line may be unnecessary
                    xs_matrix.append(new_line.split())
                else: # if the line has the expected number of outputs, proceed normally
                    xs_matrix.append(xs_info)
                # if  xs_info[0] == '6d':
                #     test=1
            last_line = xs_matrix[-1]

            
            xs_matrix = xs_matrix[:-1]  # remove last line from list it reads an extra line
            xs_list = [item for sublist in xs_matrix for item in sublist] # put all values into a single list variable
            occurences = xs_list.count('9d') # count how many times 9d occurs in the list variable
            # this loop removes all instances of 9d
            for c in range(occurences):
                xs_list.remove('9d')
            grps_done = 0    # variable to index groups done, this is needed to follow the P band format
            xs_list = [float(xs) for xs in xs_list] # convert strings to floats
            info_grp_end = [int(grp) for grp in info_grp_end]     # convert str to int
            info_grp_start = [int(grp) for grp in info_grp_start]     # convert str to int
            xs_array = []
            xs_array.append(np.zeros([172,172]))
            xs_array.append(np.zeros([172,172]))
            xs_array.append(np.zeros([172,172]))
            xs_array.append(np.zeros([172,172]))
            xs_array.append(np.zeros([172,172]))
            xs_array.append(np.zeros([172,172]))
            xs_array.append(np.zeros([172,172]))
            xs_array.append(np.zeros([172,172]))
            if info_name == 'n38':
                test=1
            # this loop builds the transfer matrix
            for arr_grp in range(len(info_grp_end)):                        # departing group
                for p in range(0,8):                                            # legendre coefficient number
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
            
            if info_name[0] != 'n*' and info_name != 'free':
               matxsr_data['transfer_mat'] += xs_array
            if info_name == 'nelas':
               matxsr_data['transfer_mat'] += xs_array 
               print('nelas')
            if info_name == 'inelas':
               matxsr_data['transfer_mat'] += xs_array 
            if info_name == 'nx':
               matxsr_data['transfer_mat'] += xs_array 
        
        
        
        
        
        
        
        
        
        
     
    cf.close()
    return matxsr_data
data = EditMatrix('tape41')




# M_sig_gp_to_g = np.zeros([172,172])
# #M_sig_gp_to_g[:,:] = test_mat[0][:,:]
  
# #======================================= Solve for psi
# S = M_sig_gp_to_g.transpose()
# #T = np.diag(v_sig_t)
# plt.matshow(np.log(S))
# plt.title('Matxsr elastic scattering')
