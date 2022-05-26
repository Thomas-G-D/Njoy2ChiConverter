# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 12:26:25 2022

@author: tdeguire
"""


import numpy as np
import matplotlib.pyplot as plt


def EditMatrix(file):
    
    
    matxsr_data = {} # initialize dictionary to hold xs data
    
    matxsr_data["sigma_t"]         = [] # done
    matxsr_data["sigma_a"]         = [] # maybe done
    matxsr_data["sigma_s"]         = [] # done
    matxsr_data["sigma_f"]         = [] # done
    matxsr_data["sigma_heat"]      = [] # done
    
    ########## still have no idea on these
    matxsr_data["nu_total"]        = []
    matxsr_data["nu_prompt"]       = []
    matxsr_data["nu_delayed"]      = []
    matxsr_data["chi_prompt"]      = []
    matxsr_data["chi_delayed"]     = []
    matxsr_data["decay_constants"] = []
    matxsr_data["gamma"]           = []
    
    #############
    matxsr_data["inv_velocity"]    = [] # done

    matxsr_data["transfer_matrices"] = [] # done
    matxsr_data["transfer_matrices_sparsity"] = [] # done

    matxsr_data["neutron_gs"] = [] # done
    matxsr_data["gamma_gs"]   = []
    matxsr_data['test_sig_s'] = np.zeros(172)
    
    
    
    
    
    
    n_to_n_transfer_initialize = False
    g__to_g_transfer_initialize = False
    n_to_g_transfer_initialize = False
    cf = open(file, 'r') # tape[i]
    structure = []   # initialize list to determine matxsr structure
    xs_matrix = []   # initialize list to contain elastic scattering information
    last_line = []   # initialize a list for containg lines at end of a specific xs
    last_line.append('start') # need this variable initialized to prevent error
    data_type = 0    # initialize the data type variable
    blocks_initialized = False
    
    ## this loop goes over every line in the file and assembles the data into a dictionary
    for line in cf:
        
        ## initialize key variables
        structure = []  # initialize variable to determine matxsr structure
        xs_matrix = []  # initialize variable to contain elastic scattering information
        
        #%% find general file structure info from card 1
        if line.find('1d') != -1:             # read card 1
            
            card_1         = line.split()     # card 1 information
            particle_types = int(card_1[1])        # type of particle
            data_types     = int(card_1[2])   # number of data types
            num_blocks     = int(card_1[-1])  # number of blocks of code
        
        #%% read out information about the data types from card 3
        if line.find('3d') != -1:
            
            card_3          = line.split()
            particle_names = []
            for i in range(particle_types):
                particle_names.append(card_3[1+i])
            data_type_names = []                         # initialize name of the data types 
            
            for i in reversed(range(data_types)):
                data_type_names.append(card_3[-(2+i)])   # extract names of data types from text
            ### want to read out until 3d but not sure how everything works yet if multiple particles!!!
            ### want to read out number energy groups from here
            structure_line = card_3
            while (structure_line[0] != '4d'):          # while still reading structure 
                structure_line = cf.readline().split()  # read structure line
                structure.append(structure_line)        # add to structure variable
            last_line    = structure[-1]
            structure    = structure[:-1]   # remove the last line that contains 5d, one line too far read by while loop
            struc_info     = [info for sublist in structure for info in sublist] # read structure into one list
            
            energy_groups = []
            for i in range(particle_types):
                energy_groups.append(int(struc_info[i]))
            initial_particle = []
            final_particle = []
            for i in range(len(data_type_names)):
                initial_particle.append(struc_info[particle_types+i])
                final_particle.append(struc_info[particle_types+len(data_type_names)+i])
            
        
        #%% Read out the energy group structure from card 4
        ### may need a last line statement depending on card 3
        if line.find('4d') != -1 or last_line[0] == '4d':            # find where energy groups start
            structure = []    
            if last_line[0] == '4d':
                    structure_line = last_line              # take the last line read from the previous loop
                    structure.append(structure_line)        # append it to the structure variable
                    #structure.append(line.split())          # then add the current line being read
            else:
                structure.append(line.split())          # add the line being read currently
                structure_line = cf.readline().split()  # read the next line
            # structure_line = line.split()    # split the line being read
            # structure.append(structure_line) # add it to the structure
            
            # this loop reads until 5d appears, which occurs at the end of the energy loops
            while (structure_line[0] != '5d'):          # while still reading structure 
                structure_line = cf.readline().split()  # read structure line
                structure.append(structure_line)        # add to structure variable
            
            # now we need to format what we just read
            last_line    = structure[-1]
            structure    = structure[:-1]   # remove the last line that contains 5d, one line too far read by while loop
            structure[0] = structure[0][1:] # remove the 4d string
            energy_bins  = [bins for sublist in structure for bins in sublist] # determine transfer group str
            occurences = energy_bins.count('4d') # count how many times 9d occurs in the list variable
            # this loop removes all instances of 4d
            for c in range(occurences):
                energy_bins.remove('4d')
            
            bins_done = 0
            #n_groups = []
            for nrg in range(particle_types):    
                matxsr_data[particle_names[nrg]+'_energy_bins'] = energy_bins[bins_done:(bins_done+energy_groups[nrg]+1)]  # save energy bins to dictionary
                bins_done = energy_groups[nrg]+1 
                
        
        #%% read out information about the data types from card 5
        ### could potentially do more with this info but may not be necessary as
        ### current method seems to work well
        if line.find('5d') != -1 or last_line[0] == '5d':
            
            data_type_struct = {}                                     # initialize dictionary
            card_5 = []
            data_blocks = []                                                
            
            for types in range(data_types):                           # for the number of data types
                data_dict = {}                                        # initialize sub dictionary  
                card_5.append(cf.readline().split())                  # red line 
                data_dict['temp']     = card_5[types][0]              # extract ambient temperature
                data_dict['dilution_fact'] = card_5[types][1]         # extract dilution factor
                data_dict['num_vecs'] = int(card_5[types][3])         # extract number of vectors
                data_dict['num_mats'] = int(card_5[types][4])         # extract number of matrices 
                data_dict['starting'] = int(card_5[types][5])         # extract starting block
                data_type_struct[data_type_names[types]] = data_dict  # add sub dictionary to dictionary
                data_blocks.append(int(card_5[types][5]))
            last_line[0] = 0                                          # reset last line or this runs twice 

        
        #%% next assemble the xs values from card 6. Card 6 contains vector structure info, card 7 contains vector values
        if line.find('6d') != -1 or last_line[0] == '6d':
            ########################### do i need to reset the structure???
            # this if loop initializes how the xs data will be reading depending
            # on how this portion of the code was entered
            if last_line[0] == '6d':
                structure_line = last_line              # take the last line read from the previous loop
                structure.append(structure_line)        # append it to the structure variable
                structure.append(line.split())          # then add the current line being read #########check this line
            else:
                structure.append(line.split())          # add the line being read currently
                structure_line = cf.readline().split()  # read the next line
                if len(structure_line) == 0:
                    structure_line = cf.readline().split()  # read structure line
                    structure.append(structure_line)        # add to structure variable
                else:
                    structure.append(structure_line)        # add it to the matrix
            
            if blocks_initialized == False:
                blocks_initialized = True
                blocks = 0
            else:
                blocks += 1
            if blocks in data_blocks:
                data_type = data_blocks.index(blocks)
                
            # this loop reads until 7d appears, which occurs at the end of xs group structure
            while (structure_line[0] != '7d'):          # while still reading structure 
                structure_line = cf.readline().split()  # read structure line
                structure.append(structure_line)        # add to structure variable
            
            # now we need to format what we just read
            last_line      = structure[-1]       # save the last line; it is needed to the xs's
            structure      = structure[:-1]      # remove the last line that contains 7d
            structure[0]   = structure[0][1:]    # remove the 6d
            struc_info     = [info for sublist in structure for info in sublist] # read structure into one list
            num_xs         = data_type_struct[data_type_names[data_type]]['num_vecs']
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
            last_line    = structure[-1]       # save the last line; it is needed for the next xs
            structure    = structure[:-1]      # remove the last line that contains 9d
            #structure[0] = structure[0][1:]    # remove the 8d
            xs_info      = [info for sublist in structure for info in sublist] # read all xs's into one list
            
            occurences = xs_info.count('7d') # count how many times 9d occurs in the list variable
            # this loop removes all instances of 4d
            for c in range(occurences):
                blocks += 1
                xs_info.remove('7d')

            # if blocks in data_blocks:
            #     data_type = data_blocks.index(blocks)
            
            # this loop assembles lists for each xs and adds them to the dictionary
            arr_position = 0 # initialize position in the xs array; needed to move from one xs to the next
            for single_xs in range(num_xs):       # for each xs
                xs_array = np.zeros(energy_groups[int(final_particle[data_type])-1])     # initialize the array to be of length 172 for 172 energy groups
                for grp in range(info_grp_start[single_xs], info_grp_end[single_xs]+1): # for each group
                    index_pos          = grp+arr_position-info_grp_start[single_xs]     # index position in the array
                    grp_nmbr           = grp-1                                          # grp number
                    xs_array[grp_nmbr] = xs_info[index_pos]                             # add value to xs array
                arr_position += (info_grp_end[single_xs]+1-info_grp_start[single_xs])   # increment starting position
                xs_array = [float(i) for i in xs_array]        # convert values from strings to floats
                matxsr_data[info_name[single_xs]] = xs_array   # add xs to matxsr dictionary
            #data_type += 1                                     # move to next data type, needed to change num_xs variable
            
            continue
            
        #%% starting reading out transfer matrices from card 8. card 8 contains matrix structure info, card 9 contains matrix values
        if line.find('8d') != -1 or last_line[0] == '8d':
            structure = []  # initialize variable to determine matxsr structure
            xs_matrix = []  # initialize variable to contain elastic scattering information
            const_block = False
            
            if blocks_initialized == False:
                blocks_initialized = True
                blocks = 0
            else:
                blocks += 1
            if blocks in data_blocks:
                data_type = data_blocks.index(blocks)

            if last_line[0] == '8d':
                structure_line = last_line             # take the last line read from the previous loop
                structure.append(structure_line)       # append it to the structure variable
                structure.append(line.split())         # then add the current line being read ########### check this line
            else:
                structure.append(line.split())          # append the structure line with the split line
                structure_line = cf.readline().split()  # read the next line in the file
                structure.append(structure_line)        # add it to the matrix     
            
            info_name = structure[0][1]                 # xs type
            init_groups = energy_groups[int(initial_particle[data_type])-1] ########## make sure this is in the right place
            final_groups = energy_groups[int(final_particle[data_type])-1] ########## make sure this is in the right place

            
            # this loop extracts the structure for the xs matrix
            while structure_line[0] != '9d':                # while still reading structure
                    structure_line = cf.readline().split()  # read structure line
                    structure.append(structure_line)        # add to structure variable
            
            last_line      = structure[-1]          # save the last line, it is needed for building the xs matrix
            structure      = structure[1:]          # remove the first line, it is unnecessary
            Pn_order       = int(structure[0][0])   # number of legendre coefficients
            const_groups   = int(structure[0][1])   # number of constant groups
            structure[0]   = structure[0][2:]       # remove unneeded values
            structure      = structure[:-1]         # remove last row (code reads one line too far)
            struc_info     = [item for sublist in structure for item in sublist] # determine transfer group structure
            info_grp_start = struc_info[-final_groups:] # determine which group is the starting group for each set of xs
            info_grp_end   = struc_info[:final_groups]  # remove those lines from structure varible
            
            # check if the transfer matrix has already been initialized
            if n_to_n_transfer_initialize == False and data_type_names[data_type] == 'nscat':
                n_to_n_transfer_mat = np.zeros([init_groups,final_groups,Pn_order]) # initialize transfer matrix
                n_to_n_transfer_initialize = True
            if g_to_g_transfer_initialize == False and data_type_names[data_type] == 'gscat':
                g_to_g_transfer_mat = np.zeros([init_groups,final_groups,Pn_order]) # initialize transfer matrix
                g_to_g_transfer_initialize = True
            if n_to_g_transfer_initialize == False and data_type_names[data_type] == 'gscat':
                n_to_g_transfer_mat = np.zeros([init_groups,final_groups,Pn_order]) # initialize transfer matrix
                n_to_g_transfer_initialize = True
            
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
            while  len(xs_info) != 0 and xs_info[0] != '8d' and xs_info[0] != '6d' and xs_info[0] != '10d':         # while still reading the transfer matrix
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
                blocks += 1
                xs_list.remove('9d')
            
            # if blocks in data_blocks:
            #     data_type = data_blocks.index(blocks)
            
            
            xs_list        = [float(xs) for xs in xs_list]         # convert strings to floats
            info_grp_end   = [int(grp) for grp in info_grp_end]    # convert str to int
            info_grp_start = [int(grp) for grp in info_grp_start]  # convert str to int
            xs_array = np.zeros([init_groups, final_groups, Pn_order])    # initialize array to store xs info    
            
            if len(last_line) != 0 and last_line[0] == '10d':
                const_block = True
                constant_mat = []
                constant_mat.append(last_line)
                constant_info = last_line
                while  len(constant_info) != 0 and constant_info[0] != '8d' and constant_info[0] != '6d' and constant_info[0] != '9d':         # while still reading the transfer matrix
                    constant_info = cf.readline().split() # read line
                    if len(constant_info) < 6:    # this if loop handles the case when negative values are present
                                            # if there are negative values, there are no spaces between numbers
                        new_line = ' '.join([str(item) for item in constant_info]) # read line into a string variable
                        for char in range(1, len(new_line)):                    # loop over each character
                            # if the character is negative and not part of scientific notation, add a space
                            if new_line[char] == '-' and new_line[char-1] != 'E' and new_line[char-1] != ' ':
                                new_line = new_line[:char] +' '+new_line[char:]
                        constant_mat.append(new_line.split())
                    else: # if the line has the expected number of outputs, proceed normally
                        constant_mat.append(constant_info)
            
                last_line     = constant_mat[-1]   # save the last line, this may be needed for another xs
                constant_mat  = constant_mat[:-1]  # remove last line from list it reads an extra line
                const_list    = [item for sublist in constant_mat for item in sublist] # put all values into a single list variable
                ### convert constant mat to a list and then remove 10d or anything else and save last line
                occurences = const_list.count('10d') # count how many times 9d occurs in the list variable
                # this loop removes all instances of 9d
                for c in range(occurences):
                    blocks +=1 
                    const_list.remove('10d')
                spec = np.zeros(final_groups)
                prod_xs = np.zeros(const_groups)
                for gp in range(len(spec)):
                    #gp = n_groups-1-g
                    spec[gp] = const_list[gp]
                for g in range(len(prod_xs)):
                    #grp = const_groups-1-g
                    prod_xs[g] = const_list[len(spec)+g]
                
                grps_done = 0    # variable to index groups done, this is needed to follow the P band format
                # this loop builds the transfer matrix
                ###### these do not line up need to re-evaluate
                start_grp = (init_groups-const_groups)
                for arr_grp in range(final_groups):                        # departing group
                    for p in range(Pn_order):                                        # legendre coefficient number
                        for dep_grp in range(start_grp,init_groups): # arrival group
                            # print('Col grp is '+str(arr_grp))
                            # print('row grp is '+str(dep_grp))
                            #row = info_grp_start[arr_grp]-1
                            xs_array[dep_grp, arr_grp, p] = spec[arr_grp]*prod_xs[dep_grp-start_grp]
                    grps_done += info_grp_end[arr_grp]    # track how many groups have already been done    
                
            
            

            grps_done = 0    # variable to index groups done, this is needed to follow the P band format
            # this loop builds the transfer matrix
            for arr_grp in range(len(info_grp_end)):                        # departing group
                for p in range(Pn_order):                                        # legendre coefficient number
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
                        xs_array[row-dep_grp, arr_grp, p] = xs_list[grps_done*Pn_order+p*info_grp_end[arr_grp]+dep_grp]
                grps_done += info_grp_end[arr_grp]    # track how many groups have already been done
            
            matxsr_data[info_name+'_mat'] = xs_array  # add the matrix containing the xs info to the dictionary
            
            # set a loop name as the corresponding dictionary name
            loop_name = info_name + '_mat'
            # add all xs matrices to a total transfer matrix except s(alpha,beta) and free gas
            ### need to see if there is a way to read this out as the s(alpha,beta) change based on name
            ### have this loop only work for data type nscat and gscat but keep them separate
            if data_type_names[data_type] == 'nscat':
                # print(data_type)
                n_transfer_mat += xs_array
            if data_type_names[data_type] == 'gscat':
                # print(data_type)
                g_transfer_mat += xs_array
    
    # make a new variable in the dicitonary to be the total transfer matrix
    max_num_moms = max([np.shape(n_transfer_mat)[2], np.shape(g_transfer_mat)[2]])
    G = len(matxsr_data['neutron_gs']) + len(matxsr_data['gamma_gs'])
    transfer_mats = []
    for m in range(0, max_num_moms):
        transfer_mats.append(np.zeros((G, G)))
        for dep_grp in range(G):
            for arr_grp in range(G):
                transfer_mats[m][dep_grp, arr_grp]
    
    ### need to see if ng is included  in this
    #### may not need the other transfer mat loop below
    # matxsr_data['transfer_mats'] = np.copy(n_transfer_mat)

    #%% apply free gas and s(alpha,beta) treatment
    # apply free gas
    # for P0 through P_max
    for p in range(np.shape(n_transfer_mat)[2]):
        # for each departing group
        for dep_grp in range(np.shape(n_transfer_mat)[0]):
            # for each arrival group
            for arr_grp in range(np.shape(n_transfer_mat)[1]):
              # if the free gas value is non-zero, updated the transfer matrix
              if abs(matxsr_data['free_mat'][dep_grp, arr_grp, p]) > 1e-18:
                matxsr_data['transfer_mats'][dep_grp, arr_grp, p] = matxsr_data['free_mat'][dep_grp, arr_grp, p]
    
    # add a new matrix to the dictionary for s(alpha,beta)
    # repeat the same process as for free gas
    # matxsr_data['transfer_mat_sab'] = np.copy(matxsr_data['transfer_mats'])
    # for p in range(Pn_order):
    #     for dep_grp in range(n_groups):
    #         for arr_grp in range(n_groups):
    #           if abs(matxsr_data['fe_mat'][dep_grp, arr_grp, p]+ matxsr_data['fe$_mat'][dep_grp, arr_grp, p]) > 1e-18:
    #             matxsr_data['transfer_mat_sab'][dep_grp, arr_grp, p] = matxsr_data['fe_mat'][dep_grp, arr_grp, p] + matxsr_data['fe$_mat'][dep_grp, arr_grp, p]
    
    # add variables that are needed to generate spectra plots
    
    for nrg in range(particle_types):    
        matxsr_data[particle_names[nrg]+'_energy_bins'].reverse()
        for grp in range(len(matxsr_data[particle_names[nrg]+'_energy_bins'])-1):
            lo_bound  = matxsr_data[particle_names[nrg]+'_energy_bins'][grp]     # low boundary of energy bin
            hi_bound  = matxsr_data[particle_names[nrg]+'_energy_bins'][grp+1]   # high boundary of energy bin
            if particle_names[nrg] == 'n':
                matxsr_data['neutron_gs'].append([float(grp), lo_bound, hi_bound])
            elif particle_names[nrg] == 'g':
                matxsr_data['gamma_gs'].append([float(grp), lo_bound, hi_bound])
    
    
    
    
    
    
   
    # need to adjust form of the transfer materices to be list of lists not 3d list
    array_shape = np.shape(matxsr_data['transfer_mats'])
    for p in range(array_shape[2]):
        matxsr_data['transfer_matrices'].append(matxsr_data['transfer_mats'][:,:,p])
    transfer_mats_nonzeros = []
    
    for m in range(0, array_shape[2]):
        mat_non_zeros = []
        for gprime in range(0, G):
            non_zeros = []
            for g in range(0, G):
                print(gprime)
                print(g)
                if np.abs(matxsr_data['transfer_matrices'][m][gprime, g]) > 1.0e-16:
                    non_zeros.append(g)
            mat_non_zeros.append(non_zeros)
        transfer_mats_nonzeros.append(mat_non_zeros)    
    matxsr_data['transfer_matrices_sparsity'] = transfer_mats_nonzeros
    
        
    matxsr_data['sigma_t']    = matxsr_data['ntot0']    # adjust naming for total xs
    matxsr_data['sigma_s']    = np.sum(matxsr_data['transfer_matrices'][0], axis = 1)
    ### siga a may not be correct may be sig_t-nelas - ninel
    # matxsr_data['sigma_a']    = matxsr_data['sigma_t'] - matxsr_data['sigma_s']
    if 'nelas' in matxsr_data.keys() and 'ninel' in matxsr_data.keys():
        matxsr_data['sigma_a']    =np.array(matxsr_data['sigma_t']) - np.array(matxsr_data['nelas']) - np.array(matxsr_data['ninel']) #- np.array(matxsr_data['nna']) -np.array(matxsr_data['nnp'])
    elif 'nelas' in matxsr_data.keys() and 'ninel' not in matxsr_data.keys():
        matxsr_data['sigma_a']    =np.array(matxsr_data['sigma_t']) - np.array(matxsr_data['nelas'])
    elif 'nelas' not in matxsr_data.keys() and 'ninel' in matxsr_data.keys():
        matxsr_data['sigma_a']    =np.array(matxsr_data['sigma_t']) - np.array(matxsr_data['ninel'])
    matxsr_data['sigma_heat'] = matxsr_data['heat']     # adjust naming for heating xs
    matxsr_data['inv_velocity']    = np.copy(matxsr_data['invel'])
    
    if 'nftot' in matxsr_data.keys():
        matxsr_data["sigma_f"] = np.copy(matxsr_data["nftot"])
        
    # matxsr_data['transfer_matrices_sab'] = []               # initialize new key
    # for p in range(8):
    #     matxsr_data['transfer_matrices_sab'].append(matxsr_data['transfer_mat_sab'][:,:,p])
    cf.close()
    
    problem_description = {}
    if matxsr_data['neutron_gs'] != []:
        if matxsr_data['gamma_gs'] != []:
            problem_description['problem_type'] = "Neutron + Gamma"
        elif matxsr_data['gamma_gs'] == []:
            problem_description['problem_type'] = "Neutron only"
    elif matxsr_data['gamma_gs'] != []:
        problem_description['problem_type'] = "Gamma only"
    else:
        problem_description['problem_type'] ="Unknown"
        print(problem_description)
        raise Exception('unknown particle type')
        
    problem_description['G_n'] = len(matxsr_data['neutron_gs'])
    problem_description['G_g'] = len(matxsr_data['gamma_gs'])
    problem_description['isotope'] = card_3[-1]
    return matxsr_data, problem_description
        
matxsr_data2, problem_description = EditMatrix('../5_9 meeting/NJOY_files/tape41_O16_w_gamma')

### need to figure out how to read card 3 with gamma outputs
### need to figure out how to handle gammas at all
### need to figure out a way to not have s(alpha,beta) hard coded

