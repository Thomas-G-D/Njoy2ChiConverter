# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 00:55:42 2022

@author: jean.ragusa
"""
import numpy as np

import sys, os
sys.path.insert(0, os.path.realpath('../chitech_composite_processor'))
import Utils_ChiTechCombiner
import Utils_Info
import Utils_NjoySpectrumPlotter
import matplotlib.pyplot as plt
from collections import defaultdict
plt.close('all')

def index(filename, lst):                                # define index function
        with open(filename, 'r') as infile:              # open file of interest
            lines = [line.split() for line in infile]    # break lines into individual words
        word2linenumbers = defaultdict(list)             # return a dictionary object with line numbers
    
        for linenumber, line in enumerate(lines, 1):     # for each line
            for word in line:                            # for each word in the line
                if word in lst:                          # if that word is in the variable lst
                    word2linenumbers[word].append(linenumber)  # add it to the dictionary object
        return word2linenumbers                          # return dictionary object   
    
def assemble_data(filename, target_value, skip_number, multiplier):
    with open(filename, 'r') as infile:              # open file of interest
        lines = [line.split() for line in infile]    # break lines into individual words
    Output_values = np.zeros((173,3))                # initialize output array
    for linenumber, line in enumerate(lines, 1):     # for each line
        if linenumber == target_value:               # if the line number matches the target value
            for x in range(0,173):                   # for each of the 172 XS
                #Output_values.append(lines[skip_number + linenumber +x*multiplier])  # add value to output
                Output_values[x][0] = lines[skip_number + linenumber +x*multiplier][0] # energy bin upper bound
                Output_values[x][1] = lines[skip_number + linenumber +x*multiplier][1] # flux tally value
                Output_values[x][2] = lines[skip_number + linenumber +x*multiplier][2] # uncertainty
    return Output_values







Isotopes = ['Al27', 'Ar36', 'Ar38', 'Ar40', 'C12', 'Cr52', 'Cu63', 'Cu65', 'Fe56', 'H1', 'Mg24', 'Mg25', 'Mg26', 'Mn55', 'N14', 'N15',
            'O16', 'O17', 'Si28', 'Si30', 'Ti48', 'Zn64', 'Zn66', 'Zn67', 'Zn68', 'Zn70', 'Cd106', 'Cd108', 'Cd110', 'Cd111', 'Cd112', 'Cd113', 'Cd114', 'Cd116']
#Isotopes = ['Al27']
difference            = {}
normalized_difference = {}
for i in Isotopes:    
    chixs_fullpath = []
    chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/' +i+'.cxs')
    N_density = []
    N_density.append(1.)
    data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath, N_density)
    print(data.keys())


    source_def = {}
    source_def["particle_type"] = 'neutron'
    source_def["energy"] = 18.49

    outp = Utils_Info.InfiniteMediumSpectrum(data, source_def, plot=True)

    Utils_NjoySpectrumPlotter.Njoy_spectrum_plotter(outp, './')

    #A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/'+i+'.txt')
    MCNP_filename = '../../outputs/ENDF-B-VIII-0/172gxs/MCNP/'+i+'XSout.txt'
    
    MCNP_tally = index(MCNP_filename, ['1tally', '14', '24'])
    
    # this loop finds the line corresponding to 1tally and 14
    tally_14_line = 0                                          # initialize value
    for k in range(0,len(MCNP_tally['1tally'])):               # for each line corresponding to 1tally
        for j in range(0,len(MCNP_tally['14'])):               # loop over all lines corresponding to 14
            if MCNP_tally['1tally'][k] == MCNP_tally['14'][j]: # if the line numbers correspond
                tally_14_line = MCNP_tally['14'][j]            # set as the output value
                
    # build a dictionary variable containing the two line numbers found above
    MCNP_dict = {'tally_14': tally_14_line}#, 'tally_24' : tally_24_line}

    # assemble the MCNP flux starting 10 lines after the tally 14 line appears
    MCNP_data = assemble_data(MCNP_filename, MCNP_dict['tally_14'],9, 1)


    mcnpE = MCNP_data[:,0]
    mcnpF = MCNP_data[:,1] * 4e9 


    dE = np.diff(np.insert(mcnpE,0,0.))
    spectrum = mcnpF/dE

    E = []
    F = []
    for g in range(MCNP_data.shape[0]-1):
        E += [mcnpE[g], mcnpE[g+1]]
        F += [spectrum[g+1], spectrum[g+1]]

    E = np.asarray(E)
    F = np.asarray(F)
    
    ### this section I addeed
    neutron_spectrum = outp[0][1] ### I added 
    #difference[i] = (F-neutron_spectrum)/F
    F_norm = F/np.sum(F)
    neutron_spectrum_norm = neutron_spectrum/np.sum(neutron_spectrum)
    #normalized_difference[i] = (F_norm-neutron_spectrum_norm)/F_norm
    #F = F/np.sum(F)  ### I added

    fig_n = plt.gcf().number
    print(fig_n)
    fig = plt.figure(fig_n)
    ax_list = fig.axes
    print(ax_list)
    ax_list[0].semilogy(E, F, label='mcnp')
    ax_list[1].loglog(E, F, label='mcnp')
    ax_list[0].legend()
    ax_list[1].legend()
    plt.show()
    #plt.savefig(i+'_normalized.png')
    plt.savefig(i+'.png')
    
    ### I added this
    E_mids = np.zeros(172)
    difs   = np.zeros(172)
    norm_difs = np.zeros(172)
    for mid in range(0, len(E_mids)):
        E_mids[mid] = mid #(E[2*mid]+E[2*mid+1])/2
        difs[mid] = (F[mid*2]-neutron_spectrum[mid*2])/F[mid*2]
        norm_difs[mid] = (F_norm[mid*2]-neutron_spectrum_norm[mid*2])/F_norm[mid*2]
        
    difference[i] = difs
    normalized_difference[i] = norm_difs 
    plt.figure()
    plt.plot(E_mids, difference[i], label = 'difference')
    plt.plot(E_mids, normalized_difference[i], label = 'normalized difference')
    plt.legend()
    plt.xlabel('Energy Bin')
    plt.ylabel(r'$\frac{mcnp - \chi}{mcnp}$')
    plt.legend()
    plt.grid()
    plt.savefig(i+'_difference.png')


