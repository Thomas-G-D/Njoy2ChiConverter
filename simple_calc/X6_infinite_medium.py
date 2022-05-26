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
import X7_Utils_MATXSR as matxsr
plt.close('all')





def index(filename, lst):                            # define index function
        with open(filename, 'r') as infile:   # open file of interest
            lines = [line.split() for line in infile]    # break lines into individual words
        word2linenumbers = defaultdict(list)             # return a dictionary object with line numbers
    
        for linenumber, line in enumerate(lines, 1):     # for each line
            for word in line:                            # for each word in the line
                if word in lst:                          # if that word is in the variable lst
                    word2linenumbers[word].append(linenumber)  # add it to the dictionary object
        return word2linenumbers                          # return dictionary object   
def assemble_data(filename, target_value, skip_number, multiplier):
    with open(filename, 'r') as infile:   # open file of interest
        lines = [line.split() for line in infile]    # break lines into individual words
    Output_values = np.zeros((173,3))                               # initialize output array
    for linenumber, line in enumerate(lines, 1):     # for each line
        if linenumber == target_value:               # if the line number matches the target value
            for x in range(0,173):                   # for each of the 172 XS
                #Output_values.append(lines[skip_number + linenumber +x*multiplier])  # add value to output
                Output_values[x][0] = lines[skip_number + linenumber +x*multiplier][0]
                Output_values[x][1] = lines[skip_number + linenumber +x*multiplier][1]
                Output_values[x][2] = lines[skip_number + linenumber +x*multiplier][2]
    return Output_values
cf = open('../5_9 meeting/test_matxsr_files/Plots/Information.txt', 'w')
# cf = open('../5_9 meeting/test_matxsr_files/Plots/Distributed_sources/Information.txt', 'w')
cf.write("{:25}".format('Isotope'))
cf.write("{:40}".format('Value'))
cf.write("{:25}".format('MATXSR'))
cf.write("{:25}".format('CXS'))
cf.write('\n')

distributed_source = True

if distributed_source is True:
    filename = ['Fe56', 'Al27', 'C12', 'Cd106', 
            'Cd108', 'Cd110', 'Cd111', 'Cd112', 
            'Cd113', 'Cd114', 'Cd116', 'Cr52', 
            'H1', 'He3', 'N14', 'Ni58',
            'O16', 'Si28', 'Ar36', 'Ar38',
            'Ar40', 'Cu63', 'Cu65', 'Mg26',
            'Mn55', 'Ti48', 'Zn64', 'Zn66',
            'Zn67', 'Zn68', 'Zn70']
else:
    filename = ['Fe56', 'Al27', 'C12', 'Cd106', 
                'Cd108', 'Cd110', 'Cd111', 'Cd112', 
                'Cd113', 'Cd114', 'Cd116', 'Cr52', 
                'H1', 'He3', 'N14', 'Ni58',
                'O16', 'Si28']

filename = ['Fe56']
for iso in range(len(filename)):
    
    print('Running for '+filename[iso])  # print out what isotope is being processed
    
    ## add a path to the CXS file
    chixs_fullpath = []
    chixs_fullpath.append('../5_9 meeting/test_matxsr_files/'+filename[iso]+'.cxs')
    
    ## add the information needed to assemble the information from the CXS file into a dictionary
    N_density = []
    N_density.append(1.)
    data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath, N_density) # dictionary of cxs information
    
    ## assemble information from matxsr file
    matxsr_data = matxsr.EditMatrix('../5_9 meeting/test_matxsr_files/tape41_'+filename[iso]) # dictionary of matxsr information
    
    ## build the source term
    source_def = {}                                # initialize
    source_def["particle_type"] = 'neutron'        # set particle type
    source_def["energy"] = 18.49                   # energy of source in MeV
    
    ## use the source information to calculate the infinite medium spectrum
    outp_cxs    = Utils_Info.InfiniteMediumSpectrum(data, source_def, distributed_source, plot=True)
    outp_matxsr = Utils_Info.InfiniteMediumSpectrum(matxsr_data, source_def, distributed_source, plot=False)
    
    
    compare = [] 
    for i in range(8):
        compare.append(data['transfer_matrices'][i] - matxsr_data['transfer_matrices'][i])
    # now update to do matxsr_sab
    # for p in range(len(matxsr_data['transfer_mat_sab'][0][0])):
    #     matxsr_data['transfer_matrices'][p] = np.copy(matxsr_data['transfer_mat_sab'][:,:,p])
    outp_matxsr_sab = Utils_Info.InfiniteMediumSpectrum(matxsr_data, source_def, plot=False)
    Utils_NjoySpectrumPlotter.Njoy_spectrum_plotter(outp_cxs, './')
    
    
    fig_n = plt.gcf().number
    fig = plt.figure(fig_n)
    ax_list = fig.axes
    #print(ax_list)
    #ax_list[0].semilogy(outp_matxsr_sab[0][0], outp_matxsr_sab[0][1], label = 'matxsr sab')
    #ax_list[1].loglog(outp_matxsr_sab[0][0], outp_matxsr_sab[0][1], label = 'matxsr sab')
    ax_list[0].semilogy(outp_matxsr[0][0], outp_matxsr[0][1], label = 'matxsr free')
    ax_list[1].loglog(outp_matxsr[0][0], outp_matxsr[0][1], label = 'matxsr free')
    plt.legend()
    #plt.show()
    # plt.savefig('../5_9 meeting/test_matxsr_files/Plots/'+filename+'_spectrum.png')
    
    MCNP_filename = 'MCNP/'+filename[iso]+'XSout.txt'
    MCNP_flux = index(MCNP_filename, ['1tally', '14', '24'])
        
    # this loop finds the line corresponding to 1tally and 14
    tally_14_line = 0                                         # initialize value
    for i in range(0,len(MCNP_flux['1tally'])):               # for each line corresponding to 1tally
        for j in range(0,len(MCNP_flux['14'])):               # loop over all lines corresponding to 14
            if MCNP_flux['1tally'][i] == MCNP_flux['14'][j]:  # if the line numbers correspond
                tally_14_line = MCNP_flux['14'][j]            # set as the output value
    MCNP_dict = {'tally_14': tally_14_line}#, 'tally_24' : tally_24_line}
    
    # assemble the MCNP flux starting 10 lines after the tally 14 line appears
    MCNP_data = assemble_data(MCNP_filename, MCNP_dict['tally_14'],9, 1)
    
    
    # MCNP_filename = 'MCNP/Distributded_src/'+filename[iso]+'.txt'
    # MCNP_flux = index(MCNP_filename, ['1tally', '104', '24'])
        
    # # this loop finds the line corresponding to 1tally and 14
    # tally_14_line = 0                                         # initialize value
    # for i in range(0,len(MCNP_flux['1tally'])):               # for each line corresponding to 1tally
    #     for j in range(0,len(MCNP_flux['104'])):               # loop over all lines corresponding to 14
    #         if MCNP_flux['1tally'][i] == MCNP_flux['104'][j]:  # if the line numbers correspond
    #             tally_14_line = MCNP_flux['104'][j]            # set as the output value
          
    # # build a dictionary variable containing the two line numbers found above
    # MCNP_dict = {'tally_104': tally_14_line}#, 'tally_24' : tally_24_line}
    
    # # assemble the MCNP flux starting 10 lines after the tally 14 line appears
    # MCNP_data = assemble_data(MCNP_filename, MCNP_dict['tally_104'],9, 1)
   
    
    
    
    
    
    
    
    # test = {}
    
    mcnpE = MCNP_data[:,0]
    mcnpF = MCNP_data[:,1] * 4e9 #4e9
    
    dE = np.diff(np.insert(mcnpE,0,0.))
    ### why do we divide by dE
    spectrum = mcnpF/dE
    
    E = []
    F = []
    for g in range(MCNP_data.shape[0]-1):
        E += [mcnpE[g], mcnpE[g+1]]
        F += [spectrum[g+1], spectrum[g+1]]
    
    E = np.asarray(E)
    F = np.asarray(F)
    #neutron_spectrum = outp_cxs[0][1]
    #test['Al27'+'norm'] = (F-neutron_spectrum)/F
    #F = F/np.sum(F)
    
    fig_n = plt.gcf().number
    #print(fig_n)
    fig = plt.figure(fig_n)
    ax_list = fig.axes
    #print(ax_list)
    ax_list[0].semilogy(E, F, label='mcnp')
    ax_list[1].loglog(E, F, label='mcnp')
    plt.legend()
    plt.suptitle('Neutron spectrum for '+filename[iso])
    #plt.show()
    # plt.savefig('../5_9 meeting/test_matxsr_files/Plots/'+filename[iso]+'_spectrum.png', bbox_inches='tight')
    # plt.savefig('../5_9 meeting/test_matxsr_files/Plots/Distributed_sources/'+filename[iso]+'_spectrum.png', bbox_inches='tight')
    
    #matxsr_diff_sab = (outp_matxsr_sab[0][1]-F)/outp_matxsr_sab[0][1]
    matxsr_diff_free = (outp_matxsr[0][1]-F)/F#outp_matxsr[0][1]
    cxs_diff = (outp_cxs[0][1]-F)/F#outp_cxs[0][1]
    plt.figure()
    plt.semilogx(E,matxsr_diff_free, c = 'C1', label = 'matxsr free')
    plt.semilogx(E,cxs_diff, c = 'C0',label = 'CXS')
    plt.axhline(y=0.0, color='r', linestyle='-')
    plt.yscale('symlog')
    plt.grid()
    plt.xlabel('energy')
    plt.ylabel('relative difference for '+filename[iso])
    plt.title('(Njoy-MCNP)/MCNP')
    plt.legend()
    # plt.savefig('../5_9 meeting/test_matxsr_files/Plots/'+filename[iso]+'_differences.png', bbox_inches='tight')
    # plt.savefig('../5_9 meeting/test_matxsr_files/Plots/Distributed_sources/'+filename[iso]+'_differences.png', bbox_inches='tight')

    diff_numerical = (outp_matxsr[0][1]-outp_cxs[0][1])/outp_matxsr[0][1]
    plt.figure()
    plt.semilogx(E,diff_numerical)
    plt.axhline(y=0.0, color='r', linestyle='-')
    plt.yscale('symlog')
    plt.grid()
    plt.xlabel('energy')
    plt.ylabel('relative difference for '+filename[iso])
    plt.title('(matxsr-cxs)/matxsr')
    plt.legend()
    # plt.savefig('../5_9 meeting/test_matxsr_files/Plots/'+filename[iso]+'_rel_differences.png', bbox_inches='tight')
    # plt.savefig('../5_9 meeting/test_matxsr_files/Plots/Distributed_sources/'+filename[iso]+'_rel_differences.png', bbox_inches='tight')
    
    mean_matxsr = np.average(matxsr_diff_free)
    median_matxsr = np.median(matxsr_diff_free)
    max_matxsr = np.max(abs(matxsr_diff_free))
    good_bins_matxsr = (abs(matxsr_diff_free) < 5e-2).sum()/2 
    ok_bins_matxsr = (abs(matxsr_diff_free) < 10e-2).sum()/2 
    
    mean_cxs = np.average(cxs_diff)
    median_cxs = np.median(cxs_diff)
    max_cxs = np.max(abs(cxs_diff))
    good_bins_cxs = (abs(cxs_diff) < 5e-2).sum()/2 
    ok_bins_cxs = (abs(cxs_diff) < 10e-2).sum()/2 
    cf.write("{:25}".format(filename[iso]))
    cf.write('\n')
    cf.write("{:25}".format(' '))
    cf.write("{:40}".format('Max relative difference'))
    cf.write("{:25}".format(str(max_matxsr)))
    cf.write("{:25}".format(str(max_cxs)))
    
    cf.write('\n')
    cf.write("{:25}".format(' '))
    cf.write("{:40}".format('Mean relative difference'))
    cf.write("{:25}".format(str(mean_matxsr)))
    cf.write("{:25}".format(str(mean_cxs)))
    
    cf.write('\n')
    cf.write("{:25}".format(' '))
    cf.write("{:40}".format('Median relative difference'))
    cf.write("{:25}".format(str(median_matxsr)))
    cf.write("{:25}".format(str(median_cxs)))
    
    cf.write('\n')
    cf.write("{:25}".format(' '))
    cf.write("{:40}".format('Bins with relative difference < 5%'))
    cf.write("{:25}".format(str(good_bins_matxsr)))
    cf.write("{:25}".format(str(good_bins_cxs)))
    
    cf.write('\n')
    cf.write("{:25}".format(' '))
    cf.write("{:40}".format('Bins with relative difference < 10%'))
    cf.write("{:25}".format(str(ok_bins_matxsr)))
    cf.write("{:25}".format(str(ok_bins_cxs)))
    cf.write('\n')
              
cf.close()