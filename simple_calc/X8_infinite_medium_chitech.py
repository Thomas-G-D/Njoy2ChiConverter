# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 00:55:42 2022

@author: jean.ragusa with edits from tdeguire
"""
import numpy as np

import sys, os
sys.path.insert(0, os.path.realpath('../chitech_composite_processor'))
import Utils_ChiTechCombiner
import Utils_Info
import Utils_NjoySpectrumPlotter
import matplotlib.pyplot as plt
from collections import defaultdict
import X8_Utils_MATXSR as matxsr  ### 9 is not currently compatible
plt.close('all')




## These functions are used to help assemble the MCNP data
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
                Output_values[x][0] = lines[skip_number + linenumber +x*multiplier][0] # energy value
                Output_values[x][1] = lines[skip_number + linenumber +x*multiplier][1] # tally result
                Output_values[x][2] = lines[skip_number + linenumber +x*multiplier][2] # uncertainty
    return Output_values

def Make1DPlotSpectrum(G : int, sigma : np.array, bounds : np.array):
    energies = np.zeros(2*G)
    sigmas = np.zeros(2*G)

    k : int = 0
    for g in range(0, G):
        elo = bounds[g]
        ehi = bounds[g+1]

        ewidth = abs(elo-ehi)

        energies[k  ] = elo 
        energies[k+1] = ehi

        sigmas[k  ] = sigma[g]/ewidth
        sigmas[k+1] = sigma[g]/ewidth

        # sigmas[k  ] = sigma[g]
        # sigmas[k+1] = sigma[g]

        k += 2

    return energies, sigmas

distributed_source = False

## prepare a file to write statistical information to
if distributed_source is True:
    cf = open('../5_9 meeting/NJOY_files/Plots/Distributed_sources/Information.txt', 'w')
else:
    cf = open('../5_9 meeting/NJOY_files/Plots/Information.txt', 'w')
cf.write("{:25}".format('Isotope'))
cf.write("{:40}".format('Value'))
cf.write("{:25}".format('MATXSR'))
cf.write("{:25}".format('CXS'))
cf.write('\n')



if distributed_source is True:
    filename = ['Al27',  'C12',   'Cd106', 'Cd108', 
                'Cd110', 'Cd111', 'Cd112', 'Cd113',
                'Cd114', 'Cd116', 'Cr52',  'Fe56',  
                'H1',    'He3',   'N14',   'Ni58',
                'O16',   'Si28',  'Ar36',  'Ar38',
                'Ar40',  'Cu63',  'Cu65',  'Mg26',
                'Mn55',  'Ti48',  'Zn64',  'Zn66',
                'Zn67',  'Zn68',  'Zn70']
else:
    filename = ['Al27',  'C12',   'Cd106', 'Cd108',
                'Cd110', 'Cd111', 'Cd112', 'Cd113',
                'Cd114', 'Cd116', 'Cr52',  'Fe56',  
                'H1',    'He3',   'N14',   'Ni58',
                'O16',   'Si28']

filename = ['Al27']
for iso in range(len(filename)):
    
    print('Running for '+filename[iso])  # print out what isotope is being processed
    
    ###############################################
    chitech = {}
    chioutfile = open(filename[iso]+'.co',"r")
    line = chioutfile.readline()
    words = line.split()
    G=172
    chispectrum = np.zeros(G)
    for g in range(0,G):
        chispectrum[g] = float(words[g])

    chitech["chispectrum"] = chispectrum
    
    
    
    
    
    
    #%% This section is for CXS and MATXSR
    ## add a path to the CXS file
    chixs_fullpath = []
    chixs_fullpath.append('../5_9 meeting/NJOY_files/'+filename[iso]+'.cxs')
    
    ## add the information needed to assemble the information from the CXS file into a dictionary
    N_density = []
    N_density.append(1.)
    data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath, N_density) # dictionary of cxs information
    
    ## assemble information from matxsr file
    matxsr_data = matxsr.EditMatrix('../5_9 meeting/NJOY_files/tape41_'+filename[iso]) # dictionary of matxsr information
    
    ## build the source term
    source_def = {}                                # initialize
    source_def["particle_type"] = 'neutron'        # set particle type
    source_def["energy"] = 18.49                   # energy of source in MeV
    
    ## use the source information to calculate the infinite medium spectrum
    outp_cxs    = Utils_Info.InfiniteMediumSpectrum(data, source_def, distributed_source, plot=True)
    outp_matxsr = Utils_Info.InfiniteMediumSpectrum(matxsr_data, source_def, distributed_source, plot=False)
    #outp_matxsr_sab = Utils_Info.InfiniteMediumSpectrum(matxsr_data, source_def, distributed_source, plot=False)
    
    ## plot the infinite medium spectrum using the CXS data
    Utils_NjoySpectrumPlotter.Njoy_spectrum_plotter(outp_cxs, './')
    
    
    ############# addomg here
    chitech["group_bounds"] = np.zeros(len(matxsr_data['energy_bins']))
    for energy_bin in range(len(matxsr_data['energy_bins'])):
        chitech["group_bounds"][energy_bin] = matxsr_data['energy_bins'][energy_bin]*1e-6
    e, sig = Make1DPlotSpectrum(G, np.flip(chitech["chispectrum"]), chitech["group_bounds"])
    
    
    ## add the MATXSR data to the plot of the infinite medium spectrum
    fig_n = plt.gcf().number
    fig = plt.figure(fig_n)
    ax_list = fig.axes
    #ax_list[0].semilogy(outp_matxsr_sab[0][0], outp_matxsr_sab[0][1], label = 'matxsr sab')
    #ax_list[1].loglog(outp_matxsr_sab[0][0], outp_matxsr_sab[0][1], label = 'matxsr sab')
    ax_list[0].semilogy(outp_matxsr[0][0], outp_matxsr[0][1], label = 'matxsr free')
    ax_list[1].loglog(outp_matxsr[0][0], outp_matxsr[0][1], label = 'matxsr free')
    ax_list[0].semilogy(e, sig, label = 'ChiTech')
    ax_list[1].loglog(e, sig, label = 'ChiTech')
    plt.legend()

    #%% This section is for MCNP
    
    ## file format and location slightly different depending on source type
    if distributed_source is True:
        
        MCNP_filename = 'MCNP/Distributded_src/'+filename[iso]+'.txt' # set the filename
        MCNP_flux = index(MCNP_filename, ['1tally', '104', '24'])     # find the lines in the file corresponding to the tallies
            
        ## this loop finds the line corresponding to tally results
        tally_14_line = 0                                          # initialize value
        for i in range(0,len(MCNP_flux['1tally'])):                # for each line corresponding to 1tally
            for j in range(0,len(MCNP_flux['104'])):               # loop over all lines corresponding to 14
                if MCNP_flux['1tally'][i] == MCNP_flux['104'][j]:  # if the line numbers correspond
                    tally_14_line = MCNP_flux['104'][j]            # set as the output value
              
        ## build a dictionary variable containing the two line numbers found above
        MCNP_dict = {'tally_104': tally_14_line}
        
        ## assemble the MCNP flux starting 10 lines after the tally 14 line appears
        MCNP_data = assemble_data(MCNP_filename, MCNP_dict['tally_104'],9, 1)
    else:
        
        MCNP_filename = 'MCNP/'+filename[iso]+'XSout.txt'         # set the filename
        MCNP_flux = index(MCNP_filename, ['1tally', '14', '24'])  # find the lines in the file corresponding to the tallies
            
        ## this loop finds the line corresponding to tally results
        tally_14_line = 0                                         # initialize value
        for i in range(0,len(MCNP_flux['1tally'])):               # for each line corresponding to 1tally
            for j in range(0,len(MCNP_flux['14'])):               # loop over all lines corresponding to 14
                if MCNP_flux['1tally'][i] == MCNP_flux['14'][j]:  # if the line numbers correspond
                    tally_14_line = MCNP_flux['14'][j]            # set as the output value
        
        ## build a dictionary variable containing the two line numbers found above
        MCNP_dict = {'tally_14': tally_14_line}
        
        ## assemble the MCNP flux starting 10 lines after the tally 14 line appears
        MCNP_data = assemble_data(MCNP_filename, MCNP_dict['tally_14'],9, 1)
    
    ## set the energies of the flux tallies and the tally results as their own variables
    mcnpE = MCNP_data[:,0]              # energy bins from MCNP
    mcnpF = MCNP_data[:,1] * 4e9        # MCNP Tally results multiplied by volume of the tally cell
    dE = np.diff(np.insert(mcnpE,0,0.)) # determine bin widths of MCNP bins
    spectrum = mcnpF/dE                 # calculate the MCNP spectrum variable
    
    ## need to convert MCNP data to a form that can be plotted as a histogram
    E = []
    F = []
    for g in range(MCNP_data.shape[0]-1):
        E += [mcnpE[g], mcnpE[g+1]]          # energy variable for histogram
        F += [spectrum[g+1], spectrum[g+1]]  # spectrum variable for histogram
    ## convert lists to arrays
    E = np.asarray(E)
    F = np.asarray(F)
    
    #%% this section is for plotting
    ## add MCNP information to the infinite medium plot and save it
    fig_n = plt.gcf().number
    fig = plt.figure(fig_n)
    ax_list = fig.axes
    ax_list[0].semilogy(E, F, label='mcnp')
    ax_list[1].loglog(E, F, label='mcnp')
    plt.legend()
    plt.suptitle('Neutron spectrum for '+filename[iso])
    # if distributed_source is True:
    #     plt.savefig('../5_9 meeting/NJOY_files/Plots/Distributed_sources/'+filename[iso]+'_spectrum.png', bbox_inches='tight')
    # else:
    #     plt.savefig('../5_9 meeting/NJOY_files/Plots/'+filename[iso]+'_spectrum.png', bbox_inches='tight')
    
    ## plot the relative differences of CXS and MATXSR to MCNP
    #matxsr_diff_sab = (outp_matxsr_sab[0][1]-F)/F
    matxsr_diff_free = (outp_matxsr[0][1]-F)/F  # relative difference of MATXSR to MCNP
    cxs_diff = (outp_cxs[0][1]-F)/F             # relative difference of CXS to MCNP
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
    # if distributed_source is True:
    #     plt.savefig('../5_9 meeting/NJOY_files/Plots/Distributed_sources/'+filename[iso]+'_differences.png', bbox_inches='tight')
    # else:
    #     plt.savefig('../5_9 meeting/NJOY_files/Plots/'+filename[iso]+'_differences.png', bbox_inches='tight')
        
    ## plot the relative differences of CXS to MATXSR
    diff_njoy = (outp_matxsr[0][1]-outp_cxs[0][1])/outp_matxsr[0][1]
    plt.figure()
    plt.semilogx(E,diff_njoy)
    plt.axhline(y=0.0, color='r', linestyle='-')
    plt.yscale('symlog')
    plt.grid()
    plt.xlabel('energy')
    plt.ylabel('relative difference for '+filename[iso])
    plt.title('(matxsr-cxs)/matxsr')
    plt.legend()
    # if distributed_source is True:
    #     plt.savefig('../5_9 meeting/NJOY_files/Plots/Distributed_sources/'+filename[iso]+'_rel_differences.png', bbox_inches='tight')
    # else:
    #     plt.savefig('../5_9 meeting/NJOY_files/Plots/'+filename[iso]+'_rel_differences.png', bbox_inches='tight')
    
    ######################################
    mcnpF = mcnpF[1:]
    compspec = np.clip(mcnpF,1.0e-10,np.max(mcnpF))
    plt.figure()
    plt.plot(np.flip(chitech["chispectrum"])/compspec, label=r"$\frac{serpent}{MCNP}$")
    plt.plot(outp_matxsr[3]/compspec, label=r"$\frac{matxsr}{MCNP}$")
    plt.axhline(y=1.0, color='r', linestyle='-')
    plt.legend()
    plt.xlabel('Bin')
    plt.ylabel('Ratio of deterministic flux to MCNP flux')
    plt.title('Comparison of results from Serpent2 and MATXSR for Al27')
    #plt.plot()
    
    #%% this section determines statistical information about the results and writes it out to text files
    
    ## statistical infor for matxsr
    mean_matxsr         = np.average(matxsr_diff_free)             # determine the mean relative difference for MATXSR
    median_matxsr       = np.median(matxsr_diff_free)              # determine the median relative difference for MATXSR
    max_matxsr          = np.max(abs(matxsr_diff_free))            # determine the max relative difference for MATXSR
    matxsr_5_diff_bins  = (abs((outp_matxsr[3]-compspec)/compspec) < 0.05).sum()#/2   # determine the number of bins with relative difference of less than 5%
    matxsr_10_diff_bins = (abs((outp_matxsr[3]-compspec)/compspec) < 0.1).sum()#/2  # determine the number of bins with relative difference of less than 10%
    # note, need to divide by 2 because the variables are in histogram for (i.e. bins will be counted twice once for each bound)
    
    ## statistical info for cxs
    mean_cxs         = np.average(cxs_diff)             # determine the mean relative difference for CXS
    median_cxs       = np.median(cxs_diff)              # determine the median relative difference for CXS
    max_cxs          = np.max(abs(cxs_diff))            # determine the max relative difference for CXS
    cxs_5_diff_bins  = (abs(cxs_diff) < 5e-2).sum()/2   # determine the number of bins with relative difference of less than 5%
    cxs_10_diff_bins = (abs(cxs_diff) < 10e-2).sum()/2  # determine the number of bins with relative difference of less than 10%
    
    ## write out the information to the text file
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
    cf.write("{:25}".format(str(matxsr_5_diff_bins)))
    cf.write("{:25}".format(str(cxs_5_diff_bins)))
    
    cf.write('\n')
    cf.write("{:25}".format(' '))
    cf.write("{:40}".format('Bins with relative difference < 10%'))
    cf.write("{:25}".format(str(matxsr_10_diff_bins)))
    cf.write("{:25}".format(str(cxs_10_diff_bins)))
    cf.write('\n')
              
cf.close()