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

chixs_fullpath = []
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/He3.cxs')
# chixs_fullpath.append('../output/testing/XMAS_172/N14_n172.csx')
# chixs_fullpath.append('../output/testing/XMAS_172/H1_n172.csx')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Ar40.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Cd106.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Cd108.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Cd110.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Cd111.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Cd112.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Cd113.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Cd114.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Cd116.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Cu63.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Mn55.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Ti48.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Zn64.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Zn67.cxs')
# chixs_fullpath.append('../../outputs/ENDF-B-VIII-0/172gxs/CXS/Zn68.cxs')
chixs_fullpath.append('Al27.cxs')
N_density = []
N_density.append(1.)
data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath, N_density)
print(data.keys())


source_def = {}
source_def["particle_type"] = 'neutron'
source_def["energy"] = 18


outp = Utils_Info.InfiniteMediumSpectrum(data, source_def, plot=True)

Utils_NjoySpectrumPlotter.Njoy_spectrum_plotter(outp, './')

# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/He3.txt')
# A = np.loadtxt('../output/testing/XMAS_172/mcnp_N14_flx_tally.txt')
# A = np.loadtxt('../output/testing/XMAS_172/mcnp_H1_flx_tally.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Ar40_low.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd106mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd108mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd110mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd111mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd112mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd113mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd114mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd116mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cu63mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Mn55mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Ti48mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Zn64mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Zn67mid.txt')
# A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Zn68mid.txt')
MCNP_filename = 'Al27XSout.txt'


MCNP_flux = index(MCNP_filename, ['1tally', '14', '24'])
    
# this loop finds the line corresponding to 1tally and 14
tally_14_line = 0                                         # initialize value
for i in range(0,len(MCNP_flux['1tally'])):               # for each line corresponding to 1tally
    for j in range(0,len(MCNP_flux['14'])):               # loop over all lines corresponding to 14
        if MCNP_flux['1tally'][i] == MCNP_flux['14'][j]:  # if the line numbers correspond
            tally_14_line = MCNP_flux['14'][j]            # set as the output value

# # this loop finds the line corresponding to 1tally and 24            
# tally_24_line = 0                                         # initialize value
# for i in range(0,len(MCNP_flux['1tally'])):               # for each line corresponding to 1tally
#     for j in range(0,len(MCNP_flux['24'])):               # loop over all lines corresponding to 24
#         if MCNP_flux['1tally'][i] == MCNP_flux['24'][j]:  # if the line numbers correspond
#             tally_24_line = MCNP_flux['24'][j]            # set as the output value
            
# build a dictionary variable containing the two line numbers found above
MCNP_dict = {'tally_14': tally_14_line}#, 'tally_24' : tally_24_line}

# assemble the MCNP flux starting 10 lines after the tally 14 line appears
MCNP_data = assemble_data(MCNP_filename, MCNP_dict['tally_14'],9, 1)
# assemble MCNP RXN rates starting 11 lines after the tally 24 line appears
#MCNP_RXN = assemble_data(MCNP_filename, MCNP_dict['tally_24'],11, 1)


# MCNP = np.zeros([len(EnergyBins)])
# MCNP[0] = 0 # need the step function to start at 0
# for i in range(0,len(MCNP_fluxes)):
#     if float(MCNP_fluxes[i][1]) != 0:
#         MCNP[i+1] = float(MCNP_RXN[i][1])/float(MCNP_fluxes[i][1])
#     else: 
#         MCNP[i+1]=0






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
neutron_spectrum = outp[0][1]
#test['Al27'+'norm'] = (F-neutron_spectrum)/F
#F = F/np.sum(F)

fig_n = plt.gcf().number
#print(fig_n)
fig = plt.figure(fig_n)
ax_list = fig.axes
print(ax_list)
ax_list[0].semilogy(E, F, label='mcnp')
ax_list[1].loglog(E, F, label='mcnp')
plt.legend()
plt.show()
plt.savefig('Zn70mid.png')