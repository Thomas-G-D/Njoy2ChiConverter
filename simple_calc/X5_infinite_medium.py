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

filename = 'Ar36'

chixs_fullpath = []

#chixs_fullpath.append('C:/Users/tdeguire/Desktop/5_9 meeting/test_matxsr_files/Fe56.cxs')
chixs_fullpath.append('../5_9 meeting/test_matxsr_files/'+filename+'.cxs')
N_density = []
N_density.append(1.)
data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath, N_density)
print(data.keys())


source_def = {}
source_def["particle_type"] = 'neutron'
source_def["energy"] = 18.49
matxsr_data = matxsr.EditMatrix('../5_9 meeting/test_matxsr_files/tape41_'+filename)

outp_cxs    = Utils_Info.InfiniteMediumSpectrum(data, source_def, plot=True)
outp_matxsr = Utils_Info.InfiniteMediumSpectrum(matxsr_data, source_def, plot=False)
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
print(ax_list)
#ax_list[0].semilogy(outp_matxsr_sab[0][0], outp_matxsr_sab[0][1], label = 'matxsr sab')
#ax_list[1].loglog(outp_matxsr_sab[0][0], outp_matxsr_sab[0][1], label = 'matxsr sab')
ax_list[0].semilogy(outp_matxsr[0][0], outp_matxsr[0][1], label = 'matxsr free')
ax_list[1].loglog(outp_matxsr[0][0], outp_matxsr[0][1], label = 'matxsr free')
plt.legend()
#plt.show()
# plt.savefig('../5_9 meeting/test_matxsr_files/Plots/'+filename+'_spectrum.png')

# MCNP_filename = 'MCNP/'+filename+'XSout.txt'
# MCNP_flux = index(MCNP_filename, ['1tally', '14', '24'])
    
# # this loop finds the line corresponding to 1tally and 14
# tally_14_line = 0                                         # initialize value
# for i in range(0,len(MCNP_flux['1tally'])):               # for each line corresponding to 1tally
#     for j in range(0,len(MCNP_flux['14'])):               # loop over all lines corresponding to 14
#         if MCNP_flux['1tally'][i] == MCNP_flux['14'][j]:  # if the line numbers correspond
#             tally_14_line = MCNP_flux['14'][j]            # set as the output value
# MCNP_dict = {'tally_14': tally_14_line}#, 'tally_24' : tally_24_line}

# # assemble the MCNP flux starting 10 lines after the tally 14 line appears
# MCNP_data = assemble_data(MCNP_filename, MCNP_dict['tally_14'],9, 1)


MCNP_filename = 'MCNP/Distributded_src/'+filename+'.txt'


MCNP_flux = index(MCNP_filename, ['1tally', '104', '24'])
    
# this loop finds the line corresponding to 1tally and 14
tally_14_line = 0                                         # initialize value
for i in range(0,len(MCNP_flux['1tally'])):               # for each line corresponding to 1tally
    for j in range(0,len(MCNP_flux['104'])):               # loop over all lines corresponding to 14
        if MCNP_flux['1tally'][i] == MCNP_flux['104'][j]:  # if the line numbers correspond
            tally_14_line = MCNP_flux['104'][j]            # set as the output value

# this loop finds the line corresponding to 1tally and 24            
tally_24_line = 0                                         # initialize value
for i in range(0,len(MCNP_flux['1tally'])):               # for each line corresponding to 1tally
    for j in range(0,len(MCNP_flux['24'])):               # loop over all lines corresponding to 24
        if MCNP_flux['1tally'][i] == MCNP_flux['24'][j]:  # if the line numbers correspond
            tally_24_line = MCNP_flux['24'][j]            # set as the output value
            
# build a dictionary variable containing the two line numbers found above
MCNP_dict = {'tally_104': tally_14_line}#, 'tally_24' : tally_24_line}

# assemble the MCNP flux starting 10 lines after the tally 14 line appears
MCNP_data = assemble_data(MCNP_filename, MCNP_dict['tally_104'],9, 1)
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
#neutron_spectrum = outp_cxs[0][1]
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
#plt.show()
plt.savefig('../5_9 meeting/test_matxsr_files/Plots/Distributed_sources/'+filename+'_spectrum.png', bbox_inches='tight')

matxsr_diff_sab = (outp_matxsr_sab[0][1]-F)/outp_matxsr_sab[0][1]
matxsr_diff_free = (outp_matxsr[0][1]-F)/F#outp_matxsr[0][1]
cxs_diff = (outp_cxs[0][1]-F)/F#outp_cxs[0][1]
plt.figure()
plt.semilogx(E,matxsr_diff_free, c = 'C1', label = 'matxsr free')
plt.semilogx(E,cxs_diff, c = 'C0',label = 'CXS')
plt.axhline(y=0.0, color='r', linestyle='-')
plt.yscale('symlog')
plt.grid()
plt.xlabel('energy')
plt.ylabel('relative difference')
plt.title('(Njoy-MCNP)/MCNP')
plt.legend()
plt.savefig('../5_9 meeting/test_matxsr_files/Plots/Distributed_sources/'+filename+'_differences.png', bbox_inches='tight')

diff_numerical = (outp_matxsr[0][1]-outp_cxs[0][1])/outp_matxsr[0][1]
plt.figure()
plt.semilogx(E,diff_numerical)
plt.axhline(y=0.0, color='r', linestyle='-')
plt.yscale('symlog')
plt.grid()
plt.xlabel('energy')
plt.ylabel('relative difference')
plt.title('(matxsr-cxs)/matxsr')
plt.legend()
plt.savefig('../5_9 meeting/test_matxsr_files/Plots/Distributed_sources/'+filename+'_rel_differences.png', bbox_inches='tight')