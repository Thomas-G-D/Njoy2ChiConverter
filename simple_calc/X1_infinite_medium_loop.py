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
plt.close('all')
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

    A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/'+i+'.txt')
    mcnpE = A[:,0]
    mcnpF = A[:,1] * 4e9

    dE = np.diff(np.insert(mcnpE,0,0.))
    spectrum = mcnpF/dE

    E = []
    F = []
    for g in range(A.shape[0]-1):
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
    plt.savefig(i+'_difference.png')

