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
import X3_Utils_MATXSR as matxsr
plt.close('all')

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
chixs_fullpath.append('Fe56.cxs')
N_density = []
N_density.append(1.)
data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath, N_density)
print(data.keys())

data2 = matxsr.EditMatrix('tape41')
data['transfer_matrices'] = data2

source_def = {}
source_def["particle_type"] = 'neutron'
source_def["energy"] = 18


outp = Utils_Info.InfiniteMediumSpectrum(data, source_def, plot=True)

Utils_NjoySpectrumPlotter.Njoy_spectrum_plotter(outp, './')

# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/He3.txt')
# # A = np.loadtxt('../output/testing/XMAS_172/mcnp_N14_flx_tally.txt')
# # A = np.loadtxt('../output/testing/XMAS_172/mcnp_H1_flx_tally.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Ar40_low.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd106mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd108mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd110mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd111mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd112mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd113mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd114mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cd116mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Cu63mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Mn55mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Ti48mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Zn64mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Zn67mid.txt')
# # A = np.loadtxt('../../outputs/ENDF-B-VIII-0/172gxs/MCNP/Zn68mid.txt')
# A = np.loadtxt('H1.txt')

# test = {}

# mcnpE = A[:,0]
# mcnpF = A[:,1] * 4e9 #4e9

# dE = np.diff(np.insert(mcnpE,0,0.))
# spectrum = mcnpF/dE

# E = []
# F = []
# for g in range(A.shape[0]-1):
#     E += [mcnpE[g], mcnpE[g+1]]
#     F += [spectrum[g+1], spectrum[g+1]]

# E = np.asarray(E)
# F = np.asarray(F)
# neutron_spectrum = outp[0][1]
# test['Al27'+'norm'] = (F-neutron_spectrum)/F
# #F = F/np.sum(F)

# fig_n = plt.gcf().number
# #print(fig_n)
# fig = plt.figure(fig_n)
# ax_list = fig.axes
# print(ax_list)
# ax_list[0].semilogy(E, F, label='mcnp')
# ax_list[1].loglog(E, F, label='mcnp')
# plt.legend()
# plt.show()
# plt.savefig('Zn70mid.png')