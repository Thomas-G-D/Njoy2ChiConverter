# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:50:57 2022

@author: tdeguire
"""

import Utils_XSWriter
import Utils_ChiTechCombiner
import X11_Utils_MATXSR_w_gamma as matxsr
filename = ['Al27',  'C12',   'Cd106', 'Cd108', 
                'Cd110', 'Cd111', 'Cd112', 'Cd113',
                'Cd114', 'Cd116', 'Cr52',  'Fe56',  
                'H1',    'He3',   'N14',   'Ni58',
                'O16',   'Si28',  'Ar36',  'Ar38',
                'Ar40',  'Cu63',  'Cu65',  'Mg26',
                'Mn55',  'Ti48',  'Zn64',  'Zn66',
                'Zn67',  'Zn68',  'Zn70', 'Cf252']


filename = ['graphite']
for iso in range(len(filename)):
    
    print('Running for '+filename[iso])  # print out what isotope is being processed
    #%% This section is for CXS and MATXSR
    ## add a path to the CXS file
    #chi_full_path = []
    chi_full_path = '../5_9 meeting/NJOY_files/matxsr_xs/'+filename[iso]+'.cxs'
    
    ## add the information needed to assemble the information from the CXS file into a dictionary
    N_density = []
    N_density.append(1.)
    #data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath, N_density) # dictionary of cxs information
    
    ## assemble information from matxsr file
    matxsr_data, problem_description = matxsr.EditMatrix('../5_9 meeting/NJOY_files/tape41_'+filename[iso]) # dictionary of matxsr information
    # problem_description = {}
    # problem_description['G_g'] = 0
    # problem_description['G_n'] = 172
    # problem_description['problem_type'] = 'Neutron only'
    # problem_description['isotope'] = filename[iso]
    #chi_full_path = 'Cf252.cxs'
    Utils_XSWriter.WriteChiTechFile(matxsr_data, chi_full_path, problem_description)
