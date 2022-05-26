# -*- coding: utf-8 -*-
"""
Created on Mon May 16 09:31:53 2022

@author: tdeguire
"""

#%% load modules
import re
import shutil

#%% names of files to be generated and the corresponding isotope ID to be placed in each file

neutron_file= ['n-026_Fe_056.endf', 'n-013_Al_027.endf', 'n-006_C_012.endf', 'n-048_Cd_106.endf',
               'n-048_Cd_108.endf', 'n-048_Cd_110.endf', 'n-048_Cd_111.endf', 'n-048_Cd_112.endf',
               'n-048_Cd_113.endf', 'n-048_Cd_114.endf', 'n-048_Cd_116.endf', 'n-024_Cr_052.endf',
               'n-001_H_001.endf', 'n-002_He_003.endf', 'n-007_N_014.endf', 'n-028_Ni_058.endf',
               'n-008_O_016.endf', 'n-014_Si_028.endf', 'n-018_Ar_036.endf', 'n-018_Ar_038.endf',
                'n-018_Ar_040.endf', 'n-029_Cu_063.endf', 'n-029_Cu_065.endf', 'n-012_Mg_026.endf',
                'n-025_Mn_055.endf', 'n-022_Ti_048.endf', 'n-030_Zn_064.endf', 'n-030_Zn_066.endf',
                'n-030_Zn_067.endf', 'n-030_Zn_068.endf', 'n-030_Zn_070.endf']
sab_file= []

output_directory= []
output_file_prefix= ['Fe56', 'Al27', 'C12', 'Cd106', 
                     'Cd108', 'Cd110', 'Cd111', 'Cd112', 
                     'Cd113', 'Cd114', 'Cd116', 'Cr52', 
                     'H1', 'He3', 'N14', 'Ni58',
                     'O16', 'Si28', 'Ar36', 'Ar38',
                     'Ar40', 'Cu63', 'Cu65', 'Mg26',
                     'Mn55', 'Ti48', 'Zn64', 'Zn66',
                     'Zn67', 'Zn68', 'Zn70']

#%%% assemble input files

## this loop is what builds each input file
for i in range(0,len(neutron_file)):                   # for each isotope ID
    FileName = 'X1_template.sh'           # set filename to template file
    
    Output_file = output_file_prefix[i]+'.sh'                 # set name of output file
    shutil.copyfile(FileName, Output_file)          # Copy the input file to the output file
    with open(Output_file, 'r+') as f:              # open the output file
        text = f.read()                             # read all of the text
        text = re.sub('Neutron_filename', neutron_file[i], text)    # replace the material with the isotope of interest
        text = re.sub('output_prefix', output_file_prefix[i], text)    # replace the material with the isotope of interest
        text = re.sub('userID', output_file_prefix[i], text)    # replace the material with the isotope of interest
        f.seek(0)                                   # set pointer to beginning of file
        f.write(text)                               # write to file
        f.truncate()                                # truncate the files size
        
        #### need to add huse
