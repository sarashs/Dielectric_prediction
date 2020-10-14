# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:04:53 2020

@author: sarashs
"""

from MD_Analyzer import MD_Analyzer as MDA
import os
from random import choice
    
def create_bash_files():
    for file_name in os.listdir("./"):
        if file_name.endswith(".in"):
            f = open(file_name[:-3] + ".sh", 'w')
            f.write("#!/bin/bash\n" + "#SBATCH --account=def-ivanov\n" + "#SBATCH --time=1-00:00\n" + "#SBATCH --nodes=2\n" + "#SBATCH --ntasks=160\n\n" + "module load intel/2019u3  intelmpi/2019u3  lammps/29Mar2019\n")
            f.write("\n" + "lmp_exec=lammps\n")
            f.write("lmp_input=" + "\"" + file_name + "\"\n")
            f.write("lmp_output=" + "\"" + file_name[:-3] + ".lammpslog\"\n")
            f.write("srun ${lmp_exec} < ${lmp_input} > ${lmp_output}\n")
            f.close()
def run_bash():
    for file_name in os.listdir("./"):
        if file_name.endswith(".sh"):
            os.system("sbatch " + file_name)
       
#percentages = [2, 5, 7, 10, 12, 15, 17]
#a = MDA('anneal_dielectricpercent_2_number_0.lammpstrj', '', columns=['ID', 'TYPE', 'X', 'Y', 'Z', 'CHARGE'])
#a.data['Types'][3] = {'mass' : 91.224, 'count' : 0, 'IDs' : []}
#for percent in percentages:
#    for num in range(10):
#        a.LAMMPS_Data_file = 'dielectric.data'
#        a.simulation_ID = 'percent_{}_number_{}'.format(percent, num)
#        while (a.data['Types'][3]['count']/a.data['Types'][2]['count']) < percent/100:
#            a.replace_atoms(3, choice(a.data['Types'][2]['IDs']))
#        a.save_as_lammps_data()
#        a.create_lammps_input('ffield.reax')
#create_bash_files()
#run_bash()
#####################################################
for file_name in os.listdir("./"):
    if file_name.endswith(".lammpstrj"):            
        a = MDA(file_name, '', columns=['ID', 'TYPE', 'X', 'Y', 'Z', 'CHARGE'])
        a.data['Types'][1]['mass'] = 15.9994
        a.data['Types'][2]['mass'] = 28.0860
        a.data['Types'][3]['mass'] = 91.224
        a.save_as_lammps_data()
        a.create_lammps_input('ffield.reax', type_of_simulation = 'fluctuation', equiliberation_duration = 100000, fluctuation_duration = 100000000)
create_bash_files()
#run_bash()