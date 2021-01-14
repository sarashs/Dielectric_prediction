# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:04:53 2020

@author: sarashs
"""

from MD_Analyzer import MD_Analyzer as MDA
import os
from random import choice
    
def create_bash_files(num_nodes = 2):
    for file_name in os.listdir("./"):
        if file_name.endswith(".in"):
            f = open(file_name[:-3] + f"_node{num_nodes}.sh", 'w')
            f.write("#!/bin/bash\n" + "#SBATCH --account=def-ivanov\n" + "#SBATCH --time=0-00:30\n" + f'#SBATCH --nodes={num_nodes}\n' + f'#SBATCH --ntasks={num_nodes*80}\n\n' + "module load intel/2019u3  intelmpi/2019u3  lammps/29Mar2019\n")
            f.write("\n" + "lmp_exec=lammps\n")
            f.write("lmp_input=" + "\"" + file_name + "\"\n")
            f.write("lmp_output=" + "\"" + file_name[:-3] + ".lammpslog\"\n")
            f.write("srun ${lmp_exec} < ${lmp_input} > ${lmp_output}\n")
            f.close()
def run_bash():
    for file_name in os.listdir("./"):
        if file_name.endswith(".sh"):
            os.system("sbatch " + file_name)
       
#Generate input files for different percentages
#for i in [160, 400, 800, 1200, 1360]:
#    a = MDA('Large_cell.data', columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
#    a.replace_atom_number(i, 2, 3)
#    a.data['Types'][3]['mass'] = 91.224
#    a.simulation_ID = str(int(i*100/8000))
#    a.save_as_lammps_data()
#    a.create_lammps_input('ffield.reax', type_of_simulation = 'anneal', equiliberation_duration = 100000, annealing_duration = 100000000)
#create_bash_files()
#run_bash()
########################
# find the fastest combination
a = MDA('Large_cell.data', columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
a.save_as_lammps_data()
a.create_lammps_input('ffield.reax', type_of_simulation = 'anneal', equiliberation_duration = 100000, annealing_duration = 500000)
for num_nodes in [1,2,3,4]:
    create_bash_files()
run_bash()            
#######################
#for file_name in os.listdir("./"):
#    if file_name.endswith('lammpstrj') and file_name.startswith('anneal'):            
#        print('here')
#        a = MDA(file_name, '', columns=['ID', 'TYPE', 'X', 'Y', 'Z', 'CHARGE']) 
#        a.restart_fluctuate(2500000)
#create_bash_files()
#run_bash()
#######################
#for file_name in os.listdir("./"):
#    if file_name.endswith(".lammpstrj"):            
#        a = MDA(file_name, '', columns=['ID', 'TYPE', 'X', 'Y', 'Z', 'CHARGE'])
#        a.data['Types'][1]['mass'] = 15.9994
#        a.data['Types'][2]['mass'] = 28.0860
#        a.data['Types'][3]['mass'] = 91.224
#        a.save_as_lammps_data()
#        a.create_lammps_input('ffield.reax', type_of_simulation = 'fluctuation', equiliberation_duration = 100000, fluctuation_duration = 100000000)
#create_bash_files()
#run_bash()
