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
            f.write("#!/bin/bash\n" + "#SBATCH --account=def-ivanov\n" + "#SBATCH --mem=16G \n" + "#SBATCH --time=1-12:00\n" + "#SBATCH --output=Bash_"+file_name[:-3]+".log \n" + "#SBATCH --cpus-per-task=36\n\n" + "module load nixpkgs/16.09  intel/2016.4  openmpi/2.1.1 lammps/20170331\n")
            f.write("\n lmp_exec=lmp_icc_openmpi\n")
            f.write("lmp_input=" + "\"" + file_name + "\"\n")
            f.write("lmp_output=" + "\"" + file_name[:-3] + ".lammpslog\"\n")
            f.write("srun ${lmp_exec} < ${lmp_input} > ${lmp_output}\n")
            f.close()
def run_bash():
    for file_name in os.listdir("./"):
        if file_name.endswith(".sh"):
            os.system("sbatch " + file_name)
       
percentages = [2]#, 5, 7, 10, 12, 15, 17]
a = MDA('atoms.data', '', columns=['ID', 'TYPE', 'CHARGE', '', 'X', 'Y', 'Z'])
a.data['Types'][3] = {'mass' : 91.224, 'count' : 0, 'IDs' : []}
for percent in percentages:
    for num in range(2):
        a.LAMMPS_Data_file = 'dielectric.data'
        a.simulation_ID = 'percent_{}_number_{}'.format(percent, num)
        while (a.data['Types'][3]['count']/a.data['Types'][2]['count']) < percent/100:
            a.replace_atoms(3, choice(a.data['Types'][2]['IDs']))
        a.save_as_lammps_data()
        a.create_lammps_input_anneal('ffield.reax')
create_bash_files()
run_bash()
