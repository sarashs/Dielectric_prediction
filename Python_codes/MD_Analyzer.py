# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 11:53:44 2020

@author: sarashs
"""
import os
from random import randint

class MD_Analyzer(object):
    """This is a python class for analyzing LAMMPS outputs
    :attribute simulation_ID: string
    """
    def __init__(self, Trajectory_file_name, simulation_ID, columns=['ID', 'TYPE', 'CHARGE', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz']):
        self.number_of_atoms = 0
        self.place_holder = {'Dimensions' : 0, 'Masses' : 0, 'Atoms' : 0}
        self.data = {'Dimensions' : (0, 0, 0), 'Types' : {1 : {'mass' : 0, 'count' : 0, 'IDs' : []}}, 'data' = [{1 : {'TYPE' : 0, 'CHARGE' : 0, 'X' : 0, 'Y' : 0, 'Z' : 0}}]} 
        self.columns = {i : j for i,j in enumerate(columns)}
        if Trajectory_file_name.endswith('.data'):
            self.LAMMPS_Data_file = Trajectory_file_name
            self.file = open(Trajectory_file_name, 'r')
            self.lines = self.file.readlines()
            for i, j in enumerate(self.lines):
                if i.endswith(' atoms'):
                    self.number_of_atoms = int(self.lines.replace(' atoms', '').replace(' ', ''))
                    self.place_holder['Dimensions'] = j + 2
                if i.startswith('Masses'):
                    self.place_holder['Masses'] = j + 2
                if i.startswith('Atoms'):
                    self.place_holder['Atoms'] = j + 2
                if self.place_holder['Masses'] < 0 and self.place_holder['Dimensions'] > 0:
                    if i[0].isdigit():
                        dim_info = i.split()
                        if 'xhi' in dim_info:
                            self.data['Dimensions'][0] = float(dim_info[0])
                        elif 'yhi' in dim_info:
                            self.data['Dimensions'][1] = float(dim_info[1])
                        elif 'zhi' in dim_info:
                            self.data['Dimensions'][2] = float(dim_info[2])
                if self.place_holder['Atoms'] < 0 and self.place_holder['Masses'] > 0:
                    if i[0].isdigit():
                        mass_info = i.split()
                        self.data['Types'][int(mass_info[0])] = {'mass' : float(mass_info[1]), 'count' : 0, 'IDs' : []}
                if self.place_holder['Atoms'] > 0 :
                    temp = i.strip().split()
                    self.data['data'] = temp
                    type_index = self.columns['TYPE']
                    charge_index = self.columns['CHARGE']
                    x_index = self.columns['X']
                    y_index = self.columns['Y']
                    z_index = self.columns['Z']
                    if 'ID' in self.columns.keys():
                        ID_index = self.columns['ID']
                        self.data['Types'][temp[type_index]]['IDs'].append(int(temp[ID_index]))
                        self.data['data'][int(temp[ID_index])] = {'TYPE' : int(temp[ID_index]), 'CHARGE' : float(temp[charge_index]), 'X' : float(temp[x_index]), 'Y' : float(temp[y_index]), 'Z' : float(temp[z_index])}
                    else:
                        self.data['Types'][temp[type_index]]['IDs'].append(j + 1 - self.place_holder['Atoms'])
            self.file.close()
        elif Trajectory_file_name.endswith('.lammpstrj'):
            self.LAMMPS_Data_file = Trajectory_file_name.replace('.lammpstrj', '.data')
        self.simulation_ID = simulation_ID
        self.updated_lines = self.lines
    def Magnetic_fluctuation(self, number_of_time_steps):
        pass
    def save_as_lammps_data(self):
        file = open(self.LAMMPS_Data_file, 'w')
        file.write('# System description #######################\n' + '#\n\n' + str(self.number_of_atoms) + ' atoms\n' + len(list(self.data['Types'].keys())) + ' atom types\n')
        file.write('0 ' + self.data['Dimensions'][0] + ' xlo xhi\n' + '0 ' + self.data['Dimensions'][1] + ' ylo yhi\n' + '0 ' + self.data['Dimensions'][2] + ' zlo zhi\n')
        file.write('#\n' + '# for a crystal:\n' + '# lx=a;  ly2+xy2=b2;  lz2+xz2+yz2=c2\n' + '# xz=c*cos(beta);  xy=b*cos(gamma)\n' + '# xy*xz+ly*yz=b*c*cos(alpha)\n' + '#\n\n')
        file.write('# Elements #################################\n\n' + 'Masses\n\n')
        for i in self.data['Types'].keys():
            file.write(i + ' ' + self.data['Types'][i]['mass'] + '\n')
        file.write('\n' + 'Atoms\n' + '# number types charges\n')
        for i in range(self.number_of_atoms):
            file.write(str(i).ljust(6) + str(self.data['data'][i]['type']).ljust(4) + str(self.data['data'][i]['CHARGE']).ljust(4) + str(self.data['data'][i]['X']).ljust(12) + str(self.data['data'][i]['Y']).ljust(12) + str(self.data['data'][i]['Z']).ljust(12))
        file.close()
    def replace_atoms(self, type_final, ID):
        current_type = self.data['data'][ID]['TYPE']
        self.data['data'][ID]['TYPE'] = type_final
        self.data['Types'][current_type]['count'] - = 1
        self.data['Types'][current_type]['IDs'].remove(ID)
        self.data['Types'][type_final]['count'] + = 1
        self.data['Types'][type_final]['IDs'].append(ID)
    def create_lammps_input_anneal(self, Input_forcefield):
        """
        This function creates the lammps input file
        :param Input_forcefield:
        """
        s=open(self.LAMMPS_Data_file.repace('.data','') + self.simulation_ID,'w')
        #for n in lists: 
        ######
        s.write('log ' + self.LAMMPS_Input_file.replace('.dat', '.log') + '\n')
        s.write('# 1.- Inizialization #######################\n')
        s.write('units real\n')
        s.write('  #mass = grams/mole\n')
        s.write('  #distance = Angstroms\n')
        s.write('  #time = femtoseconds\n')
        s.write('  #energy = kcal/mol\n')
        s.write('  #velocity = Angstroms/femtosecond\n')
        s.write('  #force = kcal/mol.Angstrom\n')
        s.write('  #torque = kcal/mole\n')
        s.write('  #temperature = degrees K\n')
        s.write('  #pressure = atmospheres (0.1013 GPa)\n')
        s.write('  #dynamic viscosity = Poise\n')
        s.write('  #charge = multiple of electron charge (+1.0 is a proton)\n')
        s.write('  #dipole = charge*Angstroms\n')
        s.write('  #electric field = volts/Angstrom\n')
        s.write('dimension 3\n')
        s.write('processors * * *\n')
        s.write('##\n')
        s.write('boundary p p p\n')
        s.write('atom_style charge\n\n# 2.- Atom definition ######################\n\n')
        s.write('atom_modify map hash\n')
        s.write('read_data   '+self.LAMMPS_Data_file+'\n')
        s.write('\n# 3.- Force-Field ##########################\n\n')
        #Forcefield params
        s.write('pair_style reax/c NULL\n')
        s.write('pair_coeff * * ' + Input_forcefield + 'H O Si Zr\n')
        #calculate number of atom types

        s.write('\n'+'fix 99 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c\n')
        s.write('neighbor        2.0 bin\n')
        s.write('neigh_modify    every 10 check yes\n\n')
        s.write('## 4.- MD & relax parameters ################\n\n')
        ######
        s.write('dump DUMP2 all custom 1000 ' + 'init_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')
        s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
        s.write('thermo 1000\n')
        #This line decides whether or not the fix is going to be in the calculated energy
        s.write('min_style cg\n')
        s.write('minimize 1.0e-5 1.0e-6 2000 2000\n')
        s.write('undump DUMP2\n')
        ##### EQUILIBRATION
        s.write('reset_timestep	0\n')
        s.write('timestep 0.1\n')
        s.write('velocity all create 300' + str(randint(1, 500000)) + 'dist gaussian\n')
        s.write('fix MD1 all nve\n')
        s.write('dump DUMP1 all custom 1000 ' + 'equilib_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')
        s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
        s.write('thermo 1000\n')
        s.write('run 100000\n')
        s.write('unfix MD1\n')
        s.write('fix MD2 all npt temp 300 300 20 aniso 1.0 1.0 100.0\n')
        s.write('run 100000\n')
        s.write('unfix MD2\n')
        s.write('undump DUMP1\n')
        ##### Anneal
        s.write('reset_timestep	0\n')
        s.write('fix MD3 all npt temp 300 1200 20 aniso 1.0 1.0 100.0\n')
        s.write('dump DUMP3 all custom 10000 ' + 'anneal_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')
        s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
        s.write('thermo 1000\n')
        s.write('run 100000\n')
        s.write('unfix MD3\n')
        s.write('fix MD4 all npt temp 1200 300 20 aniso 1.0 1.0 100.0\n')
        s.write('run 500000\n')
        s.write('unfix MD4\n')
        s.write('fix MD5 all npt temp 300 300 20 aniso 1.0 1.0 100.0\n')
        s.write('run 500000\n')
        s.write('unfix MD5\n')
        s.write('undump DUMP3\n')
        s.close()
    def consisten_plot(self):
        pass
def create_bash_files(self):
    for file_name in os.listdir("./"):
        if file_name.endswith(".in"):
            f = open(file_name[:-4] + ".sh", 'w')
            f.write("#!/bin/bash\n" + "#SBATCH --account=def-ivanov\n" + "#SBATCH --mem=16G \n" + "#SBATCH --time=5-00:00\n" + "#SBATCH --output=Bash_"+file_name[:-4]+".log \n" + "#SBATCH --cpus-per-task=36\n\n" + "module load nixpkgs/16.09  intel/2016.4  openmpi/2.1.1 lammps/20170331\n")
            f.write("\n lmp_exec=lmp_icc_openmpi\n")
            f.write("lmp_input=" + "\"" + file_name + ".in\"\n")
            f.write("lmp_output=" + "\"" + file_name[:-4] + ".lammpslog\"\n")
            f.write("srun ${lmp_exec} < ${lmp_input} > ${lmp_output}\n")
            f.close()
def run_bash(self):
    for file_name in os.listdir("./"):
        if file_name.endswith(".sh"):
            os.system("sbatch " + file_name)