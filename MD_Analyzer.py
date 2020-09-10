# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 11:53:44 2020

@author: sarashs
"""
from random import randint

class MD_Analyzer(object):
    """This is a python class for analyzing LAMMPS outputs
    :attribute simulation_ID: string
    """
    def __init__(self, Trajectory_file_name, simulation_ID, columns=['ID', 'TYPE', 'CHARGE', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz']):
        self.number_of_atoms = 0
        self.place_holder = {'Dimensions' : 0, 'Masses' : 0, 'Atoms' : 0}
        self.data = {'Dimensions' : [0, 0, 0], 'Types' : {1 : {'mass' : 0, 'count' : 0, 'IDs' : []}}, 'data' : {1 : {'TYPE' : 0, 'CHARGE' : 0, 'X' : 0, 'Y' : 0, 'Z' : 0}}} 
        self.columns = {j : i for i,j in enumerate(columns)}
        if Trajectory_file_name.endswith('.data'):
            self.LAMMPS_Data_file = Trajectory_file_name
            self.file = open(Trajectory_file_name, 'r')
            self.lines = self.file.readlines()
            for j, i in enumerate(self.lines):
                if ' atoms' in i:
                    self.number_of_atoms = int(i.replace(' atoms', '').replace(' ', ''))
                    self.place_holder['Dimensions'] = j + 2
                if i.startswith('Masses'):
                    self.place_holder['Masses'] = j + 2
                if i.startswith('Atoms'):
                    self.place_holder['Atoms'] = j + 2
                if self.place_holder['Masses'] == 0 and self.place_holder['Dimensions'] > 0:
                    if i[0].isdigit():
                        dim_info = i.split()
                        if 'xhi' in dim_info:
                            self.data['Dimensions'][0] = float(dim_info[1])
                        elif 'yhi' in dim_info:
                            self.data['Dimensions'][1] = float(dim_info[1])
                        elif 'zhi' in dim_info:
                            self.data['Dimensions'][2] = float(dim_info[1])
                if self.place_holder['Atoms'] == 0 and self.place_holder['Masses'] > 0:
                    if i[0].isdigit():
                        mass_info = i.split()
                        self.data['Types'][int(mass_info[0])] = {'mass' : float(mass_info[1]), 'count' : 0, 'IDs' : []}
                if self.place_holder['Atoms'] > 0 :
                    if i[0].isdigit():
                        temp = i.strip().split()
                        type_index = self.columns['TYPE']
                        charge_index = self.columns['CHARGE']
                        x_index = self.columns['X']
                        y_index = self.columns['Y']
                        z_index = self.columns['Z']
                        if 'ID' in self.columns.keys():
                            ID_index = self.columns['ID']
                            self.data['Types'][int(temp[type_index])]['IDs'].append(int(temp[ID_index]))
                            self.data['Types'][int(temp[type_index])]['count'] += 1
                            self.data['data'][int(temp[ID_index])] = {'TYPE' : int(temp[type_index]), 'CHARGE' : float(temp[charge_index]), 'X' : float(temp[x_index]), 'Y' : float(temp[y_index]), 'Z' : float(temp[z_index])}
                        else:
                            ID = j + 1 - self.place_holder['Atoms']
                            self.data['Types'][int(temp[type_index])]['IDs'].append(ID)
                            self.data['Types'][int(temp[type_index])]['count'] += 1
                            self.data['data'][ID] = {'TYPE' : int(temp[type_index]), 'CHARGE' : float(temp[charge_index]), 'X' : float(temp[x_index]), 'Y' : float(temp[y_index]), 'Z' : float(temp[z_index])}

            self.file.close()
        elif Trajectory_file_name.endswith('.lammpstrj'):
            self.LAMMPS_Data_file = Trajectory_file_name.replace('.lammpstrj', '.data')
        self.simulation_ID = simulation_ID
        self.updated_lines = self.lines
    def Magnetic_fluctuation(self, number_of_time_steps):
        pass
    def save_as_lammps_data(self):
        file = open(self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.data', 'w')
        file.write('# System description #######################\n' + '#\n\n' + str(self.number_of_atoms) + ' atoms\n' + str(len(list(self.data['Types'].keys()))) + ' atom types\n')
        file.write('0 ' + str(self.data['Dimensions'][0]) + ' xlo xhi\n' + '0 ' + str(self.data['Dimensions'][1]) + ' ylo yhi\n' + '0 ' + str(self.data['Dimensions'][2]) + ' zlo zhi\n')
        file.write('#\n' + '# for a crystal:\n' + '# lx=a;  ly2+xy2=b2;  lz2+xz2+yz2=c2\n' + '# xz=c*cos(beta);  xy=b*cos(gamma)\n' + '# xy*xz+ly*yz=b*c*cos(alpha)\n' + '#\n\n')
        file.write('# Elements #################################\n\n' + 'Masses\n\n')
        for i in self.data['Types'].keys():
            if  self.data['Types'][i]['count'] > 0:
                file.write(str(i) + ' ' + str(self.data['Types'][i]['mass']) + '\n')
        file.write('\n' + 'Atoms\n' + '# number types charges\n')
        key_list = list(self.data['data'].keys())
        key_list.sort()
        for i in key_list:
            file.write(str(i).ljust(6) + str(self.data['data'][i]['TYPE']).ljust(4) + str(self.data['data'][i]['CHARGE']).ljust(12) + str(self.data['data'][i]['X']).ljust(12) + str(self.data['data'][i]['Y']).ljust(12) + str(self.data['data'][i]['Z']) + '\n')
        file.close()
    def replace_atoms(self, type_final, ID):
        current_type = self.data['data'][ID]['TYPE']
        self.data['data'][ID]['TYPE'] = type_final
        self.data['Types'][current_type]['count'] -= 1
        self.data['Types'][current_type]['IDs'].remove(ID)
        self.data['Types'][type_final]['count'] += 1
        self.data['Types'][type_final]['IDs'].append(ID)
    def create_lammps_input_anneal(self, Input_forcefield):
        """
        This function creates the lammps input file
        :param Input_forcefield:
        """
        s=open(self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.in','w')
        #for n in lists: 
        ######
        #s.write('log ' + self.LAMMPS_Data_file.replace('.data', '.log') + '\n')
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
        s.write('read_data   '+ self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.data'+'\n')
        s.write('\n# 3.- Force-Field ##########################\n\n')
        #Forcefield params
        s.write('pair_style reax/c NULL\n')
        s.write('pair_coeff * * ' + Input_forcefield + ' ' + 'O Si Zr\n')
        #calculate number of atom types

        s.write('\n'+'fix 99 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c\n')
        s.write('neighbor        2.0 bin\n')
        s.write('neigh_modify    every 10 check yes\n\n')
        s.write('## 4.- MD & relax parameters ################\n\n')
        ######
        s.write('dump DUMP2 all custom 10000 ' + 'init_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')
        s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
        s.write('thermo 1000\n')
        #This line decides whether or not the fix is going to be in the calculated energy
        s.write('min_style cg\n')
        s.write('minimize 1.0e-5 1.0e-6 2000 2000\n')
        s.write('undump DUMP2\n')
        ##### EQUILIBRATION
        s.write('reset_timestep	0\n')
        s.write('timestep 0.1\n')
        s.write('velocity all create 300 ' + str(randint(1, 500000)) + ' rot yes mom yes dist gaussian\n')
        s.write('fix MD1 all nve\n')
        s.write('fix 10 all temp/rescale 1 300.0 300.0 1.0 0.5\n')
        s.write('dump DUMP1 all custom 10000 ' + 'equilib_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')
        s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
        s.write('thermo 1000\n')
        s.write('run 200000\n')
        s.write('unfix 10\n')
        s.write('unfix MD1\n')
        s.write('fix MD2 all npt temp 300 300 20 aniso 1.0 1.0 50.0\n')
        s.write('run 200000\n')
        s.write('unfix MD2\n')
        s.write('undump DUMP1\n')
        ##### Anneal
        s.write('reset_timestep	0\n')
        s.write('fix MD3 all npt temp 300 3000 20 aniso 1.0 1.0 100.0\n')
        s.write('dump DUMP3 all custom 10000 ' + 'anneal_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')
        s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
        s.write('thermo 1000\n')
        s.write('run 100000\n')
        s.write('unfix MD3\n')
        s.write('fix MD4 all npt temp 3000 300 20 aniso 1.0 1.0 100.0\n')
        s.write('run 900000\n')
        s.write('unfix MD4\n')
        s.write('fix MD5 all npt temp 300 300 20 aniso 1.0 1.0 20.0\n')
        s.write('run 100000\n')
        s.write('unfix MD5\n')
        s.write('undump DUMP3\n')
        # fluctuate
        #s.write('reset_timestep	0\n')
        #s.write('fix MD6 all nvt temp 300 300 20.0\n')
        #s.write('dump DUMP4 all custom 5000 ' + 'fluctuate_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')
        #s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
        #s.write('thermo 1000\n')
        #s.write('run 10000000\n')
        #s.write('unfix MD6\n')
        #s.write('undump DUMP4\n')
        s.close()
    def consisten_plot(self):
        pass
