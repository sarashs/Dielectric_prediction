# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 11:53:44 2020

@author: sarashs
"""
from random import randint, choice
from matplotlib import pyplot as plot
#import multiprocessing
import time
import numpy as np
from copy import deepcopy
from numba import njit
#from tqdm import tqdm

#params = {'axes.labelsize': 36}#,'axes.titlesize':24, 'font.size': 24, 'legend.fontsize': 24, 'xtick.labelsize': 24, 'ytick.labelsize': 24}
#plot.rcParams.update(params)

Kb = 1.380649e-23
ε0 = 8.8541878128e-12
e_charge = 1.602176634e-19
Angstrom = 1e-10
coef = e_charge**2 / (Angstrom)
atomic_mass = 1.66053906660e-24 #grams

class Neighbor(object):
    def __init__(self, neighbor_type, neighbor_id, neighbor_distance, vector):
        self.TYPE = neighbor_type
        self.ID = neighbor_id
        self.DISTANCE = neighbor_distance
        self.VECTOR = vector
class Atom(object):
    def __init__(self, ID, TYPE, CHARGE, X, Y, Z, Vx = 0, Vy = 0, Vz = 0):
        self.ID = ID 
        self.TYPE = TYPE
        self.CHARGE = CHARGE
        self.X = X
        self.Y = Y
        self.Z = Z
        self.Vx = Vx
        self.Vy = Vy 
        self.Vz = Vz
        self.NEIGHBORS = []

@njit
def _periodic_subtract(a, b, dimentions):
    """This function computes the subtract between vectors a and b (b-a) considering the priodic boundary condition"""
    c = np.array([0.0, 0.0 ,0.0])
    for i in range(3):
        if a[i] < b[i]:
            if b[i] - a[i] < dimentions[i]:
                c[i] = b[i] - a[i]
            else:
                c[i] = b[i] - a[i] - dimentions[i]
        else:
            if np.abs(b[i] - a[i]) < dimentions[i]:
                c[i] = b[i] - a[i]
            else:
                c[i] = b[i] - a[i] + dimentions[i]
    dist = np.linalg.norm(c)
    return [dist, c]

class simple_bunching(object):
    def __init__(self, dimensions, grid_res, sorted_list_of_atoms):
        pass
class MD_Analyzer(object):
    """This is a python class for analyzing LAMMPS outputs
    :attribute simulation_ID: string
    ...
    """
    def __init__(self, Trajectory_file_name, simulation_ID = 0, columns=['ID', 'TYPE', 'CHARGE', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz']):
        self.number_of_atoms = 0
        self.place_holder = {'Dimensions' : 0, 'Masses' : 0, 'Atoms' : 0}
        self.data = {'Dimensions' : [0, 0, 0], 'Types' : {1 : {'mass' : 0, 'count' : 0, 'IDs' : []}}, 'data' : []} #{1 : {'TYPE' : 0, 'CHARGE' : 0, 'X' : 0, 'Y' : 0, 'Z' : 0, 'NEIGHBORS' : []}} 
        """
        self.data['data'] is a dictionary where ID is keys and value is a dictionary of TYPE, CHARGE, X, Y, Z, and NEIGHBORS (which is the list of neighbors NEIGHNBORS = [neighbor_type])
        """
        self.columns = {j : i for i,j in enumerate(columns)}
        self.last_timestep_index = 0
        self.magnetic_moment = 0 #total magnetic moment of the device
        self.Trajectory_file_name = Trajectory_file_name
        self.read_data(time_step = -1)
        self.simulation_ID = str(simulation_ID)
        self.updated_lines = self.lines
        self.time_step = 0
        self.neighbor_distance = 0
    def neighbors(self, neighbor_distance):
        self._neighbor_distance = neighbor_distance
        ## for multiprocessing
        serial = True
        ##Serial processing
        if serial == True: 
            print("Serial processing")
            for item in self.data['data']:
                for j in range(item.ID, len(self.data['data'])):
                    a = [item.X, item.Y, item.Z] 
                    b = [self.data['data'][j].X, self.data['data'][j].Y, self.data['data'][j].Z]
                    vector = self.periodic_subtract(a,b)
                    dist = np.linalg.norm(vector)
                    #[dist, vector] = _periodic_subtract(np.array(a),np.array(b), np.array(self.data['Dimensions']))
                    if (dist < neighbor_distance) and (dist > 0):
                        self.data['data'][j].NEIGHBORS.append(Neighbor(item.TYPE ,item.ID, dist, -vector))  
                        item.NEIGHBORS.append(Neighbor(self.data['data'][j].TYPE ,self.data['data'][j].ID, dist, vector))  
        ##Paralel processing
        else:
            pass
        #    for item in self.data['data']:
        #        _multi_processor_neighbor(item, neighbor_distance, self.data['data'], self.data['Dimensions'])
    def density(self):
        total_mass = 0
        for item in self.data['Types'].keys():
            total_mass += self.data['Types'][item]['mass'] * self.data['Types'][item]['count']
        return (1e-6 * Angstrom**(-3) * atomic_mass) * total_mass / (self.data['Dimensions'][0] * self.data['Dimensions'][1] * self.data['Dimensions'][2]) #gr/cm^3
    def angle_statistics(self, type1, type2, type3):
        print("make sure to run neighbor() first.")
        angle_list = []
        for item1 in self.data['data']:
            if (item1.TYPE == type1):
                for item2 in item1.NEIGHBORS:
                    for item3 in self.data['data'][item2.ID-1].NEIGHBORS:
                        if (item2.TYPE == type2) and (item3.TYPE == type3) and (item3.ID != item1.ID):
                            angle_list.append(np.arccos(np.dot(-item2.VECTOR, item3.VECTOR)/(np.linalg.norm(item2.VECTOR) * np.linalg.norm(item3.VECTOR)))* 180 / np.pi)
        return angle_list
    def read_data(self, time_step = -1):
        # By default the last timestep will be read
        if self.Trajectory_file_name.endswith('.lammpstrj'):
            self.LAMMPS_Data_file = self.Trajectory_file_name.replace('.lammpstrj','.data')
            self.file = open(self.Trajectory_file_name, 'r')
            self.lines = self.file.readlines()
            last_timestep = 0
            for i in range(len(self.lines)):
                if "ITEM: TIMESTEP" in self.lines[i]:
                    if time_step == -1:
                        if int(self.lines[i+1].replace(' ','')) > last_timestep:
                            last_timestep = int(self.lines[i+1].replace(' ',''))
                            self.last_timestep_index = i
                    else:
                        if int(self.lines[i+1].replace(' ','')) == time_step:
                            last_timestep = int(self.lines[i+1].replace(' ',''))
                            self.last_timestep_index = i  
                    self.time_step = int(self.lines[self.last_timestep_index + 1].replace(' ',''))
                    self.number_of_atoms = int(self.lines[self.last_timestep_index + 3].replace(' ',''))
                    self.data['data'] = [0] * self.number_of_atoms
                    x_dim = self.lines[self.last_timestep_index + 5].split()
                    y_dim = self.lines[self.last_timestep_index + 6].split() 
                    z_dim = self.lines[self.last_timestep_index + 7].split() 
                    self.data['Dimensions'] = [float(x_dim[1]) - float(x_dim[0]), float(y_dim[1]) - float(y_dim[0]), float(z_dim[1]) - float(z_dim[0])]
            for i in range(self.last_timestep_index + 9, self.last_timestep_index + 9 + self.number_of_atoms):
                temp = self.lines[i].split()
                if int(temp[self.columns['TYPE']]) in list(self.data['Types'].keys()):
                    self.data['Types'][int(temp[self.columns['TYPE']])]['count'] += 1
                    self.data['Types'][int(temp[self.columns['TYPE']])]['IDs'].append(int(temp[self.columns['ID']]))
                else:
                    self.data['Types'][int(temp[self.columns['TYPE']])] = {'mass' : 0, 'count' : 1, 'IDs' : [int(temp[self.columns['ID']])]}
                ID = int(temp[self.columns['ID']])
                TYPE = int(temp[self.columns['TYPE']])
                CHARGE = float(temp[self.columns['CHARGE']])
                X = float(temp[self.columns['X']]) - float(x_dim[0])
                Y = float(temp[self.columns['Y']]) - float(y_dim[0])
                Z = float(temp[self.columns['Z']]) - float(z_dim[0])
                self.data['data'][ID-1] = Atom(ID, TYPE, CHARGE, X, Y, Z)
            self.file.close()
                
        if self.Trajectory_file_name.endswith('.data'):
            self.LAMMPS_Data_file = self.Trajectory_file_name
            self.file = open(self.Trajectory_file_name, 'r')
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
                            self.data['data'][ID_index-1] = Atom(int(temp[ID_index]), int(temp[type_index]), float(temp[charge_index]), float(temp[x_index]), float(temp[y_index]), float(temp[z_index]))
                        else:
                            ID = j + 1 - self.place_holder['Atoms']
                            self.data['Types'][int(temp[type_index])]['IDs'].append(ID)
                            self.data['Types'][int(temp[type_index])]['count'] += 1
                            self.data['data'][ID-1] = Atom(ID, int(temp[type_index]), float(temp[charge_index]), float(temp[x_index]), float(temp[y_index]), float(temp[z_index]))
            self.file.close()
    def recenter(self, inputs):
        output = inputs
        for i,j in enumerate(output):
            output[i] -= self.data['Dimensions'][i]/2
        return output
    def periodic_subtract(self, a, b):
        """This function computes the subtract between vectors a and b (b-a) considering the priodic boundary condition"""
        c = np.array([0.0, 0.0 ,0.0])
        for i in range(3):
            if a[i] < b[i]:
                if b[i] - a[i] < self.data['Dimensions'][i]:
                    c[i] = b[i] - a[i]
                else:
                    c[i] = b[i] - a[i] - self.data['Dimensions'][i]
            else:
                if abs(b[i] - a[i]) < self.data['Dimensions'][i]:
                    c[i] = b[i] - a[i]
                else:
                    c[i] = b[i] - a[i] + self.data['Dimensions'][i]
        return c
    def sphere_percentage(self,radius):
        ## molar percent calculator
        num_zr = 0
        num_si = 0
        ##
        for i in range(1, self.number_of_atoms + 1):
            if np.linalg.norm(self.recenter([self.data['data'][i-1].X, self.data['data'][i-1].Y, self.data['data'][i-1].Z])) < radius:
                ## molar percent calculator
                if self.data['data'][i-1].TYPE == 3:
                    num_zr += 1
                elif self.data['data'][i-1].TYPE == 2:
                    num_si += 1
                ##
        return num_zr/(num_zr + num_si)
    def sphere_set(self,radius):
        ID_list = []
        for i in range(1, self.number_of_atoms + 1):
            if np.linalg.norm(self.recenter([self.data['data'][i-1].X, self.data['data'][i-1].Y, self.data['data'][i-1].Z])) < radius:
                ID_list.append(i)
        return ID_list
    def cube_set(self,from_the_edge):
        pass
    def cube_percentage(self,from_the_edge):
        """Returns the set of atoms which are at least from_the_edge distance from the edge"""
        ## molar percent calculator
        num_zr = 0
        num_si = 0
        ##
        for i in range(1, self.number_of_atoms + 1):
            if (self.data['data'][i-1].X - from_the_edge > 0) and (self.data['data'][i-1].X + from_the_edge < self.data['Dimensions'][0]):
                if (self.data['data'][i-1].Y - from_the_edge > 0) and (self.data['data'][i-1].Y + from_the_edge < self.data['Dimensions'][1]):
                    if (self.data['data'][i-1].Z - from_the_edge > 0) and (self.data['data'][i-1].Z + from_the_edge < self.data['Dimensions'][2]):
                        ## molar percent calculator
                        if self.data['data'][i-1].TYPE == 3:
                            num_zr += 1
                        elif self.data['data'][i-1].TYPE == 2:
                            num_si += 1
                        ##
        return num_zr/(num_zr + num_si)
    def Dipole_moment(self, time_step, ID_list):
        MF = 0
        temp = [0, 0, 0]
        self.read_data(time_step)
        for i in ID_list: #range(1, self.number_of_atoms + 1):
            atom = self.recenter([self.data['data'][i-1].X, self.data['data'][i-1].Y, self.data['data'][i-1].Z]) 
            #if self.data['data'][i]['ID'] in ID_list:
            temp = [k*self.data['data'][i-1].CHARGE + temp[j] for j,k in enumerate(atom)]
        #MF = np.linalg.norm(temp)#/number_of_molecules_in_sphere
        return np.array(temp)
    def Dipole_moment_fluctuation(self, init_timestep, span_timestep, final_timestep, T, radius, shape = 'sphere'):
        if 'cube' in shape:
            ID_list = self.cube_set(radius)
        else:
            ID_list = self.sphere_set(radius)
        start =  time.time()
        current_timestep = init_timestep
        DMF = self.Dipole_moment(current_timestep, ID_list)
        #DMF = MF
        DMF2 = DMF*DMF
        DMF_list = []
        MF_list = []
        num = 1
        c = coef/(ε0 * Kb*T)
        while current_timestep < final_timestep:
            #print(current_timestep)
            num += 1
            current_timestep += span_timestep
            MF = self.Dipole_moment(current_timestep, ID_list)
            DMF += MF
            DMF2 += MF*MF
            beta = c * np.sum(DMF2 / (num) - (DMF*DMF )/ (num**2))  
            if 'cube' in shape:
                beta = beta /(3 * (self.data['Dimensions'][0] - 2 * radius) * (self.data['Dimensions'][1] - 2 * radius) * (self.data['Dimensions'][2] - 2 * radius))
            else:
                beta = beta /(4 * np.pi * radius**3)
            DMF_list.append(beta + 1)#((1 + 3 * beta) + ((1 + 3 * beta)**2 + 8)**0.5)/4)
            MF_list.append(sum(MF))
        end = time.time()
        print(f'it took {end - start} seconds.')
        return DMF_list
    def save_as_lammps_data(self, time_step = -1):
        # By default the last timestep will be read
        if time_step != -1:
            self.read_data(time_step)
        file = open(self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.data', 'w')
        file.write('# System description #######################\n' + '#\n\n' + str(self.number_of_atoms) + ' atoms\n' + str(len(list(self.data['Types'].keys()))) + ' atom types\n')
        file.write('0 ' + str(self.data['Dimensions'][0]) + ' xlo xhi\n' + '0 ' + str(self.data['Dimensions'][1]) + ' ylo yhi\n' + '0 ' + str(self.data['Dimensions'][2]) + ' zlo zhi\n')
        file.write('#\n' + '# for a crystal:\n' + '# lx=a;  ly2+xy2=b2;  lz2+xz2+yz2=c2\n' + '# xz=c*cos(beta);  xy=b*cos(gamma)\n' + '# xy*xz+ly*yz=b*c*cos(alpha)\n' + '#\n\n')
        file.write('# Elements #################################\n\n' + 'Masses\n\n')
        for i in self.data['Types'].keys():
            if  self.data['Types'][i]['count'] > 0:
                file.write(str(i) + ' ' + str(self.data['Types'][i]['mass']) + '\n')
        file.write('\n' + 'Atoms\n' + '# number types charges\n')
        key_list = self.data['data']
        key_list.sort()
        for i in key_list:
            file.write(str(i).ljust(6) + str(self.data['data'][i-1].TYPE).ljust(4) + str(round(self.data['data'][i-1].CHARGE, 8)).ljust(12) + str(round(self.data['data'][i-1].X, 8)).ljust(12) + str(round(self.data['data'][i-1].Y, 8)).ljust(12) + str(round(self.data['data'][i-1].Z, 8)) + '\n')
        file.close()
    def replace_atoms(self, type_final, ID):
        """
        Changes the type of an atom
        """
        current_type = self.data['data'][ID-1].TYPE
        self.data['data'][ID-1].TYPE= type_final
        self.data['Types'][current_type]['count'] -= 1
        self.data['Types'][current_type]['IDs'].remove(ID)
        self.data['Types'][type_final]['count'] += 1
        self.data['Types'][type_final]['IDs'].append(ID)
    def replace_atom_number(self, num_atm_replace, type_init, type_final):
        """Replace num_atm_replace atoms of type_init to type_final. It is done randomly"""
        assert (type_init in self.data['Types']), f"{type_init} is not an available data type." 
        if(type_final not in self.data['Types']):
            self.data['Types'][type_final] = {'mass' : 0, 'count' : 0, 'IDs' : []}
            print(f"Critical Warning: Type {type_final} was not available and was added to the types. Other associated features such as mass must be added manually.")
        for i in range(num_atm_replace):
             ID = choice(self.data['Types'][type_init]['IDs'])
             self.replace_atoms(type_final, ID)
    def create_lammps_input(self, Input_forcefield, type_of_simulation = 'anneal', **kwargs):
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
        if 'equiliberation_duration' in kwargs.keys():
            equiliberation_duration = kwargs['equiliberation_duration']
        else:
            equiliberation_duration = 200000
        s.write('reset_timestep	0\n')
        s.write('timestep 0.1\n')
        s.write('velocity all create 300 ' + str(randint(1, 500000)) + ' rot yes mom yes dist gaussian\n')
        s.write('fix MD1 all nve\n')
        s.write('fix 10 all temp/rescale 1 300.0 300.0 1.0 0.5\n')
        s.write('dump DUMP1 all custom 10000 ' + 'equilib_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')
        s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
        s.write('thermo 1000\n')
        s.write(f'run {equiliberation_duration}\n')
        s.write('unfix 10\n')
        s.write('unfix MD1\n')
        s.write('fix MD2 all npt temp 300 300 20 aniso 1.0 1.0 50.0\n')
        s.write(f'run {equiliberation_duration}\n')
        s.write('unfix MD2\n')
        s.write('undump DUMP1\n')
        ##### Anneal
        if type_of_simulation in 'anneal':
            if 'anneal_duration' in kwargs.keys():
                anneal_duration = kwargs['anneal_duration']
            else:
                anneal_duration = 900000
            s.write('reset_timestep	0\n')
            s.write('fix MD3 all npt temp 300 3000 20 aniso 1.0 1.0 100.0\n')
            s.write('dump DUMP3 all custom 10000 ' + 'anneal_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')
            s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
            s.write('thermo 1000\n')
            s.write('run 100000\n')
            s.write('unfix MD3\n')
            s.write('fix MD4 all npt temp 3000 300 20 aniso 1.0 1.0 100.0\n')
            s.write(f'run {anneal_duration}\n')
            s.write('unfix MD4\n')
            s.write('fix MD5 all npt temp 300 300 20 aniso 1.0 1.0 20.0\n')
            s.write('run 100000\n')
            s.write('unfix MD5\n')
            s.write('undump DUMP3\n')
        # fluctuate
        if type_of_simulation in 'fluctuation':
            if 'fluctuation_duration' in kwargs.keys():
                fluctuation_duration = kwargs['fluctuation_duration']
            else:
                fluctuation_duration = 4500000
            if 'fluctuation_thermo_duration' in kwargs.keys():
                fluctuation_thermo_duration = kwargs['fluctuation_thermo_duration']
            else:
                fluctuation_thermo_duration = 10000
            s.write('reset_timestep	0\n')
            s.write('timestep 0.5\n')
            s.write('restart 500000 ' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.restart\n')
            s.write('fix MD6 all nvt temp 300 300 20\n')
            s.write('dump DUMP4 all custom 20000 ' + 'fluctuate_' + self.LAMMPS_Data_file.replace('.data','') + self.simulation_ID + '.lammpstrj'+' id type x y z q #this size \n')        
            s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
            s.write(f'thermo {fluctuation_thermo_duration}\n')
            s.write(f'run {fluctuation_duration}\n')
            s.write('unfix MD6\n')
            s.write('undump DUMP4\n')
        s.close()
    def restart_fluctuate(self, time_step, **kwargs):
        """
        This function creates the lammps input file
        :param Input_forcefield:
        """
        s=open(self.Trajectory_file_name.replace('.lammpstrj','') + '.in','w')
        s.write('# 1.- Inizialization #######################\n')
        s.write('read_restart ' + self.Trajectory_file_name.replace('.lammpstrj','.') + 'restart.' + str(time_step) + '\n')

        #Forcefield params
        s.write('pair_style reax/c NULL\n')
        if 'Input_forcefield' in kwargs.keys():
            pass
        else:
            Input_forcefield = 'ffield.reax'
        s.write('pair_coeff * * ' + Input_forcefield + ' ' + 'O Si Zr\n')
        #calculate number of atom types

        s.write('\n'+'fix 99 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c\n')
        s.write('neighbor        2.0 bin\n')
        s.write('neigh_modify    every 10 check yes\n\n')
        s.write('## 4.- MD & relax parameters ################\n\n')
        # fluctuate
        if 'fluctuation_duration' in kwargs.keys():
            fluctuation_duration = kwargs['fluctuation_duration']
        else:
            fluctuation_duration = 4500000
        if 'fluctuation_thermo_duration' in kwargs.keys():
            fluctuation_thermo_duration = kwargs['fluctuation_thermo_duration']
        else:
            fluctuation_thermo_duration = 10000
        s.write('timestep 0.5\n')
        s.write('restart 500000 ' + self.Trajectory_file_name.replace('.lammpstrj','') + '.restart\n')
        s.write('fix MD6 all nvt temp 300 300 20.0\n')
        s.write('dump DUMP4 all custom 10000 ' + 'fluctuate_' + self.Trajectory_file_name +' id type x y z q #this size \n')
        s.write('thermo_style custom step etotal ke pe temp press pxx pyy pzz \n')
        s.write(f'thermo {fluctuation_thermo_duration}\n')
        s.write(f'run {fluctuation_duration}\n')
        s.write('unfix MD6\n')
        s.write('undump DUMP4\n')
        s.close()
def consistent_plot(X, Y, YERR, XERR, labels, formats, alphas, xlabel, ylabel, title, save_name, bar = "NO"):
    """The arguments are all lists!"""
    plot.figure(figsize=(16,12))
    FONT_SIZE = 36
    plot.figure(figsize=(16,12))
    plot.rc('font', size=FONT_SIZE)          # controls default text sizes
    plot.rc('axes', titlesize=FONT_SIZE)     # fontsize of the axes title
    plot.rc('axes', labelsize=FONT_SIZE)    # fontsize of the x and y labels
    plot.rc('xtick', labelsize=FONT_SIZE)    # fontsize of the tick labels
    plot.rc('ytick', labelsize=FONT_SIZE)    # fontsize of the tick labels
    plot.rc('legend', fontsize=32)    # legend fontsize
    plot.rc('figure', titlesize=FONT_SIZE)  # fontsize of the figure title
    for i, j in enumerate(labels):
        if len(formats) > 0:
            if "YES" in bar:
                plot.bar(X[i], Y[i], formats[i], alpha = alphas[i], markersize=5.0,linewidth=3.0,label = labels[i])
            else:
                plot.errorbar(X[i], Y[i], yerr = YERR[i], xerr = XERR[i], fmt = formats[i], alpha = alphas[i], markersize=15.0,linewidth=3.0,label = labels[i])
        else:
            if "YES" in bar:
                plot.bar(X[i], Y[i], alpha = alphas[i], markersize=5.0,linewidth=3.0,label = labels[i])
            else:
                plot.errorbar(X[i], Y[i], yerr = YERR[i], xerr = XERR[i], fmt = formats[i], alpha = alphas[i], markersize=15.0,linewidth=3.0,label = labels[i])
    plot.xlabel(xlabel)
    plot.ylabel(ylabel)
    #plot.xticks(fontsize = 36)
    #plot.yticks(fontsize = 36)
    plot.title(title)
    plot.legend()
    plot.grid()
    plot.savefig(save_name + '.pdf')
    plot.show()