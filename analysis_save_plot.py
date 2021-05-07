# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 08:35:32 2020
cat fluctuate_anneal_Large_cell17.lammpstrj fluctuate_anneal_Large_cell17_1.lammpstrj fluctuate_anneal_Large_cell17_2.lammpstrj fluctuate_anneal_Large_cell17_3.lammpstrj fluctuate_anneal_Large_cell17_4.lammpstrj fluctuate_anneal_Large_cell17_5.lammpstrj fluctuate_anneal_Large_cell17_6.lammpstrj fluctuate_anneal_Large_cell17_7.lammpstrj > large_cell17.lammpstrj

@author: sarashs
"""
from MD_Analyzer import MD_Analyzer as MDA
from MD_Analyzer import consistent_plot
from matplotlib import pyplot as plot
#from tqdm import tqdm
import pickle
from os import path
import numpy as np


#name_radius = [['fluctuate_0.lammpstrj',14, 0],['large_cell21.lammpstrj',14, 2],['large_cell51.lammpstrj',14, 5],['large_cell101.lammpstrj',14, 10],['large_cell151.lammpstrj',14, 15]]#, ['fluctuate_anneal_dielectricpercent_100.lammpstrj',15,100]]#,['fluctuate_anneal_dielectricpercent_100_number_0.lammpstrj',10, 100]]
#name_radius1 = [['fluctuate_7.lammpstrj',14, 2.8],['fluctuate_0.lammpstrj',14, 0],['fluctuate_10.lammpstrj',14, 4.6], ['fluctuate_17.lammpstrj',14,11.8]]#,['fluctuate_anneal_dielectricpercent_100.lammpstrj',15,100]]#,['fluctuate_anneal_dielectricpercent_100_number_0.lammpstrj',10, 100]]

#a = MDA('fluctuate_anneal_dielectricpercent_100.lammpstrj', 200000, columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
##a = MDA('fluctuate_anneal_dielectricpercent_0_number_0.lammpstrj', 200000, columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
#DMF = {} #dict name : (Zr_percentage, dielectric_data)
#
###Load pickle file
#with open("dielectric_data_large.pickle", 'rb') as pickle_file:
#    DMF = pickle.load(pickle_file)
####Angle analysis
#import time
#anglessiosi = {}
#anglezrosi = {}
#cutoff = 2.1
#masses = {1:15.999,2:91.224,3:91.224}
#if __name__ == '__main__':
#    for i in name_radius:
##        start = time.time()
#        a = MDA(i[0], 2000000, columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
#        for item in a.data['Types'].keys():
#            a.data['Types'][item]['mass'] = masses[item]
#        print(a.density())
#        a.neighbors(cutoff)
#        lists = a.angle_statistics(2, 1, 2)
#        anglessiosi[i[0]] = lists
#        lists = a.angle_statistics(3, 1, 2)
#        anglezrosi[i[0]] = lists
#        end = time.time()
#        print(end - start)
###Coordination analysis per sphere size
#sphere_dia = 32
#Zr_type = 3
#O_type = 1
#neighbor_distance = 2.9
#Zr_neighbors = {}
#import time
#
#start = time.time()
#for i in tqdm(name_radius[:-1]):
#    a = MDA(i[0], 200000, columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
#    ID_list = a.sphere_set(sphere_dia)
#    Zr_ID_list = {}
#    for j in ID_list:
#        if a.data['data'][j]['TYPE'] == Zr_type:
#            Zr_ID_list[j] = 0
#    for j in a.data['Types'][O_type]['IDs']:
#        for k in Zr_ID_list.keys():
#            distance = 0
#            for d in ['X','Y','Z']:
#                distance += (a.data['data'][j][d] - a.data['data'][k][d])**2
#            distance = distance ** 0.5
#            if distance <= neighbor_distance:
#                Zr_ID_list[k] += 1
#    Zr_neighbors[i[2]] = Zr_ID_list
#    
#end = time.time()
#print(end - start)
#    
#X = []
#Y = []
#for i in Zr_neighbors.keys():
#    X.append(i)
#    Y.append(sum([int(j==6)/len(Zr_neighbors[i]) for j in Zr_neighbors[i].values()]))
#plot.plot(X,Y, 'o')
# 6 coordination percentage
#2[0.03278688524590164,
#5 0.06870229007633588,
#10 0.19178082191780849,
#15 0.24561403508771928]
##best radius and percentage, the results were evauated and the name radius list was filled manually
#DMF = {}
#for i in tqdm(name_radius):    
#    a = MDA(i[0], 200000, columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
#    DMF[i[0]] = a.sphere_percentage(27)
#    print(DMF)
###
#
#DMF_list = a.Dipole_moment_fluctuation(0, 20000, 12700000, 300, 14)#, shape = 'cube')#19120000, 300, 8.7)
#a = MDA('fluctuate_anneal_dielectricpercent_0_number_0.lammpstrj', 200000, columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
#DMF_list1 = a.Dipole_moment_fluctuation(0, 20000, 17120000, 300, 14)#, shape = 'cube')#19120000, 300, 8.7)
#DMF[100] = DMF_list
##
#compute all of the results for plotting
#from multiprocessing import Pool cpu_count
#p = Pool(processes=6)
#file_name = "dielectric_data_large_27.pickle"
#DMF = {}
#if path.exists(file_name):
#  with open(file_name, 'rb') as dielectric_data:
#      DMF = pickle.load(dielectric_data)
#      dielectric_data.close()
#else:
#  with open(file_name, 'wb') as dielectric_data:
#      pickle.dump(DMF,dielectric_data)
#      dielectric_data.close()    
#for i in tqdm(name_radius):
#  if (i[2] in DMF.keys()):
#      continue
#  print(i[0])
#  a = MDA(i[0], 0, columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
#  DMF[i[2]] = a.Dipole_moment_fluctuation(0, 20000, 12000000, 300, 25)
#  print(f'i[0] is {DMF[i[2]][-10 :-1]}')
#  with open(file_name, 'wb') as dielectric_data:
#      pickle.dump(DMF,dielectric_data)
#      dielectric_data.close()       

#file_name = "dielectric_data_large32.pickle"
#DMF = {}
#if path.exists(file_name):
#  with open(file_name, 'rb') as dielectric_data:
#      DMF = pickle.load(dielectric_data)
#      dielectric_data.close()
####Plot and save
#x = [i/2000000 for i in range(1000000, 18000000, 20000)] 
#xerr = None
###y = [i/2000000 for i in range(1000000, 18000000, 20000)]
#X = [x for i in [0, 2, 5]]
#XERR = [xerr for i in [0, 2, 5, 10,15]]
#YERR = {}
#YERR ={i : [np.std(DMF[i][:j]) for j in range(850)] for i in [0, 2, 5]}
#YERR[10] = [np.std(DMF[10][:j]) for j in range(1200)]
#YERR[15] = [np.std(DMF[15][:j]) for j in range(1200)]
#file_name = "error_large32.pickle"
#with open(file_name, 'wb') as error_data:
#   pickle.dump(YERR,error_data)
#   dielectric_data.close()  
#X.append([i/2000000 for i in range(0, 24000000, 20000)]) 
#X.append([i/2000000 for i in range(0, 24000000, 20000)]) 
#Y = [DMF[i][:850] for i in [0, 2, 5]]
#Y.append(DMF[10])
#Y.append(DMF[15])  
#consistent_plot(X, Y,YERR,XERR, ['$SiO_2$']+[f'{i}% $ZrO_2$' for i in [2, 5, 10,15]], ['-','-.','--','-s','-o'], [1 for i in [0, 2, 5, 10,15]], '$Time (ns)$', '$\epsilon$', 'Dielectric constants at 300 K', 'Dielectricconstant_32')
##plot.hist(DMF_list) 
##coordination numbers
from scipy.optimize import curve_fit
#xdata = np.array([0.1875,	0.354271357,	0.52875,	0.579166667])
#ydata = np.array([sum(DMF[i][-5:-1])/4 for i in [2, 5, 10,15]])
def func(x, a, b):
    return a + b * x
#    return a + (b * x)**0.5
#xdata_fit = np.array([0.006249219,0.032659214,0.104986877,0.174144899])
#popt, pcov = curve_fit(func, xdata, ydata)
#consistent_plot([xdata, xdata_fit], [ydata, func(xdata_fit ,popt[0] ,popt[1]) ],[None, None],[None, None], ['Data']+['Curve fit'], ['o','--'], [1 for i in [0, 1]], 'Percentage of 6 coordinated Zr', 'Dielectric constant ($\epsilon$)', 'Dielectric constant vs CN', '4coordinationsvsdielectric')
xdata = np.array([2,5,10,15])
ydata = np.array([4.08125,	4.366834171,	4.71625,	4.921666667])
#def func(x, a, b):
#    return a + (b * x)**0.5
xdata_fit = np.array([0.02,0.03,0.05,0.07,0.08,0.10,0.11,0.12,0.14,0.15,0.17])
popt, pcov = curve_fit(func, xdata, ydata)
consistent_plot([np.array([0.02,0.05,0.10,0.15]), xdata_fit], [ydata, func(xdata_fit ,popt[0] ,100*popt[1]) ],[None, None],[None, None], ['Data']+['Fit'], ['o','--'], [1 for i in [0, 1]], 'x', 'Coordination Number (CN)', '', 'fig_4_coord')