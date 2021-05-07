from MD_Analyzer import MD_Analyzer as MDA

#Generate input files for different percentages
for i in [160, 400, 800, 1200, 1360]:
    a = MDA('Large_cell.data', columns=['ID', 'TYPE',  'X', 'Y', 'Z', 'CHARGE'])
    a.replace_atom_number(i, 2, 3)
    a.data['Types'][3]['mass'] = 91.224
    a.simulation_ID = str(int(i*100/8000))
    a.save_as_lammps_data()
    a.create_lammps_input('ffield.reax', type_of_simulation = 'anneal', equiliberation_duration = 100000, annealing_duration = 100000000)