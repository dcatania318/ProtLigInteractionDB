pdb_bind_file_path = '../data/PDBbind/'
pdb_bind_structure_names = []
with open(pdb_bind_file_path + 'INDEX_general_PL_name.2016', 'r') as pdb_bind_file:
    for line in pdb_bind_file.readlines():
        if not line.startswith('#'):
            line_split = line.strip('\n').split('  ')
            pdb_bind_structure_names.append(line_split[0])
print(len(pdb_bind_structure_names))

my_structure_names = []
with open('../data/PDB4/structure_names_pharmacophores', 'r') as my_file:
    for line in my_file.readlines():
        my_structure_names.append(line.strip('\n'))
print(len(my_structure_names))

crossover_structures = []
for structure_name in my_structure_names:
    if structure_name in pdb_bind_structure_names:
        crossover_structures.append(structure_name)
print(len(crossover_structures))

with open('../data/PDB4/crossover_structures', 'w') as crossover_structures_file:
    for structure_name in crossover_structures:
        crossover_structures_file.write(structure_name + '\n')
