import math
import multiprocessing as mp
from collections import defaultdict
import re

pdb_bind_file_path = '../data/PDBbind/'
pdb_bind_dict = defaultdict(list)
with open(pdb_bind_file_path + 'INDEX_general_PL_data.2016', 'r') as pdb_bind_file:
    for line in pdb_bind_file.readlines():
        if not line.startswith('#'):
            line_split = line.strip('\n').split('  ')
            match = re.findall(r'(Ki|Kd|IC50)|(\d*\.?\d+)|([a-z]*M)', line_split[4])
            if match:
                type = match[0][0]
                binding_affinity = float(match[1][1])
                measure = match[2][2]
                if (len(measure) > 0) and (measure[len(measure)-1] == 'M'):
                    if measure[0] == 'm':
                        binding_affinity *= 1e-3
                    elif measure[0] == 'u':
                        binding_affinity *= 1e-6
                    elif measure[0] == 'n':
                        binding_affinity *= 1e-9
                    elif measure[0] == 'p':
                        binding_affinity *= 1e-12
                    elif measure[0] == 'f':
                        binding_affinity *= 1e-15
                    elif measure[0] == 'M':
                        continue
                    else:
                        print(line_split[4])
                        print(type,binding_affinity)
                        print(measure)
                        break
                    binding_affinity = round(binding_affinity/1e-9, 2)
                    pdb_bind_dict[line_split[0]] = (type,binding_affinity)

pharmacophores_data = [
    # receptor_pharmacophore, ligand_pharmacophore, distance
    ['Hydrophobe', 'Hydrophobe', 4.5],
    ['Acceptor', 'Donor', 3.9],
    ['Donor', 'Acceptor', 3.9],
    ['PosIonizable', 'NegIonizable', 4.0],
    ['NegIonizable', 'PosIonizable', 4.0],
    ['Aromatic', 'Aromatic', 4.5],
    ['Aromatic', 'PosIonizable', 4.0],
    ['PosIonizable', 'Aromatic', 4.0]
]


class data:
    def __init__(self, string):
        self.id = string[0]
        self.pharmacophore = string[1]
        self.x = float(string[2])
        self.y = float(string[3])
        self.z = float(string[4])


def check_features(path, file_extension):
    hetatm_features = []
    with open(path + '/hetatm_' + file_extension + '.feats') as hetatm_feats_file:
        for line in hetatm_feats_file.readlines():
            line_split = line.strip('\n').split(',')
            hetatm_features.append(data(line_split))
    atom_features = []
    with open(path + '/atom_' + file_extension + '.feats') as atom_feats_file:
        for line in atom_feats_file.readlines():
            line_split = line.strip('\n').split(',')
            atom_features.append(data(line_split))
    pharmacophores_count = [0, 0, 0, 0, 0, 0, 0, 0]
    for atom_feature in atom_features:
        for hetatm_feature in hetatm_features:
            for i in range(0, 8):
                pharmacophore_data = pharmacophores_data[i]
                if atom_feature.pharmacophore == pharmacophore_data[0] and hetatm_feature.pharmacophore == pharmacophore_data[1]:
                    distance = math.sqrt(math.pow((atom_feature.x - hetatm_feature.x), 2) + math.pow((atom_feature.y - hetatm_feature.y), 2) + math.pow((atom_feature.z - hetatm_feature.z), 2))
                    if distance <= pharmacophore_data[2]:
                        pharmacophores_count[i] += 1
    return pharmacophores_count


def check_structure(structure_name2, q):
    path = '../data/PDB4/' + structure_name2
    file_extensions = []
    with open(path + '/file_extensions_pharmacophores', 'r') as file_extensions_file:
        for line in file_extensions_file.readlines():
            file_extensions.append(line.strip('\n'))
    binding_affinity2 = pdb_bind_dict[structure_name2]
    if int(binding_affinity2[1]) == 0:
        log_binding_affinity = float(1)
    else:
        log_binding_affinity = math.log(binding_affinity2[1])
    pharmacophores_base = [0, 0, 0, 0, 0, 0, 0, 0]
    for file_extension in file_extensions:
        pharmacophores_count = check_features(path, file_extension)
        if pharmacophores_count != pharmacophores_base:
            string_output = structure_name2 + ','
            for pharmacophore_count in pharmacophores_count:
                string_output += str(pharmacophore_count) + ','
            string_output += str(binding_affinity2[1]) + ',' + str(log_binding_affinity)
            q.put((binding_affinity2[0], string_output))
            # print(binding_affinity2[0], string_output)
            break


def listener(q):
    string = "# receptor_pharmacophore-ligand_pharmacophore counts"
    string2 = "PDB,Hydrophobe-Hydrophobe,Acceptor-Donor,Donor-Acceptor,PosIonizable-NegIonizable,NegIonizable-PosIonizable,Aromatic-Aromatic,Aromatic-PosIonizable,PosIonizable-Aromatic,"
    f_Ki = open('../datasets/protLigBindDB_Ki', 'w')
    f_Ki.write(string + '\n')
    f_Ki.write(string2 + 'Ki,log(Ki)' + '\n')
    f_Kd = open('../datasets/protLigBindDB_Kd', 'w')
    f_Kd.write(string + '\n')
    f_Kd.write(string2 + 'Kd,log(Kd)' + '\n')
    f_IC50 = open('../datasets/protLigBindDB_IC50', 'w')
    f_IC50.write(string + '\n')
    f_IC50.write(string2 + 'IC50,log(IC50)' + '\n')
    while 1:
        m = q.get()
        if m == 'kill':
            break
        if m[0] == 'Ki':
            f_Ki.write(str(m[1]) + '\n')
            f_Ki.flush()
        elif m[0] == 'Kd':
            f_Kd.write(str(m[1]) + '\n')
            f_Kd.flush()
        elif m[0] == 'IC50':
            f_IC50.write(str(m[1]) + '\n')
            f_IC50.flush()
        else:
            print(m)
    f_Ki.close()
    f_Kd.close()
    f_IC50.close()


if __name__ == '__main__':
    structure_names = []
    with open('../data/PDB4/crossover_structures', 'r') as structure_names_file:
        for line in structure_names_file.readlines():
            structure_names.append(line.strip('\n'))

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() - 1)
    watcher = pool.apply_async(listener, (q,))

    i = 0
    jobs = []
    for structure_name in structure_names:
        job = pool.apply_async(check_structure, (structure_name, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('kill')
    pool.close()
