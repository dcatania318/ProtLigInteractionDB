import multiprocessing as mp
import os

import openbabel


def pdb_to_mol(path, file_name):
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("pdb", "mol")
    mol = openbabel.OBMol()
    ob_conversion.ReadFile(mol, path + file_name + '.pdb')
    mol.AddHydrogens()
    ob_conversion.WriteFile(mol, path + file_name + '.mol')


def standardize_structure(path):
    # Proprietary code belonging to Dr Jean-Paul Ebejer has been removed.               
    if os.stat(path + '.sdf').st_size == 0:
        os.remove(path + '.sdf')
        return False
    return True


def convert_and_standardize_structure(structure_name, q):
    # print(structure_name)
    path = '../data/PDB3/' + structure_name
    file_extensions = []
    with open(path + '/file_extensions_pdb', 'r') as file_extensions_file:
        for line in file_extensions_file.readlines():
            file_extensions.append(line.strip('\n'))
    file_extensions2 = []
    for file_extension in file_extensions:
        pdb_to_mol(path, '/atom_' + file_extension)
        pdb_to_mol(path, '/hetatm_' + file_extension)
        atom_standardized = standardize_structure(path + '/atom_' + file_extension)
        hetatm_standardized = standardize_structure(path + '/hetatm_' + file_extension)
        if atom_standardized and hetatm_standardized:
            file_extensions2.append(file_extension)
        elif atom_standardized and not hetatm_standardized:
            os.remove(path + '/atom_' + file_extension + '.sdf')
        elif not atom_standardized and hetatm_standardized:
            os.remove(path + '/hetatm_' + file_extension + '.sdf')
    if len(file_extensions2) > 0:
        with open('../data/PDB3/' + structure_name + '/file_extensions_sdf',
                  'w') as file_extensions_file:
            for file_extension in file_extensions2:
                file_extensions_file.write(file_extension + '\n')
        q.put(structure_name)


def listener(q):
    f = open('../data/PDB3/structure_names_sdf', 'w')
    while 1:
        m = q.get()
        if m == 'kill':
            break
        f.write(str(m) + '\n')
        f.flush()
    f.close()


if __name__ == '__main__':
    structure_names = []
    with open('../data/PDB3/structure_names_pdb', 'r') as structure_names_file:
        for line in structure_names_file.readlines():
            structure_names.append(line.strip('\n'))

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() - 1)
    watcher = pool.apply_async(listener, (q,))

    i = 0
    jobs = []
    for structure_name in structure_names:
        job = pool.apply_async(convert_and_standardize_structure, (structure_name, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('kill')
    pool.close()
