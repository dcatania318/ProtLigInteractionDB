import multiprocessing as mp
import os

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, MolSurf, Lipinski


def check_lipisnki(structure_name, q):
    # print(structure_name)
    path = '../data/PDB3/' + structure_name
    file_extensions = []
    with open(path + '/file_extensions_sdf', 'r') as file_extensions_file:
        for line in file_extensions_file.readlines():
            file_extensions.append(line.strip('\n'))
    file_extensions2 = []
    for file_extension in file_extensions:
        if check_ligand(path + '/hetatm_' + file_extension + '.sdf'):
            file_extensions2.append(file_extension)
    if len(file_extensions2) > 0:
        with open('../data/PDB3/' + structure_name + '/file_extensions_lipinski_44_45',
                  'w') as file_extensions_file:
            for file_extension in file_extensions2:
                file_extensions_file.write(file_extension + '\n')
        q.put(structure_name)


def check_ligand(file_path):
    bool = False
    if os.path.isfile(file_path):
        suppl = Chem.SDMolSupplier(file_path)
        for mol in suppl:
            if mol is not None:
                # components of rule
                hydrogen_bond_doner = True if Lipinski.NumHDonors(mol) <= 5 else False
                hydrogen_bond_acceptors = True if Lipinski.NumHAcceptors(mol) <= 10 else False
                molecular_mass = True if Descriptors.ExactMolWt(mol) <= 500 else False
                octanol_water_partition_coefficient_logP = True if Crippen.MolLogP(mol) <= 5 else False
                components_rank = hydrogen_bond_doner + hydrogen_bond_acceptors + molecular_mass + octanol_water_partition_coefficient_logP

                # variants
                partition_coefficient_logP = True if -0.4 <= Crippen.MolLogP(mol) <= 5.6 else False
                molar_refractivity = True if 40 <= Crippen.MolMR(mol) <= 130 else False
                molecular_weight = True if 180 <= Descriptors.ExactMolWt(mol) <= 500 else False
                number_of_atoms = True if 20 <= Lipinski.HeavyAtomCount(mol) <= 70 else False
                polar_surface_area = True if MolSurf.TPSA(mol) <= 140 else False
                variants_rank = partition_coefficient_logP + molar_refractivity + molecular_weight + number_of_atoms + polar_surface_area

                if (components_rank == 4) and (variants_rank == 4 or variants_rank == 5):
                    bool = True
    return bool


def listener(q):
    f = open('../data/PDB3/structure_names_lipinski_44_45', 'w')
    while 1:
        m = q.get()
        if m == 'kill':
            break
        f.write(str(m) + '\n')
        f.flush()
    f.close()


if __name__ == '__main__':
    structure_names = []
    with open('../data/PDB3/structure_names_sdf', 'r') as structure_names_file:
        for line in structure_names_file.readlines():
            structure_names.append(line.strip('\n'))

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() - 1)
    watcher = pool.apply_async(listener, (q,))

    i = 0
    jobs = []

    for structure_name in structure_names:
        job = pool.apply_async(check_lipisnki, (structure_name, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('kill')
    pool.close()
