import multiprocessing as mp
import os

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures

pharmacophores = ['Hydrophobe', 'Acceptor', 'Donor', 'PosIonizable', 'NegIonizable',
                  'Aromatic']  # PosIonizable-Cation and NegIonizable-Anion


def extract_features(mol):
    factory = ChemicalFeatures.BuildFeatureFactory('./LigityFeatures.fdef')
    feats = factory.GetFeaturesForMol(mol)
    features = []
    for feat in feats:
        feature = feat.GetFamily()
        if feature in pharmacophores:
            id = feat.GetId()
            x, y, z = list(feat.GetPos())
            string = str(id) + ',' + feature + ',' + str(x) + ',' + str(y) + ',' + str(z)
            features.append(string)
    return features


def save_features(features, structure_name3, file_name):
    path = '../data/PDB4/' + structure_name3 + '/'
    os.makedirs(path, exist_ok=True)
    with open(path + file_name + '.feats', 'w') as features_file:
        for feature in features:
            features_file.write(feature + '\n')


def check_strucure_features(structure_name2, q):
    # print(structure_name)
    path = '../data/PDB3/' + structure_name2
    file_extensions = []
    with open(path + '/file_extensions_lipinski_44_45', 'r') as file_extensions_file:
        for line in file_extensions_file.readlines():
            file_extensions.append(line.strip('\n'))
    file_extensions2 = []
    for file_extension in file_extensions:
        atom_path = path + '/atom_' + file_extension
        hetatm_path = path + '/hetatm_' + file_extension
        if os.path.isfile(atom_path + '.sdf') and os.path.isfile(hetatm_path + '.sdf'):
            suppl_atom = Chem.SDMolSupplier(atom_path + '.sdf')
            suppl_hetatm = Chem.SDMolSupplier(hetatm_path + '.sdf')
            if len(suppl_atom) == 1 and len(suppl_hetatm) == 1:
                atom = suppl_atom[0]
                hetatm = suppl_hetatm[0]
                if atom != None and hetatm != None:
                    atom_features = extract_features(atom)
                    hetatm_features = extract_features(hetatm)
                    if len(atom_features) > 0 and len(hetatm_features) > 0:
                        save_features(atom_features, structure_name2, '/atom_' + file_extension)
                        save_features(hetatm_features, structure_name2, '/hetatm_' + file_extension)
                        file_extensions2.append(file_extension)
    if len(file_extensions2) > 0:
        with open('../data/PDB4/' + structure_name2 + '/file_extensions_pharmacophores', 'w') as file_extensions_file:
            for file_extension in file_extensions2:
                file_extensions_file.write(file_extension + '\n')
        q.put(structure_name2)


def listener(q):
    f = open('../data/PDB4/structure_names_pharmacophores', 'w')
    while 1:
        m = q.get()
        if m == 'kill':
            break
        f.write(str(m) + '\n')
        f.flush()
    f.close()


if __name__ == '__main__':
    os.makedirs('../data/PDB4/', exist_ok=True)
    structure_names = []
    with open('../data/PDB3/structure_names_lipinski_44_45', 'r') as structure_names_file:
        for line in structure_names_file.readlines():
            structure_names.append(line.strip('\n'))
    # check_strucure_features('2hu4')
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() - 1)
    watcher = pool.apply_async(listener, (q,))

    i = 0
    jobs = []
    for structure_name in structure_names:
        job = pool.apply_async(check_strucure_features, (structure_name, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('kill')
    pool.close()
