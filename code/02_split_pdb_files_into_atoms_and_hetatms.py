import multiprocessing as mp
import os
import shutil

from Bio.PDB import *


def extract_hetatms(structure_name):
    parser = PDBParser(QUIET=True)
    # structure_hetatms = parser.get_structure(structure_name, './pdb/pdb' + structure_name + '.ent')
    structure_hetatms = parser.get_structure(structure_name, '../data/PDB2/' + structure_name + '.pdb')
    for model in structure_hetatms:
        for chain in model:
            residues_to_detach = []
            for residue in list(chain):
                residue_id = residue.get_id()
                hetfield = residue_id[0]
                if not hetfield.startswith('H_'):
                    residues_to_detach.append(residue_id)
            for id in residues_to_detach:
                chain.detach_child(id)
    save_structure(structure_hetatms, structure_name, 'hetatms')

    residues_id = []
    for residue in list(structure_hetatms.get_residues()):
        residues_id.append((residue.parent.parent.id, residue.parent.id, residue.get_id()[1]))
    residues_id = list(set(residues_id))

    for residue_id in residues_id:
        temp_structure = parser.get_structure(residue_id,
                                              '../data/PDB3/' + structure_name + '/' + 'hetatms.pdb')
        models_to_detach = []
        for model in temp_structure:
            if residue_id[0] == model.get_id():
                chains_to_detach = []
                for chain in model:
                    if residue_id[1] == chain.get_id():
                        residues_to_detach = []
                        for residue in list(chain):
                            if residue_id[2] != residue.get_id()[1]:
                                residues_to_detach.append(residue.get_id())
                        for id in residues_to_detach:
                            chain.detach_child(id)
                    else:
                        chains_to_detach.append(chain.get_id())
                for id in chains_to_detach:
                    model.detach_child(id)
            else:
                models_to_detach.append(model.get_id())
        for id in models_to_detach:
            temp_structure.detach_child(id)
        save_structure(temp_structure, structure_name,
                       'hetatm_' + str(residue_id[0]) + str(residue_id[1]) + str(residue_id[2]))
    os.remove('../data/PDB3/' + structure_name + '/' + 'hetatms.pdb')
    return residues_id


def save_atoms(structure_name, residues_id):
    parser = PDBParser(QUIET=True)
    # structure_atoms = parser.get_structure(structure_name, './pdb/pdb' + structure_name + '.ent')
    structure_atoms = parser.get_structure(structure_name, '../data/PDB2/' + structure_name + '.pdb')
    for model in structure_atoms:
        for chain in model:
            residues_to_detach = []
            for residue in list(chain):
                residue_id = residue.get_id()
                hetfield = residue_id[0]
                if not hetfield.startswith(' '):
                    residues_to_detach.append(residue_id)
            for id in residues_to_detach:
                chain.detach_child(id)
    save_structure(structure_atoms, structure_name, 'atoms')

    for residue_id in residues_id:
        hetatm_structure = parser.get_structure(residue_id,
                                                '../data/PDB3/' + structure_name + '/' + 'hetatm_' + str(
                                                    residue_id[0]) + str(residue_id[1]) + str(residue_id[2]) + '.pdb')
        if len(list(hetatm_structure.get_residues())) == 1:
            ligand = list(hetatm_structure.get_residues())[0]
            structure_atoms = parser.get_structure(structure_name,
                                                   '../data/PDB3/' + structure_name + '/' + 'atoms' + '.pdb')
            models_to_detach = []
            for model in structure_atoms:
                if residue_id[0] == model.get_id():
                    chains_to_detach = []
                    for chain in model:
                        if chain.get_id() == residue_id[1]:
                            residues_to_detach = []
                            for residue in list(chain):
                                residue_close = False
                                for atom in list(residue):
                                    for hetatm in ligand:
                                        distance = atom - hetatm
                                        if distance < 12:
                                            residue_close = True
                                            break
                                    if residue_close:
                                        break
                                if not residue_close:
                                    residues_to_detach.append(residue.get_id())
                            for id in residues_to_detach:
                                chain.detach_child(id)
                            if len(chain) == 0:
                                chains_to_detach.append(chain.get_id())
                        else:
                            chains_to_detach.append(chain.get_id())
                    for id in chains_to_detach:
                        model.detach_child(id)
                else:
                    models_to_detach.append(model.get_id())
            for id in models_to_detach:
                structure_atoms.detach_child(id)
            save_structure(structure_atoms, structure_name,
                           'atom_' + str(residue_id[0]) + str(residue_id[1]) + str(residue_id[2]))
    os.remove('../data/PDB3/' + structure_name + '/' + 'atoms' + '.pdb')


def save_structure(structure, structure_name, extension):
    output = PDBIO()
    output.set_structure(structure)
    path = '../data/PDB3/' + structure_name
    os.makedirs(path, exist_ok=True)
    target_lig_file = open(path + '/' + extension + '.pdb', 'w')
    target_lig_file.write('REMARK ' + structure_name + ' ' + extension + '\n')
    output.save(target_lig_file)
    target_lig_file.close()


def extract_hetatms_and_nearby_atoms(structure_name, q):
    # print(structure_name)
    residues_id = extract_hetatms(structure_name)
    save_atoms(structure_name, residues_id)
    file_extensions = []
    for residue_id in residues_id:
        string = str(residue_id[0]) + str(residue_id[1]) + str(residue_id[2])
        atom_path = '../data/PDB3/' + structure_name + '/atom_' + string + '.pdb'
        hetatm_path = '../data/PDB3/' + structure_name + '/hetatm_' + string + '.pdb'
        if os.path.isfile(atom_path) and os.path.isfile(hetatm_path):
            file_extensions.append(string)
        else:
            if os.path.isfile(atom_path):
                os.remove(atom_path)
            if os.path.isfile(hetatm_path):
                os.remove(hetatm_path)
    if len(file_extensions) > 0:
        with open('../data/PDB3/' + structure_name + '/file_extensions_pdb',
                  'w') as file_extensions_file:
            for file_extension in file_extensions:
                file_extensions_file.write(file_extension + '\n')
        q.put(structure_name)
    else:
        shutil.rmtree('../data/PDB3/' + structure_name)


def listener(q):
    f = open('../data/PDB3/structure_names_pdb', 'w')
    while 1:
        m = q.get()
        if m == 'kill':
            break
        f.write(str(m) + '\n')
        f.flush()
    f.close()


if __name__ == '__main__':
    structure_names = []
    with open('../data/PDB2/Contents', 'r') as structure_names_file:
        for line in structure_names_file.readlines():
            structure_names.append(line.strip('\n'))

    # with open('../data/PDB3/structure_names', 'w') as f:
    #     for structure_name in structure_names:
    #         if os.path.isdir('../data/PDB3/'+structure_name):
    #             f.write(structure_name+'\n')
    #             f.flush()

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() - 1)
    watcher = pool.apply_async(listener, (q,))

    i = 0
    jobs = []

    for structure_name in structure_names:
        job = pool.apply_async(extract_hetatms_and_nearby_atoms, (structure_name, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('kill')
    pool.close()
