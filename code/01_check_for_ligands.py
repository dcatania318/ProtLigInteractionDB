import multiprocessing as mp
from Bio.PDB import *
import gzip
import os


def check_if_contains_ligands(code, file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(code, file)
    for model in structure:
        for chain in model:
            for residue in list(chain):
                residue_id = residue.get_id()
                hetfield = residue_id[0]
                if hetfield.startswith("H_"):
                    return True
    return False


def check_file(pdb_code, q):
    code = pdb_code.decode("utf-8").lower()
    archive_fn = "pdb%s.ent.gz" % code
    file_path = "../data/PDB/%s/" % code[1:3]
    file_name = file_path + archive_fn
    final_file = "../data/PDB2/%s.pdb" % code

    if os.path.isfile(file_name):
        gz = gzip.open(file_name, 'rb')
        with open(final_file, 'wb') as out:
            out.writelines(gz)
        gz.close()

        if check_if_contains_ligands(code, final_file):
            q.put(code)
        else:
            os.remove(final_file)
    # else:
    #     print("Error %s does not exist" % file_name)


def listener(q):
    f = open("../data/PDB2/Contents", 'w')
    while 1:
        m = q.get()
        if m == 'kill':
            break
        f.write(str(m) + '\n')
        f.flush()
    f.close()


if __name__ == '__main__':
    os.makedirs('../data/PDB2/', exist_ok=True)
    pdbl = PDBList()
    all_pdb_files = pdbl.get_all_entries()

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() - 1)
    watcher = pool.apply_async(listener, (q,))

    i = 0
    jobs = []
    for file in all_pdb_files:
        job = pool.apply_async(check_file, (file, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('kill')
    pool.close()
