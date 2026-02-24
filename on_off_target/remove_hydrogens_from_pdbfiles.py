
import os
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
from tqdm import tqdm
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb_dir', type=str, default='./input_pdbs')
    args = parser.parse_args()

    all_files = [filename for filename in os.listdir(args.pdb_dir) if filename.endswith('.pdb')]

    for filename in tqdm(all_files):

        # Read PDB
        structure = pdb.PDBFile.read(os.path.join(args.pdb_dir, filename)).get_structure()[0]

        is_hydrogen_mask = (structure.element == 'H')

        pdb_file = pdb.PDBFile()
        pdb.set_structure(pdb_file, structure[~is_hydrogen_mask])
        pdb_file.write(os.path.join(args.pdb_dir, filename))
