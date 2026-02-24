

#!/usr/bin/env python3
"""
Offset residue numbers in a PDB file so that each chain's minimum resnum is 1,
but only if the chain contains residue numbers less than 1.
"""

import sys
import biotite.structure.io.pdb as pdb
import biotite.structure as struc


def fix_resnums(input_path: str, output_path: str):
    pdb_file = pdb.PDBFile.read(input_path)
    structure = pdb.get_structure(pdb_file, model=1)

    chain_ids = struc.get_chains(structure)

    for chain_id in chain_ids:
        chain_mask = structure.chain_id == chain_id
        chain_resnums = structure.res_id[chain_mask]
        min_resnum = chain_resnums.min()

        if min_resnum < 1:
            offset = 1 - min_resnum
            structure.res_id[chain_mask] += offset
            print(f"Chain {chain_id}: offset by +{offset} (min was {min_resnum})")
        else:
            print(f"Chain {chain_id}: no change needed (min is {min_resnum})")

    out_file = pdb.PDBFile()
    pdb.set_structure(out_file, structure)
    out_file.write(output_path)
    print(f"Written to {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.pdb> <output.pdb>")
        sys.exit(1)

    fix_resnums(sys.argv[1], sys.argv[2])