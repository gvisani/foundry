#!/usr/bin/env python
"""
Fill in residue number gaps within each chain of a PDB file.
For example, if chain L has resnums [1,2,9,11,12], this produces [1,2,3,4,5].

Assumes no insertion codes (run renumber_pdb.py first if needed).

Usage:
    python fill_resnum_gaps.py input.pdb [-o OUTPUT_PDB] [-m MAPPING_JSON]
"""

import argparse
import json
from pathlib import Path
from collections import OrderedDict

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb


def fill_gaps(structure: struc.AtomArray):
    """
    Renumber residues so that there are no gaps within each chain.
    Residues are kept in their original order and assigned consecutive
    numbers starting from the first residue's resnum.

    Returns
    -------
    new_structure : struc.AtomArray
        Copy with updated res_id.
    mapping : dict
        {chain_id: [{"old_resnum": int, "new_resnum": int}, ...]}
    """
    new_structure = structure.copy()
    new_res_id = new_structure.res_id.copy()

    starts = struc.get_residue_starts(structure)

    chain_to_starts: dict[str, list[int]] = OrderedDict()
    for s in starts:
        ch = structure.chain_id[s]
        chain_to_starts.setdefault(ch, []).append(s)

    mapping: dict[str, list[dict]] = {}

    for chain, res_starts in chain_to_starts.items():
        # Collect (res_id, start_idx) — already in file order
        residues = [(int(structure.res_id[s]), s) for s in res_starts]
        # Sort by res_id to be safe
        residues.sort(key=lambda x: x[0])

        chain_mapping = []
        # Start numbering from the first residue's resnum
        current = residues[0][0]

        for rid, start in residues:
            res_mask = struc.get_residue_masks(structure, np.array([start]))[0]
            new_res_id[res_mask] = current

            chain_mapping.append({
                "old_resnum": rid,
                "new_resnum": current,
            })
            current += 1

        mapping[chain] = chain_mapping

    new_structure.res_id = new_res_id
    return new_structure, mapping


def main():
    parser = argparse.ArgumentParser(
        description="Fill in residue number gaps in a PDB file."
    )
    parser.add_argument("input_pdb", help="Path to input PDB file")
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output PDB path (default: <input>_filled.pdb)",
    )
    parser.add_argument(
        "-m", "--mapping",
        default=None,
        help="Output JSON mapping path (default: <input>_filled_mapping.json)",
    )
    args = parser.parse_args()

    input_path = Path(args.input_pdb)
    stem = input_path.stem

    output_pdb = args.output or str(input_path.parent / f"{stem}_filled.pdb")
    output_json = args.mapping or str(input_path.parent / f"{stem}_filled_mapping.json")

    pdb_file = pdb.PDBFile.read(str(input_path))
    structure = pdb_file.get_structure(model=1)

    new_structure, mapping = fill_gaps(structure)

    n_changed = sum(
        1
        for chain_entries in mapping.values()
        for e in chain_entries
        if e["old_resnum"] != e["new_resnum"]
    )
    n_total = sum(len(v) for v in mapping.values())

    out_pdb_file = pdb.PDBFile()
    out_pdb_file.set_structure(new_structure)
    out_pdb_file.write(output_pdb)

    with open(output_json, "w") as f:
        json.dump(mapping, f, indent=2)

    print(f"Renumbered {n_changed}/{n_total} residues across {len(mapping)} chain(s).")
    print(f"  PDB  → {output_pdb}")
    print(f"  JSON → {output_json}")


if __name__ == "__main__":
    main()
