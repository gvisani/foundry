#!/usr/bin/env python
"""
Rename chains in a PDB file according to user-provided mappings.

Usage:
    python rename_chains.py input.pdb A:H L:K -o output.pdb
    python rename_chains.py input.pdb A:H L:K B:B   # B:B is a no-op but explicit

Mappings are given as OLD:NEW pairs. Chains not mentioned are kept as-is.
If remapping would create a collision with an existing (unmapped) chain,
the script will error out.
"""

import argparse
import json
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb


def parse_mappings(mapping_args: list[str]) -> dict[str, str]:
    """Parse 'OLD:NEW' mapping strings into a dict."""
    mappings = {}
    for m in mapping_args:
        parts = m.split(":")
        if len(parts) != 2 or not parts[0] or not parts[1]:
            raise ValueError(
                f"Invalid mapping '{m}'. Expected format OLD:NEW (e.g. A:H)"
            )
        old, new = parts
        if old in mappings:
            raise ValueError(f"Duplicate mapping for chain '{old}'")
        mappings[old] = new
    return mappings


def rename_chains(structure: struc.AtomArray, mappings: dict[str, str]):
    """
    Apply chain renaming to a structure.

    Parameters
    ----------
    structure : struc.AtomArray
    mappings : dict
        {old_chain_id: new_chain_id}

    Returns
    -------
    new_structure : struc.AtomArray
    applied : dict
        The full effective mapping for all chains (including identity).
    """
    new_structure = structure.copy()
    existing_chains = set(structure.chain_id)

    # Build the full mapping (identity for unmapped chains)
    full_mapping = {}
    for ch in existing_chains:
        full_mapping[ch] = mappings.get(ch, ch)

    # Check for chains in the mapping that don't exist
    for old in mappings:
        if old not in existing_chains:
            raise ValueError(
                f"Chain '{old}' not found in PDB. "
                f"Available chains: {sorted(existing_chains)}"
            )

    # Check for collisions in the output
    new_ids = list(full_mapping.values())
    if len(new_ids) != len(set(new_ids)):
        seen = {}
        for old, new in full_mapping.items():
            if new in seen:
                raise ValueError(
                    f"Collision: chains '{seen[new]}' and '{old}' "
                    f"would both map to '{new}'"
                )
            seen[new] = old

    # Apply
    new_chain_id = new_structure.chain_id.copy()
    for old, new in mappings.items():
        mask = structure.chain_id == old
        new_chain_id[mask] = new
    new_structure.chain_id = new_chain_id

    return new_structure, full_mapping


def main():
    parser = argparse.ArgumentParser(
        description="Rename chains in a PDB file.",
        epilog="Example: python rename_chains.py input.pdb A:H L:K",
    )
    parser.add_argument("input_pdb", help="Path to input PDB file")
    parser.add_argument(
        "mappings",
        nargs="+",
        help="Chain mappings as OLD:NEW pairs (e.g. A:H L:K)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output PDB path (default: <input>_renamed.pdb)",
    )
    parser.add_argument(
        "-m", "--mapping",
        default=None,
        help="Output JSON mapping path (default: <input>_chain_mapping.json)",
    )
    args = parser.parse_args()

    input_path = Path(args.input_pdb)
    stem = input_path.stem

    output_pdb = args.output or str(input_path.parent / f"{stem}_renamed.pdb")
    output_json = args.mapping or str(
        input_path.parent / f"{stem}_chain_mapping.json"
    )

    mappings = parse_mappings(args.mappings)

    pdb_file = pdb.PDBFile.read(str(input_path))
    structure = pdb_file.get_structure(model=1)

    new_structure, full_mapping = rename_chains(structure, mappings)

    out_pdb_file = pdb.PDBFile()
    out_pdb_file.set_structure(new_structure)
    out_pdb_file.write(output_pdb)

    with open(output_json, "w") as f:
        json.dump(full_mapping, f, indent=2)

    changed = [f"{old}→{new}" for old, new in mappings.items() if old != new]
    print(f"Renamed {len(changed)} chain(s): {', '.join(changed)}")
    print(f"  PDB  → {output_pdb}")
    print(f"  JSON → {output_json}")


if __name__ == "__main__":
    main()