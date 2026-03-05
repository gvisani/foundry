
import numpy as np
import torch
from biotite import structure
from biotite.structure import AtomArray
from atomworks.io import parse

from typing import Dict, List, Tuple


def rigid_align(struc_from: AtomArray, struc_to: AtomArray, mask_from: np.ndarray, mask_to: np.ndarray) -> AtomArray:
    """Aligns struc_from to struc_to alongside atoms defined by masks."""
    R, center_from, center_to = rigid_align_params(struc_from, struc_to, mask_from, mask_to)
    return apply_rigid_align(struc_from, R, center_from, center_to)


def rigid_align_params(struc_from: AtomArray, struc_to: AtomArray, mask_from: np.ndarray, mask_to: np.ndarray):
    assert struc_from.shape[0] == mask_from.shape[0]
    assert struc_to.shape[0] == mask_to.shape[0]
    assert np.sum(mask_from) == np.sum(mask_to), f"{np.sum(mask_from)} != {np.sum(mask_to)}"

    coords_from = struc_from[mask_from].coord  # [N, 3]
    coords_to = struc_to[mask_to].coord        # [N, 3]

    center_from = coords_from.mean(axis=0)  # [3]
    center_to = coords_to.mean(axis=0)      # [3]

    centered_from = coords_from - center_from  # [N, 3]
    centered_to = coords_to - center_to        # [N, 3]

    C = centered_from.T @ centered_to  # [3, 3]
    U, S, Vt = np.linalg.svd(C)

    # Ensure proper rotation (det = +1)
    d = np.sign(np.linalg.det(U @ Vt))
    F_det = np.diag([1.0, 1.0, d])
    R = U @ F_det @ Vt  # [3, 3]

    return R, center_from, center_to


def apply_rigid_align(struc: AtomArray, R: np.ndarray, center_from: np.ndarray, center_to: np.ndarray) -> AtomArray:
    struc = struc.copy()
    struc.coord = (struc.coord - center_from) @ R + center_to
    return struc


def extract_sequence_from_structure(struc_file: str) -> Dict[str, str]:

    result_dict = parse(
        filename=struc_file,
        build_assembly=["1"],
        add_missing_atoms=True,
        remove_waters=True,
        hydrogen_policy="remove",
        model=1,
    )

    chain_to_seq = {}
    for chain in result_dict['chain_info']:
        chain_pure_str = str(chain) # from np.str_ to str
        chain_to_seq[chain_pure_str] = result_dict['chain_info'][chain]['processed_entity_non_canonical_sequence']
    
    return chain_to_seq


def make_mask(atom_array: AtomArray, elements: str | List[str | Tuple[str, int]]) -> np.ndarray:
    '''
    elements can be:
    - str, in which case it is interpreted as a chain_id
    - List with elements:
        - str, interpreted as a chain
        - tuple(str, int), interpreted as (chain, res_id)
    '''
    mask = np.zeros(atom_array.shape[0], dtype=bool)

    if isinstance(elements, str):
        elements = [elements]

    for elem in elements:
        if isinstance(elem, str):
            mask_update = atom_array.chain_id == elem
        elif isinstance(elem, tuple):
            chain_id = elem[0]
            res_id = elem[1]
            mask_update = np.logical_and(atom_array.chain_id == chain_id, atom_array.res_id == res_id)
        else:
            raise ValueError(f"Invalid type in input 'elements', must be str or tuple, got {type(elem)}")
        mask = np.logical_or(mask, mask_update)
    
    return mask



if __name__ == '__main__':

    pdbpath = '/gscratch/stf/gvisan01/foundry/design/campaigns/ha_stem_mabs/input_pdbs/4fqi_ab_align_renum_filled.pdb'
    chain_to_seq = extract_sequence_from_structure(pdbpath)
    for chain, seq in chain_to_seq.items():
        print(f"{chain}: {seq}")

    # result_dict = parse(
    #     filename="/gscratch/stf/gvisan01/foundry/design/campaigns/ha_stem_mabs/input_pdbs/4fqi_ab_align_renum.pdb",
    #     build_assembly=["1"],
    #     add_missing_atoms=True,
    #     remove_waters=True,
    #     hydrogen_policy="remove",
    #     model=1,
    # )
    # struc = result_dict["assemblies"]["1"][0]
    # print(type(struc))
    # print(len(struc))

