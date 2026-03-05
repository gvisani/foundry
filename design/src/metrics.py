
import biotite.structure as struc

import numpy as np
from atomworks.io.utils.io_utils import load_any

from biotite.structure import AtomArray

from utils import make_mask

from typing import List, Dict, Tuple


def rmsd_of_chains_ca_only(reference: AtomArray, subject: AtomArray, chains_to_rmsd: List[str]) -> float:

    reference_chain = reference[np.logical_or.reduce([reference.chain_id == chain_id for chain_id in chains_to_rmsd], axis=0)]
    subject_chain = subject[np.logical_or.reduce([subject.chain_id == chain_id for chain_id in chains_to_rmsd], axis=0)]

    reference_chain_ca = reference_chain[reference_chain.atom_name == 'CA']
    subject_chain_ca = subject_chain[subject_chain.atom_name == 'CA']

    return float(struc.rmsd(reference_chain_ca, subject_chain_ca))


def refolding_rmsd(pdbpath_1: str, pdbpath_2: str, chains_to_align: List[str], chains_to_rmsd: List[str]) -> float:

    struc_1 = load_any(pdbpath_1, extra_fields="all")[0] # AtomArray, first model
    struc_2 = load_any(pdbpath_2, extra_fields="all")[0] # AtomArray, first model

    struc_1_to_align = struc_1[np.logical_or.reduce([struc_1.chain_id == chain_id for chain_id in chains_to_align], axis=0)]
    struc_2_to_align = struc_2[np.logical_or.reduce([struc_2.chain_id == chain_id for chain_id in chains_to_align], axis=0)]

    struc_1_to_align_ca = struc_1_to_align[struc_1_to_align.atom_name == 'CA']
    struc_2_to_align_ca = struc_2_to_align[struc_2_to_align.atom_name == 'CA']

    # compute transformation
    _, transformation = struc.superimpose(struc_2_to_align_ca, struc_1_to_align_ca)

    # apply transformation
    struc_1_aligned = transformation.apply(struc_1)

    return rmsd_of_chains_ca_only(struc_1_aligned, struc_2, chains_to_rmsd)


def distance_statistics(pdbpath: str, elements_1: str | List[str | Tuple[str, int]], elements_2: str | List[str | Tuple[str, int]]) -> Dict[str, float]:

    struc = load_any(pdbpath, extra_fields="all")[0] # AtomArray, first model
    struc = struc[struc.atom_name == 'CA'] # ca only

    mask_1 = make_mask(struc, elements_1)    
    mask_2 = make_mask(struc, elements_2)    

    return distance_statistics_from_coords(struc[mask_1].coord, struc[mask_2].coord)


def distance_statistics_from_coords(coord1_N3: np.ndarray, coord2_M3: np.ndarray, k: int = 10) -> Dict[str, float]:

    dist_matrix_NM = np.linalg.norm(
        coord1_N3[:, np.newaxis, :] - coord2_M3[np.newaxis, :, :],
        axis=-1,
    )

    sorted_distances_NxM = np.sort(dist_matrix_NM.ravel())

    return {
        'min': float(sorted_distances_NxM[0]),
        f'min_{k}': float(np.mean(sorted_distances_NxM[:k])),
        'avg': float(np.mean(sorted_distances_NxM))
    }


if __name__ == '__main__':

    ## testing
    refolding_rmsd(
        '/gscratch/stf/gvisan01/foundry/on_off_target/inference_outputs/magea3_fix_bb_fix_side/0/magea3_fix_bb_fix_side_magea3_0_model_0.cif',
        '/gscratch/stf/gvisan01/foundry/on_off_target/folding_output/0__magea3_fix_bb_fix_side_magea3_0_model_0__templating_is_True/0__magea3_fix_bb_fix_side_magea3_0_model_0__templating_is_True_model.cif',
        ['A'],
        ['D'],
    )




