



import os
import json
from omegaconf import OmegaConf
import argparse
from typing import Tuple, List, Dict, Optional
import numpy as np
from biotite.structure import AtomArray, array
from biotite.structure.io.pdbx import CIFFile, set_structure
from atomworks.io import parse
from utils import extract_sequence_from_structure
from mpnn_sample_and_prep_rf3_inputs import make_rf3_input_json
from pathlib import Path



def sample_sequences_with_proteinmpnn_multi_target(
        pdbfile_1: str,
        pdbfile_2: str,
        chain_to_design_1: str,
        chain_to_design_2: str,
        num_samples: int,
        temperature: float,
        temp_dir: str
    ) -> List[str]:

    os.makedirs(temp_dir, exist_ok=True)

    temp_pdbfile = os.path.join(temp_dir, f"{Path(pdbfile_1).stem}__{Path(pdbfile_2).stem}.cif")

    struc_2_chain_map = merge_structures(pdbfile_1, pdbfile_2, temp_pdbfile)
    chain_to_design_2 = struc_2_chain_map[chain_to_design_2]

    # if fasta file already exists, delete it
    fasta_file = os.path.join(temp_dir, f'{Path(temp_pdbfile).stem}.fa')
    if os.path.exists(fasta_file):
        os.remove(fasta_file)
        
    command = f"python /gscratch/stf/gvisan01/foundry/models/mpnn/src/mpnn/inference.py"
    command += f" --model_type protein_mpnn"
    command += f" --checkpoint_path /gscratch/stf/gvisan01/foundry/models/mpnn/proteinmpnn_v_48_020.pt"
    command += f" --is_legacy_weights True"
    command += f" --out_directory {temp_dir}"
    command += f" --write_fasta True"
    command += f" --write_structures False"
    command += f" --structure_path {temp_pdbfile}"
    command += f" --temperature {temperature}"
    command += f" --batch_size {num_samples}"
    command += f" --number_of_batches 1"
    command += f" --designed_chains {','.join([chain_to_design_1, chain_to_design_2])}"
    command += f" --homo_oligomer_chains {"\'[[" + ",".join([f'"{chain_to_design_1}"', f'"{chain_to_design_2}"']) + "]]\'"}"
    os.system(command)

    chain_to_seq = extract_sequence_from_structure(temp_pdbfile)
    chain_to_start_and_end = get_start_and_end_of_chains_in_seq([chain_to_design_1, chain_to_design_2], chain_to_seq)

    ## read fasta, extract sampled sequences of designed chain
    chain_to_sampled_sequences = {}

    for chain in [chain_to_design_1, chain_to_design_2]:

        sampled_sequences = []
        with open(fasta_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('>'):
                    continue
                else:
                    sampled_sequences.append(line.strip()[chain_to_start_and_end[chain][0] : chain_to_start_and_end[chain][1]])
        
        chain_to_sampled_sequences[chain] = sampled_sequences
    
    assert chain_to_sampled_sequences[chain_to_design_1][0] == chain_to_sampled_sequences[chain_to_design_2][0]
    
    return list(set(chain_to_sampled_sequences[chain_to_design_1])) # remove duplicates


def get_start_and_end_of_chains_in_seq(chains: List[str], chain_to_seq: Dict[str, str]) -> Dict[str, Tuple[int, int]]:
    
    # compute ends of designed chain in the global sequence, to slice proteinmpnn's output
    chain_to_start_and_end = {}

    all_chains = sorted(list(chain_to_seq.keys()))

    for chain in chains:

        start_i = 0
        for c in all_chains:
            if c != chain:
                start_i += len(chain_to_seq[c])
            else:
                end_i = start_i + len(chain_to_seq[c])
                break
            
        chain_to_start_and_end[chain] = (start_i, end_i)
    
    return chain_to_start_and_end


def max_distance(structure1: AtomArray, structure2: AtomArray) -> float:
    """
    Returns the maximum pairwise distance between any atom in structure1
    and any atom in structure2.
    """
    coords1 = structure1.coord  # shape (N, 3)
    coords2 = structure2.coord  # shape (M, 3)

    # Compute all pairwise distances using broadcasting
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]  # (N, M, 3)
    distances = np.sqrt(np.sum(diff ** 2, axis=-1))               # (N, M)

    return float(np.max(distances))


def merge_structures(pdbfile_1: str, pdbfile_2: str, pdbfile_out: str) -> Dict[str, str]:

    struc_1 = parse(
        filename=pdbfile_1,
        build_assembly=["1"],
        add_missing_atoms=True,
        remove_waters=True,
        hydrogen_policy="remove",
        model=1,
    )["assemblies"]["1"][0]

    struc_2 = parse(
        filename=pdbfile_2,
        build_assembly=["1"],
        add_missing_atoms=True,
        remove_waters=True,
        hydrogen_policy="remove",
        model=1,
    )["assemblies"]["1"][0]

    ## renumber chains in struc_2 to follow struc_1 alphabetically
    chains_1 = sorted(list(map(lambda x: str(x), list(set(struc_1.chain_id)))))
    chains_2 = sorted(list(map(lambda x: str(x), list(set(struc_2.chain_id)))))
    all_letters = [chr(ord('A') + i) for i in range(26)]
    available = [c for c in all_letters if c not in chains_1]
    chain_map = {old: available[i] for i, old in enumerate(chains_2)}
    new_chain_ids = struc_2.chain_id.copy()
    for old, new in chain_map.items():
        new_chain_ids[struc_2.chain_id == old] = new
    struc_2.chain_id = new_chain_ids

    ## separate the two structure by at least 101 Angstroms

    # compute maximum atom-atom distance between the two structures
    max_dist = max_distance(
        struc_1[struc_1.atom_name == "CA"],
        struc_2[struc_2.atom_name == "CA"]
    )

    # add 101 + max_dist to all coordinates of struc_2
    struc_2.coord = struc_2.coord + 101 + max_dist

    ## merge and save
    merged = array(list(struc_1) + list(struc_2))
    cif_file = CIFFile()
    set_structure(cif_file, merged)
    cif_file.write(pdbfile_out)

    return chain_map





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    args = parser.parse_args()

    cfg = OmegaConf.load(args.config)


    os.makedirs(cfg.folding_dir, exist_ok=True)

    name_to_chain_to_seq = {}
    name_to_designed_structure = {}
    name_to_chain_to_seq_alt = {}
    name_to_designed_structure_alt = {}
    for filename in os.listdir(cfg.diffusion_dir):
        if filename.endswith('.cif') and not filename.endswith('_alt.cif'):

            design_name = filename[:-4] # remove .cif
            design_name_alt = design_name + '_alt'
            
            designed_structure = os.path.join(cfg.diffusion_dir, design_name + '.cif')
            designed_structure_alt = os.path.join(cfg.diffusion_dir, design_name_alt + '.cif')
            

            chain_to_seq = extract_sequence_from_structure(designed_structure)
            chain_to_seq_alt = extract_sequence_from_structure(designed_structure_alt)
            
            name_to_chain_to_seq[f'{design_name}__seq_rfd3'] = chain_to_seq
            name_to_designed_structure[f'{design_name}__seq_rfd3'] = designed_structure
            name_to_chain_to_seq_alt[f'{design_name_alt}__seq_rfd3'] = chain_to_seq_alt
            name_to_designed_structure_alt[f'{design_name_alt}__seq_rfd3'] = designed_structure_alt

            if cfg.num_proteinmpnn_samples > 0:
                # sample protein sequences with proteinmpnn, add to name_to_chain_to_seq with new name for each sample, and also save designed structure path for each sample to compute metrics later
                sampled_sequences = sample_sequences_with_proteinmpnn_multi_target(
                    designed_structure,
                    designed_structure_alt,
                    cfg.designed_chain,
                    cfg.designed_chain,
                    cfg.num_proteinmpnn_samples,
                    cfg.proteinmpnn_temperature,
                    cfg.proteinmpnn_temp_dir
                )
                for i, sampled_seq in enumerate(sampled_sequences):
                    # NOTE: sampled seq it's the sam for main and lat targets, that's the whole point of this
                    name_to_chain_to_seq[f'{design_name}__seq_mpnn_{i}'] = {**chain_to_seq, cfg.designed_chain: sampled_seq}
                    name_to_designed_structure[f'{design_name}__seq_mpnn_{i}'] = designed_structure
                    name_to_chain_to_seq_alt[f'{design_name_alt}__seq_mpnn_{i}'] = {**chain_to_seq_alt, cfg.designed_chain: sampled_seq}
                    name_to_designed_structure_alt[f'{design_name_alt}__seq_mpnn_{i}'] = designed_structure_alt

    make_rf3_input_json(name_to_chain_to_seq, cfg.input_json_for_folding, cfg.templating if cfg.use_templating else None)

    os.makedirs(cfg.folding_alt_target_dir, exist_ok=True)
    make_rf3_input_json(name_to_chain_to_seq_alt, cfg.input_json_for_folding_alt_target, cfg.templating_alt_target if cfg.use_templating else None)

    with open(cfg.name_to_designed_structure_json, 'w+') as f:
        json.dump(name_to_designed_structure, f, indent=4)
    
    with open(cfg.name_to_designed_structure_json_alt_target, 'w+') as f:
        json.dump(name_to_designed_structure_alt, f, indent=4)


