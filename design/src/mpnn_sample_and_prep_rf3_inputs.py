

import os
import json
from omegaconf import OmegaConf
import argparse
from typing import Tuple, List, Dict, Optional

from utils import extract_sequence_from_structure



def sample_sequences_with_proteinmpnn(input_structure: str, chain_to_design: str, chain_start_and_end: Tuple[int, int], num_samples: int, temperature: float, temp_dir: str) -> List[str]:

    os.makedirs(temp_dir, exist_ok=True)

    # if fasta file already exists, delete it
    fasta_file = os.path.join(temp_dir, f'{os.path.basename(input_structure)[:-4]}.fa')
    if os.path.exists(fasta_file):
        os.remove(fasta_file)
    
    command = f"python /gscratch/stf/gvisan01/foundry/models/mpnn/src/mpnn/inference.py"
    command += f" --model_type protein_mpnn"
    command += f" --checkpoint_path /gscratch/stf/gvisan01/foundry/models/mpnn/proteinmpnn_v_48_020.pt"
    command += f" --is_legacy_weights True"
    command += f" --out_directory {temp_dir}"
    command += f" --write_fasta True"
    command += f" --write_structures False"
    command += f" --structure_path {input_structure}"
    command += f" --temperature {temperature}"
    command += f" --batch_size {num_samples}"
    command += f" --number_of_batches 1"
    command += f" --designed_chains {chain_to_design}"
    os.system(command)

    ## read fasta, extract sampled sequences of designed chain
    sampled_sequences = []
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('>'):
                continue
            else:
                sampled_sequences.append(line.strip()[chain_start_and_end[0] : chain_start_and_end[1]])
    
    sampled_sequences = list(set(sampled_sequences)) # remove duplicates
    
    return sampled_sequences

def make_rf3_input_json(name_to_chain_to_seq: Dict[str, Dict[str, str]], outpath: str, templating: Optional[Dict[str, List[str]]] = None):

    '''
    Assuming templates are the same across all "names" (each name is a distinct design)
    '''

    data = [
        {
            "name": name,
            "components": [
                {
                    "chain_id": chain_id,
                    "seq": seq,
                }
                for chain_id, seq in chain_to_seq.items()
            ]
        }
        for name, chain_to_seq in name_to_chain_to_seq.items()
    ]

    if templating is not None:
        
        chains_with_template = set([sel[0] for sel in templating["template_selection"]])

        new_data = []
        for data_dict in data:
            name = data_dict['name']

            # keep only components without template
            new_components = []
            for component in data_dict['components']:
                if component['chain_id'] not in chains_with_template:
                    new_components.append(component)
            # add paths as extra components
            for path in templating["paths"]:
                new_components.append({"path": path})

            new_data.append({
                "name": name,
                "components": new_components,
                "template_selection": list(templating["template_selection"]),
            })

        data = new_data


    with open(outpath, 'w+') as f:
        json.dump(data, f, indent=4)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    args = parser.parse_args()

    cfg = OmegaConf.load(args.config)


    os.makedirs(cfg.folding_dir, exist_ok=True)

    name_to_chain_to_seq = {}
    name_to_designed_structure = {}
    for filename in os.listdir(cfg.diffusion_dir):
        if filename.endswith('.cif'):
            design_name = filename[:-4] # remove .cif
            designed_structure = os.path.join(cfg.diffusion_dir, filename)
            folding_name = f'{design_name}'

            chain_to_seq = extract_sequence_from_structure(designed_structure)

            # compute ends of designed chain in the global sequence, to slice proteinmpnn's output
            chains = sorted(list(chain_to_seq.keys()))
            start_i = 0
            for chain in chains:
                if chain != cfg.designed_chain:
                    start_i += len(chain_to_seq[chain])
                else:
                    end_i = start_i + len(chain_to_seq[chain])
                    break
            
            name_to_chain_to_seq[f'{folding_name}__seq_rfd3'] = chain_to_seq
            name_to_designed_structure[f'{folding_name}__seq_rfd3'] = designed_structure

            if cfg.num_proteinmpnn_samples > 0:
                # sample protein sequences with proteinmpnn, add to name_to_chain_to_seq with new name for each sample, and also save designed structure path for each sample to compute metrics later
                sampled_sequences = sample_sequences_with_proteinmpnn(designed_structure, cfg.designed_chain, [start_i, end_i], cfg.num_proteinmpnn_samples, cfg.proteinmpnn_temperature, cfg.proteinmpnn_temp_dir)
                for i, sampled_seq in enumerate(sampled_sequences):
                    name_to_chain_to_seq[f'{folding_name}__seq_mpnn_{i}'] = {**chain_to_seq, cfg.designed_chain: sampled_seq}
                    name_to_designed_structure[f'{folding_name}__seq_mpnn_{i}'] = designed_structure

    make_rf3_input_json(name_to_chain_to_seq, cfg.input_json_for_folding, cfg.templating if cfg.use_templating else None)

    with open(cfg.name_to_designed_structure_json, 'w+') as f:
        json.dump(name_to_designed_structure, f, indent=4)
    
    print('-'* 40)
    print('\tsaved input json to: ' + cfg.input_json_for_folding)
    print('\tsaved name to designed structure json to: ' + cfg.name_to_designed_structure_json)
    print('-'* 40)

    ## change the target
    if "alternative_target" in cfg or "templating_alt_target" in cfg:
        name_to_chain_to_seq_alt_target = {}
        for name, chain_to_seq in name_to_chain_to_seq.items():
            alt_chain_to_seq = {}
            for chain_id, seq in chain_to_seq.items():
                if "alternative_target" in cfg and chain_id in cfg.alternative_target:
                    alt_chain_to_seq[chain_id] = cfg.alternative_target[chain_id]
                else:
                    alt_chain_to_seq[chain_id] = seq
            name_to_chain_to_seq_alt_target[name] = alt_chain_to_seq
        
        os.makedirs(cfg.folding_alt_target_dir, exist_ok=True)
        make_rf3_input_json(name_to_chain_to_seq_alt_target, cfg.input_json_for_folding_alt_target, cfg.templating_alt_target if cfg.use_templating else None)

        print('-'* 40)
        print('\tsaved input json with alternative target to: ' + cfg.input_json_for_folding_alt_target)
        print('-'* 40)


