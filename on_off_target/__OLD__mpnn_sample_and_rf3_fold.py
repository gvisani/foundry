

import os
from atomworks.io import parse
import json
from tqdm import tqdm

from metrics import refolding_rmsd

from typing import Tuple, List, Dict, Optional


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

def make_rf3_input_json(name_to_chain_to_seq: Dict[str, Dict[str, str]], outpath: str, chain_to_template_path: Optional[Dict[str, str]] = None):

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

    if chain_to_template_path is not None:
        new_data = []
        for data_dict in data:
            name = data_dict['name']

            # delete components associated with chains
            new_components = []
            templates_used = set()
            for component in data_dict['components']:
                if component['chain_id'] not in chain_to_template_path:
                    new_components.append(component)
                else:
                    template_path = chain_to_template_path[component['chain_id']]
                    if template_path not in templates_used:
                        new_components.append({"path": template_path})
                        templates_used.add(template_path)
                
            # save template chains
            template_selection = list(chain_to_template_path.keys())

            new_data.append({
                "name": name,
                "components": new_components,
                "template_selection": template_selection,
            })

        data = new_data


    with open(outpath, 'w+') as f:
        json.dump(data, f, indent=4)




if __name__ == '__main__':

    designed_chain = "D"

    use_templating = True

    ## fixed across designs to the same on-target
    chain_to_template_path = {
        'A': '/gscratch/stf/gvisan01/foundry/on_off_target/input_pdbs/5brz_mhc.pdb',
        'B': '/gscratch/stf/gvisan01/foundry/on_off_target/input_pdbs/5brz_mhc.pdb',
    }
    to_align_chains = list(chain_to_template_path.keys())
    designed_chains = [designed_chain]

    num_proteinmpnn_samples = 10
    proteinmpnn_temperature = 0.2
    proteinmpnn_temp_dir = '/gscratch/stf/gvisan01/foundry/on_off_target/proteinmpnn_temp_dir/'


    ## variable designs
    design_spec = 'magea3_fix_bb_fix_side' # magea3_fix_bb_fix_side, magea3_fix_bb_flex_side


    # make input json, and also save paths to designed structure to compute metrics later
    input_dir = f'./inference_outputs/{design_spec}/'
    out_dir = f'./folding_output/{design_spec}/'

    os.makedirs(out_dir, exist_ok=True)
    name_to_chain_to_seq = {}
    name_to_designed_structure = {}
    for design_folder in os.listdir(input_dir): # list all the .cif files
        for filename in os.listdir(os.path.join(input_dir, design_folder)):
            if filename.endswith('.cif'):
                design_name = filename[:-4] # remove .cif
                designed_structure = os.path.join(input_dir, design_folder, filename)
                folding_name = f'{design_folder}__{design_name}__templating_is_{use_templating}'

                chain_to_seq = extract_sequence_from_structure(designed_structure)

                # compute ends of designed chain in the global sequence, to slice proteinmpnn's output
                chains = sorted(list(chain_to_seq.keys()))
                start_i = 0
                for chain in chains:
                    if chain != designed_chain:
                        start_i += len(chain_to_seq[chain])
                    else:
                        end_i = start_i + len(chain_to_seq[chain])
                        break

                # sample protein sequences with proteinmpnn, add to name_to_chain_to_seq with new name for each sample, and also save designed structure path for each sample to compute metrics later
                sampled_sequences = sample_sequences_with_proteinmpnn(designed_structure, designed_chain, [start_i, end_i], num_proteinmpnn_samples, proteinmpnn_temperature, proteinmpnn_temp_dir)

                name_to_chain_to_seq[f'{folding_name}__seq_rfd3'] = chain_to_seq
                name_to_designed_structure[f'{folding_name}__seq_rfd3'] = designed_structure

                for i, sampled_seq in enumerate(sampled_sequences):
                    name_to_chain_to_seq[f'{folding_name}__seq_mpnn_{i}'] = {**chain_to_seq, designed_chain: sampled_seq}
                    name_to_designed_structure[f'{folding_name}__seq_mpnn_{i}'] = designed_structure

    input_json = os.path.join(out_dir, f'{design_spec}__inputs.json')
    make_rf3_input_json(name_to_chain_to_seq, input_json, chain_to_template_path if use_templating else None)

    name_to_designed_structure_json = os.path.join(out_dir, f'{design_spec}__name_to_designed_structure.json')
    with open(name_to_designed_structure_json, 'w+') as f:
        json.dump(name_to_designed_structure, f, indent=4)
    
    print('-'* 40)
    print('\tsaved input json to: ' + input_json)
    print('-'* 40)

    ## run folding
    folding_command = f"rf3 fold ckpt_path='/gscratch/stf/gvisan01/_foundry_models/rf3_foundry_01_24_latest_remapped.ckpt' inputs='{input_json}' out_dir='{out_dir}' num_steps=50"
    print('-'* 40)
    print('Running folding command:')
    print('\t' + folding_command)
    print('-'* 40)
    os.system(folding_command)

    ## compute design consistency metrics
    print('Computing design consistency metrics...')

    with open(name_to_designed_structure_json, 'r') as f:
        name_to_designed_structure = json.load(f)


    for folding_name, designed_structure in tqdm(name_to_designed_structure.items()):

        folded_structure = os.path.join(out_dir, folding_name, f'{folding_name}_model.cif')

        if not os.path.exists(folded_structure):
            print(f'Folded structure not found for {folding_name}, skipping metric computation for this design.')
            continue

        ref_rmsd = refolding_rmsd(designed_structure, folded_structure, to_align_chains, [designed_chain]) ## TODO this might be wrong...

        with open(os.path.join(out_dir, folding_name, f'{folding_name}_consistency_metrics.json'), 'w+') as f:
            json.dump(
                {
                    "refolding_rmsd_of_designed_chains": ref_rmsd,
                },
                f, indent = 4
            )


