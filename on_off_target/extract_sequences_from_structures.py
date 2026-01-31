
from atomworks.io import parse
import json

from typing import Dict


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


def make_rf3_input_json(name_to_chain_to_seq: Dict[str, Dict[str, str]], outpath: str):

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

    with open(outpath, 'w+') as f:
        json.dump(data, f, indent=4)



if __name__ == '__main__':

    chain_to_seq = extract_sequence_from_structure('/gscratch/stf/gvisan01/foundry/on_off_target/inference_outputs/magea3_fix_bb_fix_side/1/magea3_fix_bb_fix_side_magea3_0_model_3.cif')
    
    make_rf3_input_json({'magea3_fix_bb_fix_side_magea3_0_model_3': chain_to_seq}, 'test.json')



