

import os
import json
from tqdm import tqdm
from omegaconf import OmegaConf
import argparse
from metrics import refolding_rmsd


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    parser.add_argument('--run_alt_target', type=int, default=0)    
    args = parser.parse_args()

    cfg = OmegaConf.load(args.config)

    if args.run_alt_target:
        out_dir = cfg.folding_alt_target_dir
    else:
        out_dir = cfg.folding_dir
    

    ## compute design consistency metrics if running on the original target used for diffusion design
    if args.run_alt_target:
        print('Computing design consistency metrics...')

        with open(cfg.name_to_designed_structure_json_alt_target, 'r') as f:
            name_to_designed_structure = json.load(f)

        for folding_name, designed_structure in tqdm(name_to_designed_structure.items()):

            folded_structure = os.path.join(out_dir, folding_name, f'{folding_name}_model.cif')

            if not os.path.exists(folded_structure):
                print(f'Folded structure not found for {folding_name}, skipping metric computation for this design.')
                continue

            ref_rmsd = refolding_rmsd(designed_structure, folded_structure, sorted(list(set([sel[0] for sel in cfg.templating.template_selection]))), [cfg.designed_chain]) ## TODO: check validity of this

            with open(os.path.join(out_dir, folding_name, f'{folding_name}_consistency_metrics.json'), 'w+') as f:
                json.dump(
                    {
                        "refolding_rmsd_of_designed_chains": ref_rmsd,
                    },
                    f, indent = 4
                )
    else:
        print('Computing design consistency metrics...')

        with open(cfg.name_to_designed_structure_json, 'r') as f:
            name_to_designed_structure = json.load(f)

        for folding_name, designed_structure in tqdm(name_to_designed_structure.items()):

            folded_structure = os.path.join(out_dir, folding_name, f'{folding_name}_model.cif')

            if not os.path.exists(folded_structure):
                print(f'Folded structure not found for {folding_name}, skipping metric computation for this design.')
                continue

            ref_rmsd = refolding_rmsd(designed_structure, folded_structure, sorted(list(set([sel[0] for sel in cfg.templating.template_selection]))), [cfg.designed_chain]) ## TODO: check validity of this

            with open(os.path.join(out_dir, folding_name, f'{folding_name}_consistency_metrics.json'), 'w+') as f:
                json.dump(
                    {
                        "refolding_rmsd_of_designed_chains": ref_rmsd,
                    },
                    f, indent = 4
                )
