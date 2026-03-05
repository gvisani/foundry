

import os
import json
import numpy as np
from pathlib import Path
from tqdm import tqdm
from omegaconf import OmegaConf
import argparse
from atomworks.io.utils.io_utils import load_any

import sys
sys.path.append('../../src/')
from metrics import distance_statistics_from_coords
from utils import make_mask, rigid_align_params, apply_rigid_align

import biotite.structure.io.pdbx as pdbx



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    args = parser.parse_args()

    cfg = OmegaConf.load(args.config)

    ## CA-CA distances between designed chain and full other chains **of aligned alt target**

    if cfg.design_spec == '4fqi':
        main_target_pdbfile = '/gscratch/stf/gvisan01/foundry/design/campaigns/ha_stem_mabs/input_pdbs/4fqi_ab_align_renum_filled_renamed.pdb'
        alt_target_pdbfile = '/gscratch/stf/gvisan01/foundry/design/campaigns/ha_stem_mabs/input_pdbs/5kan_ab_align_renum_renamed.pdb'
    elif cfg.design_spec == '5kan':
        main_target_pdbfile = '/gscratch/stf/gvisan01/foundry/design/campaigns/ha_stem_mabs/input_pdbs/5kan_ab_align_renum_renamed.pdb'
        alt_target_pdbfile = '/gscratch/stf/gvisan01/foundry/design/campaigns/ha_stem_mabs/input_pdbs/4fqi_ab_align_renum_filled_renamed.pdb'
    else:
        raise ValueError(f"invalid design_spec {cfg.design_spec}")

    main_target_struc = load_any(main_target_pdbfile, extra_fields=None)[0]
    main_target_struc = main_target_struc[main_target_struc.atom_name == "CA"]

    alt_target_struc_orig = load_any(alt_target_pdbfile, extra_fields=None)[0]
    alt_target_struc_orig = alt_target_struc_orig[alt_target_struc_orig.atom_name == "CA"]

    designed_chain = cfg.designed_chain
    other_chains = list(cfg.alternative_target.keys()) # TODO more robust way of specifying this

    mask_alt_target = make_mask(alt_target_struc_orig, other_chains)

    # gather all the names, for which there should be both a .cif file (structure) and a .json file (metrics)
    names = set()
    for filename in os.listdir(cfg.diffusion_dir):
        names.add(Path(filename).stem)
    names = list(names)

    # names = ['rfd3__4fqi_4fqi_5kan_0_model_0']

    for name in tqdm(names):
        cif_path = os.path.join(cfg.diffusion_dir, name+'.cif')
        json_path = os.path.join(cfg.diffusion_dir, name+'.json')

        try:
            with open(json_path, 'r') as f_metadata:
                metadata = json.load(f_metadata)

                design_struc = load_any(cif_path, extra_fields="all")[0]
                design_struc = design_struc[design_struc.atom_name == "CA"]

                # align alt_target_struc to design_struc based upon main_target_struc, to which it is already aligned
                params = rigid_align_params(main_target_struc, design_struc, make_mask(main_target_struc, other_chains), make_mask(design_struc, other_chains))
                alt_target_struc = apply_rigid_align(alt_target_struc_orig, *params)

                # ## outputting this structure to check it for mistakes
                # cif_file = pdbx.CIFFile()
                # pdbx.set_structure(cif_file, alt_target_struc)
                # cif_file.write(os.path.join(cfg.diffusion_dir, f"{name}__alt_target_aligned_like_main_with_design.cif"))

                mask_design = make_mask(design_struc, designed_chain)

                dist_to_other_chains = distance_statistics_from_coords(
                    design_struc[mask_design].coord, alt_target_struc[mask_alt_target].coord
                )

                metadata["metrics"] = metadata["metrics"] | {"ca_distance_of_binder_to_fixed_seq_of_aligned_alt_target": dist_to_other_chains}

            with open(json_path, 'w+') as f_metadata:
                json.dump(metadata, f_metadata, indent=4)
            
        except json.decoder.JSONDecodeError:
            print(f'json decode error with {name}')
            continue

    

