

import os
import json
from pathlib import Path
from tqdm import tqdm
from omegaconf import OmegaConf
import argparse
from atomworks.io.utils.io_utils import load_any
from metrics import distance_statistics



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    args = parser.parse_args()

    cfg = OmegaConf.load(args.config)

    ## CA-CA distances between designed chain and hotspot
    ## CA-CA distances between designed chain and full other chains

    designed_chain = cfg.designed_chain
    other_chains = list(cfg.alternative_target.keys()) # TODO more robust way of specifying this

    # gather all the names, for which there should be both a .cif file (structure) and a .json file (metrics)
    names = set()
    for filename in os.listdir(cfg.diffusion_dir):
        names.add(Path(filename).stem)

    for name in tqdm(names):
        cif_path = os.path.join(cfg.diffusion_dir, name+'.cif')
        json_path = os.path.join(cfg.diffusion_dir, name+'.json')

        try:
            with open(json_path, 'r') as f_metadata:
                metadata = json.load(f_metadata)
                dist_to_other_chains = distance_statistics(cif_path, designed_chain, other_chains)
                metadata["metrics"] = metadata["metrics"] | {"ca_distance_of_binder_to_fixed_seq": dist_to_other_chains}

                # hotspots, if they exist
                if "select_hotspots" in metadata["specification"]:

                    def convert_chain(chain): # just to deal with current nuisance of chain mismatch... TODO make this more robust
                        if chain == 'H':
                            return 'A'
                        elif chain == 'L':
                            return 'B'
                        return chain
                    
                    hotspots = [(convert_chain(elem[0]), int(elem[1:])) for elem in metadata["specification"]["select_hotspots"].keys()]
                    dist_to_hotspots = distance_statistics(cif_path, designed_chain, hotspots)
                    metadata["metrics"] = metadata["metrics"] | {"ca_distance_of_binder_to_hotspots": dist_to_hotspots}

            with open(json_path, 'w+') as f_metadata:
                json.dump(metadata, f_metadata, indent=4)
        except json.decoder.JSONDecodeError:
            print(f'json decode error with {name}')
            continue

    

