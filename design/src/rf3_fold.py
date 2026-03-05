
import os
import json
from tqdm import tqdm
from omegaconf import OmegaConf
import argparse



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    parser.add_argument('--run_alt_target', type=int, default=0)    
    args = parser.parse_args()

    cfg = OmegaConf.load(args.config)

    if args.run_alt_target:
        input_json = cfg.input_json_for_folding_alt_target
        out_dir = cfg.folding_alt_target_dir
    else:
        input_json = cfg.input_json_for_folding
        out_dir = cfg.folding_dir


    ## run folding
    folding_command = f"rf3 fold ckpt_path='/gscratch/stf/gvisan01/_foundry_models/rf3_foundry_01_24_latest_remapped.ckpt' inputs='{input_json}' out_dir='{out_dir}' num_steps=50"
    print('-'* 40)
    print('Running folding command:')
    print('\t' + folding_command)
    print('-'* 40)
    os.system(folding_command)


