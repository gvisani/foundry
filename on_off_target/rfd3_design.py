
import os
from omegaconf import OmegaConf
import argparse

import os
import gzip
import shutil

def uncompress(filepath: str) -> None:
    '''
    Uncompresses .gz file

    Deletes compressed file afterwards
    '''
    with gzip.open(filepath, 'rb') as f_in:
        with open(filepath[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    os.remove(filepath) # delete compressed file


def uncompress_in_dir(dirpath: str) -> None:
    '''
    Recursively because it's fun :)
    '''
    for file_or_dir in os.listdir(dirpath):

        file_or_dir_path = os.path.join(dirpath, file_or_dir)

        if os.path.isdir(file_or_dir_path):
            uncompress_in_dir(file_or_dir_path)
        else:
            if file_or_dir_path.endswith('.gz'):
                uncompress(file_or_dir_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    args = parser.parse_args()

    cfg = OmegaConf.load(args.config)

    if "alt_scale" in cfg:
        alt_scale = cfg.alt_scale
    else:
        alt_scale = 0.5

    command = f"rfd3 design"
    # command = f"python /gscratch/stf/gvisan01/foundry/models/rfd3/src/rfd3/run_inference.py"
    command += f" out_dir={cfg.diffusion_dir}"
    command += f" ckpt_path=/gscratch/stf/gvisan01/_foundry_models/rfd3_latest.ckpt"
    command += f" inputs={cfg.rfd3_json_path}"
    command += f" n_batches={cfg.n_diffusion_batches_of_8}"
    command += f" inference_sampler.step_scale=3"
    command += f" inference_sampler.gamma_0=0.2"
    command += f" inference_sampler.alt_scale={alt_scale}"
    command += f" inference_sampler.use_classifier_free_guidance=False" # CFG adds the unconditional velocity!

    print('-'* 40)
    print('Running RFD3 design with command:')
    print('\t' + command)
    print('-'* 40)

    os.system(command)

    uncompress_in_dir(cfg.diffusion_dir)

