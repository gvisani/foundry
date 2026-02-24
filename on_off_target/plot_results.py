
import os
from tqdm import tqdm
import numpy as np
import json
import matplotlib.pyplot as plt


if __name__ == '__main__':

    ## plot distribution of refolding RMSD, ipTM, and binder-peptide PAE

    for design_spec in ['test_with_alt', 'test_with_alt__alt']: # ['magea3_fix_bb_fix_side', 'magea3_fix_bb_fix_side__alt']:
        
        input_dir = f'./folding/{design_spec}/'

        ref_rmsd_list_rfd3 = []
        iptm_list_rfd3 = []
        pep_binder_min_pae_list_rfd3 = []

        ref_rmsd_list = []
        iptm_list = []
        pep_binder_min_pae_list = []

        for pred_dir in tqdm(os.listdir(input_dir)):
            pred_name = pred_dir
            pred_dir = os.path.join(input_dir, pred_dir)

            if os.path.isdir(pred_dir):

                if os.path.exists(os.path.join(pred_dir, f'{pred_name}_consistency_metrics.json')):
                    with open(os.path.join(pred_dir, f'{pred_name}_consistency_metrics.json'), 'r') as f:
                        metrics = json.load(f)
                        ref_rmsd = metrics['refolding_rmsd_of_designed_chains']
                else:
                    ref_rmsd = None

                with open(os.path.join(pred_dir, f'{pred_name}_summary_confidences.json'), 'r') as f:
                    metrics = json.load(f)
                    iptm = metrics['iptm']
                    pep_binder_min_pae = metrics['chain_pair_pae_min'][2][3]
            
                if pred_name.endswith('rfd3'):
                    if ref_rmsd is not None:
                        ref_rmsd_list_rfd3.append(ref_rmsd)
                    iptm_list_rfd3.append(iptm)
                    pep_binder_min_pae_list_rfd3.append(pep_binder_min_pae)
                else:
                    if ref_rmsd is not None:
                        ref_rmsd_list.append(ref_rmsd)
                    iptm_list.append(iptm)
                    pep_binder_min_pae_list.append(pep_binder_min_pae)
        

        fig, axs = plt.subplots(figsize=(13, 4), ncols=3, nrows=1,sharey=True)

        ax = axs[0]
        ax.hist(ref_rmsd_list)
        ax.hist(ref_rmsd_list_rfd3)
        ax.set_title(design_spec)
        ax.set_ylabel('Count')
        ax.set_xlabel('Refolding RMSD of binder\naligned along the MHC')

        ax = axs[1]
        ax.hist(iptm_list)
        ax.hist(iptm_list_rfd3)
        ax.set_title(design_spec)
        ax.set_ylabel('Count')
        ax.set_xlabel('ipTM')

        ax = axs[2]
        ax.hist(pep_binder_min_pae_list, label='RFD3+pMPNN')
        ax.hist(pep_binder_min_pae_list_rfd3, label='RFD3')
        ax.set_title(design_spec)
        ax.set_ylabel('Count')
        ax.set_xlabel('Min PAE between peptide and binder')
        ax.legend()

        plt.tight_layout()
        plt.savefig(f'{design_spec}.png')
        plt.close()
