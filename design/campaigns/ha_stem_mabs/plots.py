
import os
import numpy as np
import matplotlib.pyplot as plt
import json
from pathlib import Path



if __name__ == '__main__':

    t1 = '4fqi'
    t2 = '5kan'
    t1_t2 = '4fqi_5kan'

    
    min_distances = {}
    proportion_loop = {}

    for (name, main, alt) in [(t1, t1, t2),
                              (t2, t2, t1),
                              (t1_t2, t1, t2),
                              (t1_t2, t1, t2),
                              ('4fqi_neg05_5kan', t1, t2),
                              ('4fqi_neg07_5kan', t1, t2),
                              ('4fqi_neg10_5kan', t1, t2),
                              ('4fqi_neg15_5kan', t1, t2),
                              ('4fqi_neg25_5kan', t1, t2),
                              ('4fqi_neg50_5kan', t1, t2)]:

        min_distances[name] = {main: [], alt: []}
        proportion_loop[name] = []

        for filename in os.listdir(f"./diffusion/{name}"):
            if filename.endswith('.json'):
                try:
                    with open(f"./diffusion/{name}/{filename}", 'r') as f:
                        metrics = json.load(f)["metrics"]
                except json.decoder.JSONDecodeError:
                    continue
                if filename.endswith('_alt.json'):
                    min_distances[name][alt].append(metrics["ca_distance_of_binder_to_fixed_seq"]['min'])
                else:
                    min_distances[name][main].append(metrics["ca_distance_of_binder_to_fixed_seq"]['min'])
                    proportion_loop[name].append(metrics["loop_fraction"])
                    if '_' not in name:
                        min_distances[name][alt].append(metrics["ca_distance_of_binder_to_fixed_seq_of_aligned_alt_target"]['min'])


    fontsize = 15

    ncols = 3
    nrows = 1
    fig, axs = plt.subplots(figsize=(ncols*3.5, nrows*4.0), ncols=ncols, nrows=nrows, sharex=True, sharey=True)

    ax = axs[0]
    ax.hist(min_distances['4fqi']['4fqi'], alpha=0.4, color='blue', label='target = 4fqi')
    ax.hist(min_distances['4fqi']['5kan'], alpha=0.4, color='red', label='target = 5kan')
    ax.set_title('diffusion target\n4fqi', fontsize=fontsize)
    ax.set_xlabel('min CA-CA dist\nbinder-target (Ang)', fontsize=fontsize)
    ax.set_ylabel('count', fontsize=fontsize)
    ax.tick_params(labelsize=fontsize-2)
    ax.legend(fontsize=fontsize-2)

    ax = axs[1]
    ax.hist(min_distances['5kan']['4fqi'], alpha=0.4, color='blue', label='target = 4fqi')
    ax.hist(min_distances['5kan']['5kan'], alpha=0.4, color='red', label='target = 5kan')
    ax.set_title('diffusion target\n5kan', fontsize=fontsize)
    ax.set_xlabel('min CA-CA dist\nbinder-target (Ang)', fontsize=fontsize)
    ax.tick_params(labelsize=fontsize-2)
    ax.legend(fontsize=fontsize-2)

    ax = axs[2]
    ax.hist(min_distances['4fqi_5kan']['4fqi'], alpha=0.4, color='blue', label='target = 4fqi')
    ax.hist(min_distances['4fqi_5kan']['5kan'], alpha=0.4, color='red', label='target = 5kan')
    ax.set_title('diffusion target\n4fqi_5kan', fontsize=fontsize)
    ax.set_xlabel('min CA-CA dist\nbinder-target (Ang)', fontsize=fontsize)
    ax.tick_params(labelsize=fontsize-2)
    ax.legend(fontsize=fontsize-2)

    plt.tight_layout()
    plt.savefig('min_dist_to_targets.png', dpi=300)
    plt.close()


    w_list = ['05', '07', '10', '15', '25'] #, '50']

    ncols = len(w_list)
    nrows = 1
    fig, axs = plt.subplots(figsize=(ncols*3.5, nrows*4.0), ncols=ncols, nrows=nrows, sharex=True, sharey=True)

    for i, w in enumerate(w_list):

        name = f'4fqi_neg{w}_5kan'

        ax = axs[i]
        ax.hist(min_distances[name]['4fqi'], alpha=0.4, color='blue', label='target = 4fqi')
        ax.hist(min_distances[name]['5kan'], alpha=0.4, color='red', label='target = 5kan')
        ax.set_title(f'diffusion target\n{name}', fontsize=fontsize)
        ax.set_xlabel('min CA-CA dist\nbinder-target (Ang)', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize-2)
        ax.legend(fontsize=fontsize-2)

    plt.tight_layout()
    plt.savefig('min_dist_to_targets_neg.png', dpi=300)
    plt.close()


    w_list = ['05', '07', '10', '15', '25', '50']

    ncols = len(w_list)
    nrows = 1
    fig, axs = plt.subplots(figsize=(ncols*3.5, nrows*4.0), ncols=ncols, nrows=nrows, sharex=True, sharey=True)

    for i, w in enumerate(w_list):

        name = f'4fqi_neg{w}_5kan'

        ax = axs[i]
        ax.hist(min_distances[name]['4fqi'], alpha=0.4, color='blue', label='target = 4fqi')
        ax.hist(min_distances[name]['5kan'], alpha=0.4, color='red', label='target = 5kan')
        ax.set_title(f'diffusion target\n{name}', fontsize=fontsize)
        ax.set_xlabel('min CA-CA dist\nbinder-target (Ang)', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize-2)
        ax.legend(fontsize=fontsize-2)

    plt.tight_layout()
    plt.savefig('min_dist_to_targets_neg_with_50.png', dpi=300)
    plt.close()


    ## plot loopiness
    names = ['4fqi', '5kan', '4fqi_5kan', '4fqi_neg05_5kan', '4fqi_neg07_5kan', '4fqi_neg10_5kan', '4fqi_neg15_5kan', '4fqi_neg25_5kan', '4fqi_neg50_5kan']
    
    ncols = len(names)
    nrows = 1
    fig, axs = plt.subplots(figsize=(ncols*3.5, nrows*4.0), ncols=ncols, nrows=nrows, sharex=True, sharey=True)

    for i, name in enumerate(names):

        ax = axs[i]
        ax.hist(proportion_loop[name], alpha=0.4, color='purple')
        ax.set_title(f'diffusion target\n{name}', fontsize=fontsize)
        ax.set_xlabel('proportion of loops in binder', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize-2)

    plt.tight_layout()
    plt.savefig('proportion_loops.png', dpi=300)
    plt.close()


    ## plot distribution of refolding RMSD, ipTM, and binder-peptide PAE

    design_spec_to_rf3_metrics = {}

    for design_spec in ['4fqi_5kan', '4fqi_5kan__alt', '4fqi', '4fqi__alt', '5kan', '5kan__alt']:
        
        input_dir = f'./folding/{design_spec}/'

        ref_rmsd_list_rfd3 = []
        iptm_list_rfd3 = []

        ref_rmsd_list = []
        iptm_list = []

        for pred_dir in os.listdir(input_dir):
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
            
                if pred_name.endswith('rfd3'):
                    if ref_rmsd is not None:
                        ref_rmsd_list_rfd3.append(ref_rmsd)
                    iptm_list_rfd3.append(iptm)
                else:
                    if ref_rmsd is not None:
                        ref_rmsd_list.append(ref_rmsd)
                    iptm_list.append(iptm)
        
        design_spec_to_rf3_metrics[design_spec] = {
            'iptm': iptm_list,
            'iptm_rfd3': iptm_list_rfd3,
            'ref_rmsd': ref_rmsd_list,
            'ref_rmsd_rfd3': ref_rmsd_list_rfd3
        }
    
    for design_spec in ['4fqi_5kan', '4fqi', '5kan']:
        
        fig, axs = plt.subplots(figsize=(9, 4), ncols=2, nrows=1, sharey=True)
        fontsize = 14

        if design_spec == '5kan':
            main_color = 'red'
            alt_color = 'blue'
            main_target = '5kan'
            alt_target = '4fqi'
        else:
            main_color = 'blue'
            alt_color = 'red'
            main_target = '4fqi'
            alt_target = '5kan'

        metrics = design_spec_to_rf3_metrics[design_spec]
        metrics_alt = design_spec_to_rf3_metrics[design_spec + '__alt']

        ax = axs[0]
        ax.hist(metrics['ref_rmsd'], color=main_color, alpha=0.4)
        ax.hist(metrics_alt['ref_rmsd'], color=alt_color, alpha=0.4)
        ax.set_title(design_spec, fontsize=fontsize)
        ax.set_ylabel('Count', fontsize=fontsize)
        ax.set_xlabel('Refolding RMSD of binder\naligned along the target', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize-2)

        ax = axs[1]
        ax.hist(metrics['iptm'], color=main_color, alpha=0.4, label=f'target = {main_target}')
        ax.hist(metrics_alt['iptm'], color=alt_color, alpha=0.4, label=f'target = {alt_target}')
        ax.set_title(design_spec, fontsize=fontsize)
        ax.set_xlabel('ipTM', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize-2)
        ax.legend()

        plt.tight_layout()
        plt.savefig(f'{design_spec}.png')
        plt.close()


