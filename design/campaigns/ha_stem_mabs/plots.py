
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


