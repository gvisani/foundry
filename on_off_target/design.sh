
identifier=magea3_fix_bb_fix_side  #  magea3_fix_bb_fix_side  magea3_fix_bb_flex_side


output_dir=./inference_outputs

if [[ -d "$output_dir/$identifier" ]] && [[ -n "$(ls -A "$output_dir/$identifier" 2>/dev/null)" ]]; then
    # Find max index among existing directories
    max_idx=$(ls -1 "$output_dir/$identifier" | grep -E '^[0-9]+$' | sort -n | tail -1)
    next_idx=$((max_idx + 1))
else
    next_idx=0
fi

rfd3 design \
    out_dir=$output_dir/$identifier/$next_idx  \
    ckpt_path=/gscratch/stf/gvisan01/_foundry_models/rfd3_latest.ckpt \
    inputs=./configs/$identifier.json \
    n_batches=10 \
    inference_sampler.step_scale=3 \
    inference_sampler.gamma_0=0.2

