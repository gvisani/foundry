
# config=./configs/magea3_fix_bb_fix_side_cfg.yaml

config=./configs/test_with_alt.yaml


# python rfd3_design.py --config $config

python mpnn_sample_and_prep_rfd3_inputs.py --config $config

## The RF3 step is very slow! Since we are folding num_proteinmpnn_samples many more  omplexes than we are designing
python rf3_fold.py --config $config --run_alt_target 0
python compute_additional_metrics.py --config $config --run_alt_target 0

python rf3_fold.py --config $config --run_alt_target 1
python compute_additional_metrics.py --config $config --run_alt_target 1
