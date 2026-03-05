

config=./configs/4fqi_neg50_5kan.yaml # 4fqi_5kan


python ../../src/rfd3_design.py --config $config


python ../../src/mpnn_sample_and_prep_rf3_inputs.py --config $config

## The RF3 step is very slow! Since we are folding num_proteinmpnn_samples many more  omplexes than we are designing
python ../../src/rf3_fold.py --config $config --run_alt_target 0
python ../../src/compute_additional_metrics_rf3.py --config $config --run_alt_target 0

python ../../src/rf3_fold.py --config $config --run_alt_target 1
python ../../src/compute_additional_metrics_rf3.py --config $config --run_alt_target 1
