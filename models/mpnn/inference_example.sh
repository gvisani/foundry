

## the main difference between this new implementation and the original one is that,
# if you're asking to design a specific chain, this implementation will still output the *whole* protein sequence for all chains,
# and without even providing chain breaks for ease of parsing


outdir='./example'
mkdir -p $outdir

python src/mpnn/inference.py \
        --model_type protein_mpnn \
        --checkpoint_path /gscratch/stf/gvisan01/foundry/models/mpnn/proteinmpnn_v_48_020.pt \
        --is_legacy_weights True \
        --out_directory $outdir \
        --write_fasta True \
        --write_structures False \
        --structure_path /gscratch/stf/gvisan01/ProteinMPNN-finetuning-TCRpMHC/data/pdbs/1ao7.pdb \
        --temperature 0.5 \
        --batch_size 1 \
        --number_of_batches 10 \
        --designed_chains B

## TODO add some code here that only keeps chain B in the output file

