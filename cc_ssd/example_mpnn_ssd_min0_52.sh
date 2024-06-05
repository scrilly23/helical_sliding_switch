#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/kortemme/scrilly/ProteinMPNN/vanilla_proteinmpnn/my_mpnn
#$ -e /wynton/home/kortemme/scrilly/ProteinMPNN/vanilla_proteinmpnn/my_mpnn
#$ -q gpu.q
#$ -l mem_free=2G
#$ -m bea
#$ -l h_rt=12:00:00

conda activate mlfold_gpu_2

folder_with_pdbs="/wynton/home/kortemme/scrilly/helix_sliding/20230710_min_bbs_for_mpnn/min0_52_threads"

output_dir="/wynton/home/kortemme/scrilly/helix_sliding/20230710_single_state_seqs/min0_52_mpnn_seqs"

if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi


path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_tied_positions=$output_dir"/tied_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_positions.jsonl"

chains_to_design="A B"
#pos_neg_chain_list="A"
#chain_betas="1.0"

fixed_positions=",24 25 26 27 28 29 30 31 32 33 34 35 36"

python /wynton/home/kortemme/scrilly/ProteinMPNN/vanilla_proteinmpnn/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

#python ../helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

#python ../helper_scripts/make_pos_neg_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --homooligomer 1 --pos_neg_chain_list $pos_neg_chain_list --pos_neg_chain_betas $chain_betas

python /wynton/home/kortemme/scrilly/ProteinMPNN/vanilla_proteinmpnn/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python /wynton/home/kortemme/scrilly/ProteinMPNN/vanilla_proteinmpnn/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 10 \
        --sampling_temp "0.3" \
        --batch_size 1 \
        --omit_AAs 'CHX' \
        --save_probs 1
        #--chain_id_jsonl $path_for_assigned_chains \
        #--tied_positions_jsonl $path_for_tied_positions \
