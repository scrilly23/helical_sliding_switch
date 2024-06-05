#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q gpu.q
#$ -pe smp 1
#$ -l mem_free=10G
#$ -l h_rt=01:00:00
#$ -t 1-3070
#$ -l h='qb3-atgpu*'
#$ -l gpu_mem=10000M

#global stuff
TASK_ID=$((${SGE_TASK_ID}-1))
#TASK_ID=$((1-1))

conda activate /wynton/home/kortemme/scrilly/localcolabfold/colabfold-conda

export CUDA_VISIBLE_DEVICES=$SGE_GPU
#module purge
module load Sali cuda

#INPUTS
#total number of files to loop through 
NUM_JOBS=3070
END=3070
FASTA_DIR="/wynton/home/kortemme/scrilly/helix_sliding/20230822_0321_min0_52_ssd_sep_fa"
OUTPUT_FOLDER="/wynton/home/kortemme/scrilly/helix_sliding/20230827_min0_52_ssd_colabfold_array"
#############################

#make output dir if doesn't exist

if [ ! -d $OUTPUT_FOLDER]; then
  mkdir -p $OUTPUT_FOLDER
fi

INPUT_FILES=($FASTA_DIR/*.fa)

INPUTFILENAME="${INPUT_FILES[$SGE_TASK_ID - 1]}"
echo $INPUTFILENAME

colabfold_batch \
    --num-recycle 3 \
    --num-models 5 \
    --msa-mode single_sequence \
    --rank plddt \
    --use-dropout \
    --model-type alphafold2_multimer_v3 \
    --disable-unified-memory \
    $INPUTFILENAME \
    $OUTPUT_FOLDER



