#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -m bea
#$ -j y
#$ -pe smp 4
#$ -l h_vmem=8G
#$ -l h_rt=2:00:00

module load CBI trimgalore

set -euo pipefail

###INPUTS####
INPUT_DIR="${1:-$PWD}"
TRIMMED_DIR="${INPUT_DIR}/trimmed"
MERGED_DIR="${INPUT_DIR}/merged"
DESIGNS_CSV="./r3_hs_all_protein_seqs_final.csv"
####

####TRIM AND MERGE####
mkdir -p "$TRIMMED_DIR" "$MERGED_DIR" logs

mapfile -t R1_FILES < <(find "$INPUT_DIR" -maxdepth 1 -name "*_R1_*.fastq.gz" | sort)

for R1 in "${R1_FILES[@]}"; do
    R2="${R1/_R1_/_R2_}"
    [[ ! -f "$R2" ]] && echo "WARNING: No R2 found for $R1, skipping." && continue

    SAMPLE=$(basename "$R1" | sed 's/_R1_.*//')

    trim_galore \
        --quality 20 \
        --fastqc \
        --gzip \
        --cores 4 \
        -o "$TRIMMED_DIR" \
        --paired "$R1" "$R2"

    TRIMMED_R1=$(find "$TRIMMED_DIR" -maxdepth 1 -name "${SAMPLE}*_val_1.fq.gz" | head -1)
    TRIMMED_R2=$(find "$TRIMMED_DIR" -maxdepth 1 -name "${SAMPLE}*_val_2.fq.gz" | head -1)

    /wynton/home/kortemme/scrilly/NGmerge-master/NGmerge \
        -1 "$TRIMMED_R1" \
        -2 "$TRIMMED_R2" \
        -o "${MERGED_DIR}/${SAMPLE}" \
        -n 4
done

####

####EXTRACT MATCHING READS####
conda activate /wynton/home/kortemme/scrilly/proteindesign/proteindesign

python extract_matching_reads.py "$MERGED_DIR" "${INPUT_DIR}/mapped_seqs_outdir" --save_seqs --designs "$DESIGNS_CSV"

####
