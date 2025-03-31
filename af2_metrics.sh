#!/bin/bash
#$ -S /bin/bash
#$ -m bea
#$ -pe smp 20
#$ -cwd
#$ -l mem_free=1G
#$ -l h_rt=48:00:00

conda activate /wynton/home/kortemme/scrilly/proteindesign/proteindesign

DIR=$(pwd)

python af2_utils.py \
	--colaboutdir $DIR \
	--outdir $DIR \
	--min0_bb_dir "/wynton/home/kortemme/scrilly/helix_sliding/20230710_min_bbs_for_mpnn/min2_53_threads" \
	--single_state  