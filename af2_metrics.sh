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
	--outdir $DIR