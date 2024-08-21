#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -m bea
#$ -l mem_free=1G
#$ -l h_rt=4:00:00

DIR=$(pwd)
PDB_HEADER="13632_"
FILE_HEADER="MSD_13632_07144"

python socket_utils.py \
	--indir $DIR \
	--pdbheader $PDB_HEADER \
	--outdir $DIR \
	--fileheader $FILE_HEADER
