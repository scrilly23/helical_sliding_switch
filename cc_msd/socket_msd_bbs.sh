#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -t 1-5000
#$ -l mem_free=1G
#$ -l h_rt=00:05:00

#global stuff
TASK_ID=$((${SGE_TASK_ID}-1))
#TASK_ID=$((1-1))

#INPUTS
#total number of files to loop through 
NUM_JOBS=5000
END=5000
INDIR="/wynton/home/kortemme/scrilly/helix_sliding/20240517_mpnn_msd_t_01_cf_array"
#############################

INPUT_FILES=($INDIR/*.pdb)
INPUTFILENAME="${INPUT_FILES[$SGE_TASK_ID - 1]}"

NAME=$(basename $INPUTFILENAME .pdb)
PDB_PATH="$NAME.pdb"
DSSP_PATH="$NAME.dssp"

socket2 -f $PDB_PATH -s $DSSP_PATH -c '7.0' -v -l -o ${NAME}_socket.txt -r ${NAME}_socket.rasmol
