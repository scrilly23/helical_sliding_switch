#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mem_free=2G
#$ -l h_rt=8:00:00

module load Sali dssp

for filename in /wynton/home/kortemme/scrilly/helix_sliding/20240517_mpnn_msd_t_01_cf_array/*.pdb; do
    echo "$filename"
    mkdssp "$filename"  "/wynton/home/kortemme/scrilly/helix_sliding/20240517_mpnn_msd_t_01_cf_array/$(basename "$filename" .pdb).dssp"
done
