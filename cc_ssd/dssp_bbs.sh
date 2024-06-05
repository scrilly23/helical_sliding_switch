#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l mem_free=2G
#$ -l h_rt=8:00:00

module load Sali dssp

for filename in /wynton/home/kortemme/scrilly/helix_sliding/20230321_bgs_sample_z1off_omega1_idealcc_valine/*.pdb; do
    echo "$filename"
    mkdssp "$filename"  "/wynton/home/kortemme/scrilly/helix_sliding/20230620_dssp_0321_bbs/dssp_files/$(basename "$filename" .pdb).dssp"
done
