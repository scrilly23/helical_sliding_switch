#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 100
#$ -l mem_free=2G
#$ -l h_rt=12:00:00

for filename in /wynton/home/kortemme/scrilly/helix_sliding/20230321_bgs_sample_z1off_omega1_idealcc_valine/*.pdb; do
    VAR=basename $filename .pdb
    for filename2 in /wynton/home/kortemme/scrilly/helix_sliding/20230620_dssp_0321_bbs/dssp_files/*.dssp; do
        VAR2=basename $filename2 .dssp
        if [[ $VAR == *$VAR2* ]]; then
            socket2 -f $filename -s $filename2 -c '8.0' -v -l -o "/wynton/home/kortemme/scrilly/helix_sliding/20230620_dssp_0321_bbs/socket_output/$(basename "$filename" .pdb)_socket.txt" -r "/wynton/home/kortemme/scrilly/helix_sliding/20230620_dssp_0321_bbs/socket_rasmol_output/$(basename "$filename" .pdb)_socket.rasmol"
        fi
    done
done