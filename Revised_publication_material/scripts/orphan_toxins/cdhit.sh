#!/bin/bash
#PBS -N quast
#PBS -S /bin/bash
#PBS -q workq
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -o /home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_toxins/cd-hit_log.txt
cd $PBS_O_WORKDIR
np=`wc -l < $PBS_NODEFILE`
/home/nasls/miniconda3/envs/cd-hit/bin/cd-hit -i all_toxins_set1.fasta -o all_toxins_clustered_set1.fasta -c 0.9 -d 0




