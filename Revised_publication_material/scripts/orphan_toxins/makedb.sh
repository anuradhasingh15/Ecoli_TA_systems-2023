#!/bin/bash
#PBS -N makedb
#PBS -S /bin/bash
#PBS -q workq
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -o /home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_toxins/makedblog.txt
cd $PBS_O_WORKDIR
np=`wc -l < $PBS_NODEFILE`
makeblastdb -in /home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_toxins/all_toxins_clustered_set1.fasta -dbtype nucl -out /home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_toxins/blast_db/blast_db
