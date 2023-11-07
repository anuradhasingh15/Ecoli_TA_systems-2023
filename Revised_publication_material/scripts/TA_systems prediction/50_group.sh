#!/bin/bash
#PBS -N 50_group
#PBS -S /bin/bash
#PBS -q workq
#PBS -l nodes=1:ppn=40
#PBS -j oe
#PBS -o /home/nasls/ANU/phd_project/TA-SYSTEMS/anaylsis-re/sling/qsub_logfiles/50_group.txt
cd $PBS_O_WORKDIR
np=`wc -l < $PBS_NODEFILE`
/home/nasls/miniconda3/envs/sling/bin/sling group ID_1 toxins -u -it -c 40 -o /home/nasls/ANU/phd_project/TA-SYSTEMS/anaylsis-re/sling/output



