#!/bin/bash
 
#$ -S /bin/bash
#$ -N networkAnalysis
 
#$ -o /work/$USER/Nets/logs/$JOB_NAME-$JOB_ID.out
#$ -e /work/$USER/Nets/logs/$JOB_NAME-$JOB_ID.err
#$ -m beas
#$ -M christina.lachmann@ufz.de
#$ -l h_rt=144:00:00
#$ -l h_vmem=8G
# $ -l highmem
#$ -binding linear:1

DIV="$1"
LINK_MODE="$2"
INPUT_REAL="$3"
INPUT_NM="$4"

module load R/3.5.1-2

echo "network analysis task $SGE_TASK_ID ..."
Rscript /work/$USER/Nets/scripts/08b_HPC_bipartite_magic7_nestedness_with_nm.R $DIV $LINK_MODE $INPUT_REAL $INPUT_NM
echo "done"
