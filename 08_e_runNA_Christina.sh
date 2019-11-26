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
INPUT="$3"

module load R/3.4.3-1

echo "network analysis task $SGE_TASK_ID ..."
Rscript /work/$USER/Nets/scripts/08_HPC_bipartite_magic7.R $DIV $LINK_MODE $INPUT
echo "done"
