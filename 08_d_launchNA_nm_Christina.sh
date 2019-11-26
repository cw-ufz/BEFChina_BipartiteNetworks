#!/bin/bash

#  launchNA_Christina.sh
#  
#
#  Created by Anna Heintz-Buschart on 04.07.19.
#

for div in 1 2 3
do
for lm in p20 p40 p60 p80 p100
do

input_real=/gpfs1/work/$USER/Nets/input/05_tree_combinations_samples.RDS
input_nm=/gpfs1/work/$USER/Nets/input/05_tree_combinations_samples_nm.RDS

shortinput=$(basename $input_real .RDS)

#CMD="qsub -cwd -N na.$div.$lm.$shortinput /gpfs1/work/$USER/Nets/scripts/runNA_Christina.sh $div $lm $input_real $input_nm"
#echo $CMD

qsub -cwd -N na.$div.$lm.$shortinput /gpfs1/work/$USER/Nets/scripts/runNA_nm_Christina.sh $div $lm $input_real $input_nm


done
done
