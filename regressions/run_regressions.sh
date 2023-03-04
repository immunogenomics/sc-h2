##!/bin/bash

filename='regressions_d12.py'
threshold=25
peaks_subset_num=1
family='poisson'
date='083122'
queue_type='normal'
cpus=3

for CT in T myeloid B fibroblast endothelial
do
    command="python ${filename} --CT ${CT} --i ${i} --threshold $threshold --peaks_subset_num $peaks_subset_num --family $family --date $date"
    bsub -J ${CT}[1-10] -q $queue_type -R 'select[hname!=cn001]' -R 'select[hname!=cn002]' -R 'select[hname!=cn003]' -R 'select[hname!=cn004]' -R 'select[hname!=cn005]' -n $cpus -e $%J.err $command
done

#vshort = 15min
#short = 1h
#medium = 24h
#normal = 3d