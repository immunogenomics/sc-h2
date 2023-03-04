##!/bin/bash

date='083122'
queue='short'

for cell_type in T myeloid B fibroblast endothelial
do
    for annot_type in OPEN DYNAMIC INVARIANT
    do
        command="python bed_to_annot.py --cell_type ${cell_type} --annot ${annot_type}  --date ${date}"
            
        bsub -J ${cell_type}${annot_type}[1-10] -q $queue -R 'select[hname!=cn001]' -R 'select[hname!=cn002]' -R 'select[hname!=cn003]' -R 'select[hname!=cn004]' -R 'select[hname!=cn005]' -e $%J.err $command
    done
done