##!/bin/bash

# (1) regression scores for each peak -> bed files

date='083122'
peak_window=500

for cell_type in T B myeloid fibroblast endothelial 
do
    command="python scores_to_bed.py --cell_type ${cell_type} --date $date --peak_window $peak_window"

    bsub -J ${cell_type}[1-10] -q vshort -R 'select[hname!=cn001]' -R 'select[hname!=cn002]' -R 'select[hname!=cn003]' -R 'select[hname!=cn004]' -R 'select[hname!=cn005]' -e $%J.err $command
done

##################################################

# (2) bed version hg38 -> hg19

module load liftover
date='083122_'
inp_path="/data/srlab/agupta/data/h2/bed/d12/hg38/"
out_path="/data/srlab/agupta/data/h2/bed/d12/hg19/"
chain_path="/data/srlab/agupta/code/stat_gen/hg38ToHg19.over.chain"

for cell_type in T B myeloid fibroblast endothelial
do
    for annot_type in OPEN DYNAMIC INVARIANT
    do
        command="liftOver $inp_path$date${cell_type}_${annot_type}_hg38.bed $chain_path $out_path$date${cell_type}_${annot_type}_hg19.bed $inp_path$date${cell_type}_${annot_type}_UNMAPPED_hg38.bed"

        bsub -J ${cell_type}${annot_type}[1-10] -q vshort -R 'select[hname!=cn001]' -R 'select[hname!=cn002]' -R 'select[hname!=cn003]' -R 'select[hname!=cn004]' -R 'select[hname!=cn005]' -e $%J.err $command
    done
done