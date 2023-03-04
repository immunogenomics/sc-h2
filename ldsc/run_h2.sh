# run h2 checks for traits of interest

source activate ldsc
cd /data/srlab/agupta/code/ldsc
date='083122'
annot_base_folder='/data/srlab/agupta/data/h2'
outfolder='h2_output/09_2022_traits/'

ld_chr_path="/data/srlab/external-data/LDSCORE/data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
baselineLD_annot="/data/srlab/amariuta/jobs/ClassifierRegMap/runLDSC/BL/customized_baselineLD_cts."
frqfile_path="/data/srlab/external-data/LDSCORE/data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq/1000G.EUR.QC."

queue='short'
cpus=1

while IFS="," read -r trait fname remainder; do

    h2_path="/data/srlab/agupta/data/all_PASS_traits/${fname}.sumstats"
    trait_output_name=${trait}
    
    for cell_type in T myeloid B fibroblast endothelial
    do
    
        D_annot="$annot_base_folder/annotations/d12/${cell_type}/${date}_${cell_type}_DYNAMIC_hg19."
        I_annot="$annot_base_folder/annotations/d12/${cell_type}/${date}_${cell_type}_INVARIANT_hg19."

        out_file="$annot_base_folder/$outfolder/072622.${trait_output_name}.${cell_type}_DI:COMBINED"

        command="python ldsc.py --h2 $h2_path --w-ld-chr $ld_chr_path --ref-ld-chr $D_annot,$I_annot,$baselineLD_annot --overlap-annot --frqfile-chr $frqfile_path --out $out_file --print-coefficients"

        bsub -J "${trait_output_name}_${annot}" -q $queue -n $cpus -R 'select[hname!=cn001]' -R 'select[hname!=cn002]' -R 'select[hname!=cn003]' -R 'select[hname!=cn004]' -R 'select[hname!=cn005]' -e /data/srlab/agupta/code/%$J.err -o /data/srlab/agupta/code/%$J.out $command
        
    done
done < "/data/srlab/agupta/data/ref_files/AI_diseases_h2g.csv"