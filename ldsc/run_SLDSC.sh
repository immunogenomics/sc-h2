## run S-LDSC to generate the supporting 4 files for each .annot.gz file

cd /data/srlab/agupta/code/ldsc/
source activate ldsc
queue='short'
cpus=1
date='083122_'
basepath="/data/srlab/agupta/data/h2/annotations/d12"
bfile="/data/srlab/external-data/LDSCORE/data.broadinstitute.org/alkesgroup/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC."
snps_path="/data/srlab/ssakaue/workspace/scATAC/S-LDSC/snp/weights.hm3_noMHC."

for annot in T-0 T-1 T-2 T-3 T-4 T-5 T-6 T-7 T-8 T-9 T-10 T-11 T-12 T-13 T-14 T-15 T-16 T-17 T-18 T-20 T-21 T-22 T-23
do
    fname="$date${annot}_5kPeaks_binary_hg19."
    annot_path="$basepath/cell_subtypes/$fname"

    for i in {1..22}
    do
        command="python ldsc.py --l2 --bfile $bfile${i} --ld-wind-cm 1 --annot $annot_path${i}.annot.gz --out $annot_path${i} --print-snps $snps_path${i}.snps --thin-annot"
        bsub -J "${annot}LD"[1-10] -q $queue -n $cpus -R 'select[hname!=cn001]' -R 'select[hname!=cn002]' -R 'select[hname!=cn003]' -R 'select[hname!=cn004]' -R 'select[hname!=cn005]' $command
    done
done