import pandas as pd
import numpy as np
import os
from intervaltree import Interval, IntervalTree
import argparse
import sys
import glob
import matplotlib.pyplot as plt

import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--annot', type=str)
parser.add_argument('--cell_type', type=str)
parser.add_argument('--date', type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

#make annot files from .bed
#bed files should include chr, start, stop, weight for their columns

def make_annots(annots_dest, name, bed_folder, basedir):
    common_SNPS_scores = []
    
    for i in range(1, 23): #loop through chromosomes
        os.chdir(annots_dest)
        print(i)
        if (os.path.isfile("%s.%s.annot.gz" %(name,i))==False):
            os.chdir(basedir+'bed/') #bed folder
            my_df = pd.read_csv(bed_folder+"%s.bed" %name, sep='\t', names=['CHR','START','STOP','WEIGHT'])
            my_df = my_df[my_df['CHR'] == 'chr%s' %i]
            my_df = my_df[['CHR','START','STOP','WEIGHT']]

            #get bim file
            bfile = ("/data/srlab/external-data/LDSCORE/data.broadinstitute.org/alkesgroup/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.%s" %i)
            df_bim = pd.read_csv(bfile + '.bim',delim_whitespace=True, usecols = [0,1,2,3],
                                     names = ['CHR','SNP','CM','BP'])
            sbim = pd.Series(df_bim['BP'])

            #get frequency file (to determine whether SNPs are common, MAF>5%, or not)
            print('getting frq file')
            frqfile = ("/data/srlab/external-data/LDSCORE/data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq/1000G.EUR.QC.%s" %i)
            df_frq = pd.read_csv(frqfile + '.frq',delim_whitespace=True, usecols = [0,1,4]) #cols: CHR, SNP, MAF
            df_frq['BP'] = df_bim['BP'] #to map rsID to BP positions (that's what our scores are mapped to)
            #to keep common SNPs
            common_SNPS_in_chr = set(df_frq[df_frq['MAF']>0.05]['BP'])

            print('checking overlap between peaks and SNPs')
            all_snps_present = []
            all_weights = []
            for row in range(len(my_df)):
                this_genes_weight = my_df.iloc[row,:]['WEIGHT']
                start, stop = my_df.sort_values("START").iloc[row,:]['START'],my_df.sort_values("START").iloc[row,:]['STOP']
                snps_in_range = list(sbim[sbim.between(start,stop)])
                if len(snps_in_range)>0:
                    all_snps_present.extend([int(i) for i in snps_in_range])
                    all_weights.extend([this_genes_weight]*len(snps_in_range))

            # to get all SNPs with annotation 1 based on which peaks they are in
            annot_df = pd.DataFrame(np.array([all_snps_present,all_weights]).T,columns=['BP','ANNOT'])

            # for any SNP that isn't within an open/annotated peak, assign an annotation of 0:
            print("merging BP list with scores")
            merged_df = pd.merge(df_bim,annot_df,how='left',on='BP').fillna(0)

            # TAKING THE MAX TO ASSIGN A BINARY SCORE OF 1 TO ANY REGION THAT IS MAPPED AS 1 THROUGH ANY PEAK
            print('taking max annot for a given SNP, ie if peaks overlap')
            merged_df_averaged = merged_df.groupby("BP").max()#mean()
            merged_df_averaged['BP'] = merged_df_averaged.index

            annotfile_final = ("%s.%s.annot.gz" %(name,i))
            os.chdir(annots_dest)
            print(len(merged_df_averaged), merged_df_averaged.head())
            merged_df_averaged = merged_df_averaged.drop(["CHR","BP","CM","SNP"], axis = 1)
            merged_df_averaged.to_csv(annotfile_final,sep='\t',compression='gzip',index=False)

####################################################################################################   
#run make_annots

# bed input file
bed_folder = '/data/srlab/agupta/data/h2/bed/d12/hg19/'

# name of input and output identifier
name = args.date+'_'+args.annot+'_hg19'

# annot output folder
basedir = '/data/srlab/agupta/data/h2/'
annots_dest = basedir+"annotations/d12/cell_score_bins/"

make_annots(annots_dest, name, bed_folder, basedir)