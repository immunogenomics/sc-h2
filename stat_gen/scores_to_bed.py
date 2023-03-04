import pandas as pd
import numpy as np
import scipy as sp
import anndata as adata
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--cell_type', type=str)
parser.add_argument('--date', type=str)
parser.add_argument('--peak_window', type=int)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

## inputting peak-based regression scores
head_folder = "/data/srlab/agupta/data/peak_gene_scores/d12_scores"
fname = args.date+"_poisson_regression_all_outputs.csv"                 
atac_fn = head_folder+'/'+args.cell_type+'/'+fname

peak_scores = pd.DataFrame(pd.read_csv(atac_fn,index_col=0))
peak_scores = peak_scores[peak_scores.index.str.contains('chrX')==False]
peak_scores = peak_scores[peak_scores.index.str.contains('chrY')==False]

peak_scores['DYNAMIC'] = [int(i) for i in peak_scores['model_LRT_qval']<0.05]
peak_scores['INVARIANT'] = 1-peak_scores['DYNAMIC']
dynamic_peaks = list(peak_scores[peak_scores['DYNAMIC']==1].index)
print('num dynamic peaks:', len(dynamic_peaks))
peak_scores['OPEN'] = [1]*len(peak_scores)

#update zcores to be 0 / pvals to be 1 where the peaks are NOT DYNAMIC
peak_scores.loc[peak_scores["model_LRT_qval"] >= 0.05, ['PC1 zscore','PC2 zscore','PC3 zscore','PC4 zscore','PC5 zscore',
                                                       'PC6 zscore','PC7 zscore','PC8 zscore','PC8 zscore','PC10 zscore']] = 0
peak_scores.loc[peak_scores["model_LRT_qval"] >= 0.05, ['PC1 pval','PC2 pval','PC3 pval','PC4 pval','PC5 pval',
                                                       'PC6 pval','PC7 pval','PC8 pval','PC9 pval','PC10 pval']] = 1

## binarize PC-based scores
zscores = peak_scores[['PC1 zscore','PC2 zscore','PC3 zscore','PC4 zscore','PC5 zscore','PC6 zscore','PC7 zscore','PC8 zscore','PC8 zscore','PC10 zscore']]
pvals = peak_scores[['PC1 pval','PC2 pval','PC3 pval','PC4 pval','PC5 pval','PC6 pval','PC7 pval','PC8 pval','PC9 pval','PC10 pval']]
pval_binary = pvals[pvals<0.05].fillna(0)
pval_binary[pval_binary!=0] = 1

#pos
pos_peaks_df = zscores[zscores>0].fillna(0)
pos_peaks_df[pos_peaks_df!=0]=1
pos_peaks_df.columns = ['PC'+str(i)+'_POS' for i in range(1,11)]
# #neg
neg_peaks_df = zscores[zscores<0].fillna(0)
neg_peaks_df[neg_peaks_df!=0]=1
neg_peaks_df.columns = ['PC'+str(i)+'_NEG' for i in range(1,11)]
print('num pos and neg PC peaks', pos_peaks_df.sum(),neg_peaks_df.sum())

# build DF with all annotations that we want to convert to .bed
full_df = pd.concat([peak_scores[['OPEN','DYNAMIC','INVARIANT']],pos_peaks_df,neg_peaks_df],axis=1)


## CONVERT TO INDIVIDUAL .BED FILES

chromosome = [str(i).split(":")[0] for i in full_df.index]
peak_start_bp = [str(i).split(":")[1].split("-")[0] for i in full_df.index]
peak_end_bp = [str(i).split(":")[1].split("-")[1] for i in full_df.index]

for annot_type in full_df.columns:
    
    peaks_for_bed_df = pd.DataFrame([chromosome, peak_start_bp, peak_end_bp, full_df[annot_type]]).T
    peaks_for_bed_df.columns = ['chr','start','end','SCORE']

    # #to make the window around peaks wider
    peaks_for_bed_df_extended = peaks_for_bed_df.copy()
    peaks_for_bed_df_extended['start'] = [int(i) - args.peak_window for i in list(peaks_for_bed_df_extended['start'])]
    peaks_for_bed_df_extended['end'] = [int(i) + args.peak_window for i in list(peaks_for_bed_df_extended['end'])]
    
    # #output without index or header
    bed_folder = '/data/srlab/agupta/data/h2/bed/d12/hg38/' ##output to the hg38 folder

    peaks_for_bed_df_extended.to_csv(bed_folder+args.date+"_"+args.cell_type+"_"+annot_type+"_hg38.bed",
                                     header=False,
                                     index=False,
                                     sep='\t')