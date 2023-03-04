from regression_dependencies import *

CT = 'T'

peaks = sc.read_h5ad("/data/srlab/agupta/data/d12_matrices/h5ads/peaksXcells_"+CT+".h5ad")
peaks = peaks[:,full_df.index]
#get total counts for each peak across cells
peaks.var['peak_total_counts'] = peaks.X.sum(axis=0).T
peaks.var['peak_mean_counts'] = peaks.X.mean(axis=0).T


peaks_df = peaks[:,list(invariant_peaks.index)+list(dynamic_peaks.index)].to_df().T
peaks_df['peak_type'] = ['invariant']*len(invariant_peaks.index) + ['dynamic']*len(dynamic_peaks.index)
peaks_df['log1p_mean_counts'] = [np.log1p(i) for i in list(peaks.var['peak_mean_counts'])]
peaks_df['log1p_total_counts'] = [np.log1p(i) for i in list(peaks.var['peak_total_counts'])]
peaks_df['peak_mean_counts'] = list(peaks.var['peak_mean_counts'])

## TEST FOR INDEPENDENCE OF THE TWO MEANS (INVARIANT VS DYNAMIC PEAK COUNTS)

from scipy.stats import ttest_ind, wilcoxon, ranksums

invariant_peak_counts = peaks_df[peaks_df['peak_type']=='invariant']['peak_mean_counts']
dynamic_peak_counts = peaks_df[peaks_df['peak_type']=='dynamic']['peak_mean_counts']

print('ttest:', ttest_ind(dynamic_peak_counts, invariant_peak_counts, alternative='two-sided'))
#whereas wilcoxon requires same sample sizes, ranksum allows for unequal sample sizes
print('wilcoxon rank sum:', ranksums(dynamic_peak_counts, invariant_peak_counts, alternative='two-sided'))

plt.figure(figsize=(5,5))
sns.violinplot(data=peaks_df,x='peak_type',y='peak_mean_counts',palette=['black','grey'])
sns.despine()
plt.show()

# Number of 1) accessible and 2) dynamic peaks in each cell type

open_peaks_ct, dyn_peaks_ct = [], []
for CT in ['B','T','myeloid','fibroblast','endothelial']:
    ct_reg_df = pd.DataFrame(pd.read_csv('/data/srlab/agupta/data/peak_gene_scores/d12_scores/'+CT+'/083122_poisson_regression_10PC_all_outputs.csv',index_col=0))
    print(CT)
    print('total open peaks in ',CT, ':', len(ct_reg_df))
    print('num dynamic peaks in ',CT,':', len(ct_reg_df[ct_reg_df['model_LRT_qval']<0.05]))
    print(len(ct_reg_df[ct_reg_df['model_LRT_qval']<0.05]) / len(ct_reg_df))
    open_peaks_ct.append(len(ct_reg_df))
    dyn_peaks_ct.append(len(ct_reg_df[ct_reg_df['model_LRT_qval']<0.05]))
    
plt.figure(figsize=(4,3))
plt.bar(np.arange(1,6,1), open_peaks_ct,color='black',label='open')
plt.bar(np.arange(1,6,1), dyn_peaks_ct,color='grey',label='dynamic')
plt.ylabel("Number of peaks")
# plt.legend()
# plt.xticks(['B','T','myeloid','fibroblast','endothelial'])
sns.despine()