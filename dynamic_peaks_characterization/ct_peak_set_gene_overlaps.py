from regression_dependencies import *

## identifying which genes' promoters different peak sets overlap

CT = 'endothelial'

head_folder = "/data/srlab/agupta/data/peak_gene_scores/d12_scores"
fname = "poisson_regression_10PC_all_outputs.csv"

atac_fn = head_folder+'/'+CT+'/'+fname
peak_scores = pd.DataFrame(pd.read_csv(atac_fn,index_col=0))
peak_scores = peak_scores[peak_scores.index.str.contains('chrX')==False]
peak_scores = peak_scores[peak_scores.index.str.contains('chrY')==False]

peak_scores['DYNAMIC'] = [int(i) for i in peak_scores['model_LRT_qval']<0.05]
peak_scores['INVARIANT'] = 1-peak_scores['DYNAMIC']
dynamic_peaks = list(peak_scores[peak_scores['DYNAMIC']==1].index)

dynamic_peaks = peak_scores[peak_scores['DYNAMIC']==1]
invariant_peaks = peak_scores[peak_scores['INVARIANT']==1]
print('num dynamic peaks:', len(dynamic_peaks))
print('num invariant peaks:', len(invariant_peaks))
peak_scores['OPEN'] = [1]*len(peak_scores)

DYN_peak_scores = peak_scores.loc[:, peak_scores.columns != 'INVARIANT']

#only keep those peaks whose individual pvalues are <0.05:
sig_zscores = zscores[pval_binary==1].fillna(0)

DYN_peak_scores['INVARIANT'] = 1-DYN_peak_scores['DYNAMIC']
full_df = DYN_peak_scores[['DYNAMIC','INVARIANT']]

"""### biological functional enrichment of each cluster of dynamic peaks (across the cell types)"""

RNA_all = sc.read_h5ad("/data/srlab/agupta/data/d12_matrices/h5ads/genesXcells_"+CT+".h5ad")
sc.pp.filter_genes(RNA_all, min_cells=50)

#get genes with TSS info
gene_TSS = pd.DataFrame(pd.read_csv("/data/srlab/agupta/data/ref_files/genes_with_TSS.txt",index_col=3,sep='\t'))
# expand the window around TSS by some bps to get roughly 'promoter' regions
window_size = 1e3
window_starts, window_ends = [], []
for gene in gene_TSS.index:
    # if gene_TSS.loc[gene,'strand'] == '+': ## DECIDE WHETHER TO SWITCH ORIENTATION
    window_starts.append(max(0, int(gene_TSS.loc[gene,'TSS']-window_size)))
    window_ends.append(max(0,int(gene_TSS.loc[gene,'TSS']+window_size)))

gene_TSS['promoter_window_starts'], gene_TSS['promoter_window_ends'] = window_starts, window_ends
gene_TSS_for_overlaps = gene_TSS[['chr','promoter_window_starts','promoter_window_ends']]

## use all genes but only those expressed in at least N (ie 50) cells (across the 5 cell types)
genes_keep = list(set(RNA_all.var_names).intersection(set(gene_TSS_for_overlaps.index)))
genes_keep_adata = RNA_all[:,genes_keep]
genes_keep_adata = genes_keep_adata[:,genes_keep_adata.var['n_cells']>50]
final_genes_check = list(genes_keep_adata.var_names)

gene_TSS_for_overlaps_final = gene_TSS_for_overlaps.loc[final_genes_check,:]
gene_TSS_for_overlaps_final

## mapping genes to peak set of interest

peak_gene_mapping_dict = dict()

for ptype in ['DYNAMIC','INVARIANT']:
    print(ptype)
    this_peak_set = full_df[full_df[ptype]==1]
    cluster = CT+'_'+ptype
    ###########################################################

    # peaks within this cluster, structured properly
    cluster_peaks_df = pd.DataFrame(columns=['chr','start','end'])
    cluster_peaks_df['chr'] = [str(i).split(':')[0] for i in list(this_peak_set.index)]
    cluster_peaks_df['start'] = [int(str(i).split(':')[1].split('-')[0])-500 for i in list(this_peak_set.index)]
    cluster_peaks_df['end'] = [int(str(i).split(':')[1].split('-')[1])+500 for i in list(this_peak_set.index)]

    #for each cluster, check overlap of peaks with genes (ie get the set of genes whose promoters overlap peaks)
    peak_set = cluster_peaks_df

    overlapping_genes_set = set()

    for chrnum in np.arange(1,23,1):
        chrname = 'chr'+str(chrnum)
        print(chrname)
        chr_DYN_peaks_df = peak_set[peak_set['chr']==chrname]
        chr_TSS_df = gene_TSS_for_overlaps_final[gene_TSS_for_overlaps_final['chr']==chrname]

        #loop through peak windows (for dynamic peaks for now) -- to build up an interval tree with all peaks
        print('looping through peaks in '+cluster)
        peaks_tree = IntervalTree()
        for peak_i in range(len(chr_DYN_peaks_df)):
            peakstart,peakend = chr_DYN_peaks_df.iloc[peak_i,:]['start'],chr_DYN_peaks_df.iloc[peak_i,:]['end']
            peaks_tree[peakstart:peakend]=chrname+":"+str(peakstart)+"-"+str(peakend)

        #for ACTIVE promoters (expressed genes)
        ##loop through gene promoter regions (TSS +/- 1kb --> 2kb around TSS)
        for gene_prom_i in range(len(chr_TSS_df)):
            pstart,pend=chr_TSS_df.iloc[gene_prom_i,:]['promoter_window_starts'], chr_TSS_df.iloc[gene_prom_i,:]['promoter_window_ends']
            #get all peaks that fall within this gene's TSS window
            if len(sorted(peaks_tree[pstart:pend]))>0:
                overlapping_genes_set.add(chr_TSS_df.index[gene_prom_i])

    peak_gene_mapping_dict[cluster] = overlapping_genes_set

genes_forPkSet = pd.DataFrame()
genes_forPkSet = pd.concat([genes_forPkSet,
                           pd.DataFrame([peak_gene_mapping_dict[CT+"_DYNAMIC"],peak_gene_mapping_dict[CT+"_INVARIANT"]])])
genes_forPkSet = genes_forPkSet.T
genes_forPkSet.columns = ([CT+'_DYNAMIC',CT+'_INVARIANT'])
genes_forPkSet