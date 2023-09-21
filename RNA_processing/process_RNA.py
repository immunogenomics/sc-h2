from regression_dependencies import *
!conda info --envs

### read in raw data ###

CT = 'T' #'B' 'myeloid' 'fibroblast' 'endothelial'
peaks = sc.read_h5ad("/data/srlab/agupta/data/d12_matrices/h5ads/peaksXcells_"+CT+".h5ad")
genes = sc.read_h5ad("/data/srlab/agupta/data/d12_matrices/h5ads/genesXcells_"+CT+".h5ad")
RNA = genes.copy()
ATAC = peaks.copy()

# assign donor ID to cells
RNA.obs['donor'] = [str.split(cname,"_")[0] for cname in RNA.obs_names]
ATAC.obs['donor'] = [str.split(cname,"_")[0] for cname in RNA.obs_names]

### filtering ###
sc.pl.highest_expr_genes(RNA, n_top=20)
sc.pp.filter_genes(RNA, min_cells=50)
sc.pp.filter_cells(RNA, min_genes=500)

RNA.obs['total_counts'] = RNA.X.sum(axis=1)
ATAC.obs['total_counts'] = ATAC.X.sum(axis=1)

# filter out cells with a high mitochondrial percentage
RNA.var['mt'] = RNA.var_names.str.startswith('MT-')
total_MT_counts = RNA[:,RNA.var['mt']].X.sum(axis=1).A.T[0]
MT_percent = [total_MT_counts[i]/list(RNA.obs['total_counts'])[i] for i in range(len(RNA.obs_names))]
RNA.obs['percent_mito'] = MT_percent
RNA = RNA[RNA.obs['percent_mito']<0.2]

# filter out MT and RP genes
noMTRP = [i for i in RNA.var_names if not (i.startswith("MT-") or i.startswith("RP"))]
RNA = RNA[:,noMTRP]

### RNA-based normalizations and scaling, prior to defining PCs ###

# normalize counts (CP10K)
sc.pp.normalize_total(RNA, target_sum=1e4)

# log transform
sc.pp.log1p(RNA, base=2)

# keep only highly variable genes, for defining PCs
sc.pp.highly_variable_genes(RNA, batch_key='donor', subset=True, n_top_genes=5000)
sc.pl.highly_variable_genes(RNA)

### harmonize across donors ###

import harmonypy as hm

numPCs = 20
covar_df = pd.DataFrame(RNA.obs['donor'])
covar_df['total_counts'] = RNA.obs['total_counts']

sc.tl.pca(RNA, n_comps=numPCs)

donor_dict = dict(zip(list(set(list(RNA.obs['donor']))), np.arange(0,12,1))) 
donor_nums = [donor_dict[d] for d in list(RNA.obs['donor'])]
RNA.obs['donor_num'] = donor_nums

# run Harmony
RNA_PCs_input = RNA.obsm['X_pca']
ho = hm.run_harmony(RNA_PCs_input, covar_df, ['donor'], theta=1)

# Write the Harmony-adjusted PCs to a new file
res = pd.DataFrame(ho.Z_corr).T

#add donor and total counts as covariates
res['donor']=list(RNA.obs['donor'])
res['donor_num'] = donor_nums
res['total_counts']= list(RNA.obs['total_counts'])
res.index = RNA.obs_names
RNA.obsm['X_harmonized_pca'] = res.iloc[:,:-3]

new_cols=[]
for col in res.columns:
    if type(col)==int:
        new_cols.append("PC"+str(col+1))
    else:
        new_cols.append(str(col))
res.columns=new_cols

# post-harmony UMAP checks
sc.pp.neighbors(RNA, n_neighbors=200, n_pcs=numPCs,
               use_rep='X_harmonized_pca')
sc.tl.umap(RNA)

# CELL SUBTYPE #
plt.figure(figsize=(5,6))
sc.pl.umap(RNA, color=['donor'],
           palette='Set3',size=80, frameon=False)
sc.pl.umap(RNA, color=['ct_subtype'],
           palette='tab20',size=80, frameon=False)
sc.pl.umap(RNA, color=['total_counts','percent_mito','n_genes'],
          color_map='cool', use_raw=False, frameon=False, wspace=0.1)

# res has the harmony-adjusted PCs and new cell-hPC values, which we can use in our regressions
