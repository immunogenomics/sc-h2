from regression_dependencies import *

### Reading in and processing (public) raw PBMC data as inputs into regressions
# only one donor, so no donor-based adjustments or covariates

CT = 'T'
genes = sc.read_h5ad('/data/srlab/agupta/data/HB_PBMCs/'+CT+'_RNA.h5ad'
RNA = genes.copy()

## filtering and normalization

sc.pl.highest_expr_genes(RNA, n_top=20)
sc.pp.filter_genes(RNA, min_cells=50)
sc.pp.filter_cells(RNA, min_genes=500)
print("RNA ADATA shape post filtering", RNA.shape)

RNA.obs['total_counts'] = RNA.X.sum(axis=1)

# filter out cells with a high MT percentage
RNA.var['mt'] = RNA.var_names.str.startswith('MT-')
total_MT_counts = RNA[:,RNA.var['mt']].X.sum(axis=1).A.T[0]
print("filtering out cells with >20% of counts from MT")
MT_percent = [total_MT_counts[i]/list(RNA.obs['total_counts'])[i] for i in range(len(RNA.obs_names))]
RNA.obs['percent_mito'] = MT_percent
RNA = RNA[RNA.obs['percent_mito']<0.2]
print(len(total_MT_counts), len(RNA.obs_names), len(RNA.obs['total_counts']))

# filter out MT and RP genes
noMTRP = [i for i in RNA.var_names if not (i.startswith("MT-") or i.startswith("RP"))]
RNA = RNA[:,noMTRP]
print(RNA.shape)

# normalize counts (CP10K)
sc.pp.normalize_total(RNA, target_sum=1e4, inplace=False, exclude_highly_expressed=True)

# log transform
sc.pp.log1p(RNA, base=2)

# keep only highly variable genes, for defining PCs
sc.pp.highly_variable_genes(RNA, subset=True, n_top_genes=5000) ### CHECK THIS
sc.pl.highly_variable_genes(RNA)
print(len(RNA.var_names))
sc.pl.highest_expr_genes(RNA, n_top=20)

## PRE HARMONY
import harmonypy as hm

numPCs = 20
# prep RNA PC-based input for harmony
print('calculating top',numPCs,'PCs')
sc.tl.pca(RNA, n_comps=numPCs)

# UMAP
print('calculating UMAP')
sc.pp.neighbors(RNA, n_neighbors=100, n_pcs=numPCs,
               use_rep='X_pca')
sc.tl.umap(RNA)
RNA.uns['ct_subtype_colors'] = subtypecolors

# inputting Azimuth cell state labels

RNA_for_AZ = sc.read_h5ad('/data/srlab/agupta/data/HB_PBMCs/'+CT+'_MY_processed_RNA.h5ad') #_processed
AZ_cs = pd.read_csv('/data/srlab/agupta/data/HB_PBMCs/'+CT+'_PBMC_Azimuth_cellStates.tsv',sep='\t',index_col=0)
for cs in list(set(AZ_cs['predicted.celltype.l2'])):
    print(cs, len(AZ_cs[AZ_cs['predicted.celltype.l2']==cs]))

if CT=='T':
    AZ_ctcolordict = {'CD4 Naive':'#61B435',
                      'CD4 TCM':'#B983FF',
                      'CD4 TEM':'#A4A608',
                      'CD8 Naive':'#00BCD8',
                      'CD8 TCM':'#B79F00',
                      'CD8 TEM':'#00AFF6',
                      'MAIT':'#0ABF64',
                      'Treg':'#EE7F49',
                      'dnT':'#00C097',
                      'gdT':'#9590FF'}
if CT=='myeloid':
    AZ_ctcolordict = {'cDC2':'#F68B86',
                      'pDC':'#00C0AF',
                      'HSPC':'#FF67A4',
                      'CD14 Mono':'#FD6EC0',
                      'CD16 Mono':'#D376FF'}

if CT=='B':
        AZ_ctcolordict = {'B naive':'#FE6F94',
                      'Plasmablast':'#07BDD9',
                      'B intermediate':'#00B7E9',
                      'B memory':'#00BA38'}

AZ_subtype_colors = list(AZ_ctcolordict.values())

RNA_for_AZ = RNA_for_AZ[AZ_cs.index,:]
## map Azimuth cell state labels to our PBMC cells
RNA_for_AZ.obs['AZ_cs'] = AZ_cs['predicted.celltype.l2']
## subset to only include cell states that have at least 10 cells:
RNA_for_AZ = RNA_for_AZ[RNA_for_AZ.obs['AZ_cs'].isin(AZ_ctcolordict.keys())]
RNA_for_AZ

#reorder cell state-color mappings
if CT=="B":
    AZ_subtype_colors = ['#00B7E9', '#00BA38','#FE6F94', '#07BDD9']
if CT=='myeloid':
    AZ_subtype_colors = ['#FD6EC0', '#D376FF', '#FF67A4','#F68B86', '#00C0AF']

# plot PBMCs for this cell type
sc.pl.umap(RNA_for_AZ, color='AZ_cs', palette=AZ_subtype_colors, size=200,
           frameon=False)

plt.figure(figsize=(8,8))
for st in set(RNA_for_AZ.obs['AZ_cs']):
    color=AZ_ctcolordict[st]
    centroid_x = np.median(list(RNA_for_AZ[RNA_for_AZ.obs['AZ_cs']==st].obsm['X_umap'][:,0]))
    centroid_y = np.median(list(RNA_for_AZ[RNA_for_AZ.obs['AZ_cs']==st].obsm['X_umap'][:,1]))
    plt.scatter([centroid_x], [centroid_y],c=color,s=800, linewidths=4, edgecolors='black')
plt.axis('off')
sns.despine()

# prep res file for output (what will go in as our input covariates df for regressions:
meta = pd.read_csv("/data/srlab/agupta/data/HB_PBMCs/buildingblocks/PBMC_meta_"+CT+".csv")\
res = pd.DataFrame(RNA.obsm['X_pca'][:,:10],index=RNA.obs_names)
res.columns=['PC'+str(i) for i in np.arange(1,11,1)]

## testing regressions code

inp_folder = 'srlab'
threshold=50
peaks_subset_num=10000

print('reading in covariates df')
res = pd.read_csv("/data/"+inp_folder+"/agupta/data/HB_PBMCs/PCloadings/"+CT+"_covariates.csv",index_col=0) ##ARG
# read in peaks X cells DF
ATAC_filepath = "/data/"+inp_folder+"/agupta/data/HB_PBMCs/"+CT+"_ATAC.h5ad"
print('reading in ATAC cell type subset matrix')
ATAC_CT = read_ATAC(ATAC_filepath)
#ENSURE SAME CELLS, IN SAME ORDER
ATAC_CT_subsetted = ATAC_CT[res.index,:]
ATAC_CT_subsetted.obs['ATAC_total_counts'] = ATAC_CT[res.index,:].X.sum(axis=1)
print('subsetted peaks X cells matrix shape: ', ATAC_CT_subsetted.shape, res.shape)

## SUBSET to peaks open in at least threshold of all cells in cell type
ATAC_BINARIZED = anndata.AnnData(X = (ATAC_CT_subsetted.X>0)*1)
num_cells_open = list(np.array(ATAC_BINARIZED.X.sum(axis=0)))
peaks_open_enough = [idx for idx, val in enumerate(num_cells_open[0]) if val > threshold]
print('number of peaks open in at least', str(threshold), 'cells: ', len(peaks_open_enough))

# making a df of covariates
reg_input_df_copy = res.copy()

## REGRESSIONS ##
peaks_subset_num = peaks_subset_num
peaks = ATAC_CT_subsetted[:,peaks_open_enough].var_names[::peaks_subset_num]
print("num peaks we are about to run regressions on:", len(peaks))

chisq_pvals = {}
regression_pvals, regression_zscores = {}, {}

for peak in peaks:
    reg_input_df_copy['peak_raw_counts'] = ATAC_CT_subsetted[:,peak].X.A
    reg_input_df_copy['ATAC_total_counts'] = list(ATAC_CT_subsetted.obs['ATAC_total_counts'])

    print('fitting & predicting null and full models')
    #null model
    null_formula = "peak_raw_counts ~  offset(log(ATAC_total_counts))"
    #full model
    full_formula = "peak_raw_counts ~  offset(log(ATAC_total_counts)) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"

    #yhat_null
    nullmodelfit, yhat_null, _ = get_glm_fit_predict(null_formula, reg_input_df_copy, 'poisson')
    #yhat_alt
    ### run PCs-peak multiple regression; coef_df: to get each covariate's info
    altmodelfit, yhat_alt, coef_df = get_glm_fit_predict(full_formula, reg_input_df_copy, 'poisson')

    #compare null from full model
    print('comparing null vs full models for model fit')
    # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova.glm
    chisq_pval = get_LRT(nullmodelfit, altmodelfit)
    chisq_pvals[peak] = chisq_pval

    #get PC-peak pvals and zscores
    PC_peak_pvals = list(coef_df['Pr(>|z|)'][-10:])
    regression_pvals[peak] = PC_peak_pvals

    PC_peak_Zscores = list(coef_df['z value'][-10:])
    regression_zscores[peak] = PC_peak_Zscores
    print('done with peak', peak)

"""regression outputs to bed input dfs

### regression outputs to bed input dfs
all_df = pd.DataFrame()
for CT in ['T','B','myeloid']:
    print(CT)
    head_folder = "/data/srlab/agupta/data/HB_PBMCs/peak_gene_scores/"+CT+"_PBMC_regressions.csv"
    OG_df = pd.DataFrame(pd.read_csv(head_folder,index_col=0))
    #dynamic peaks for this CT
    print('dyn')
    dyn_df = OG_df[OG_df['model_LRT_FDR']<0.05]
    dyn_df = dyn_df[dyn_df.index.str.contains("chrX") == False]
    dyn_df[CT+"_D"] = [1]*len(dyn_df)
    dyn_df = dyn_df[[CT+"_D"]]
    #invariant peaks for this CT
    print('inv')
    inv_df = OG_df[OG_df['model_LRT_FDR']>=0.05]
    inv_df = inv_df[inv_df.index.str.contains("chrX") == False]
    inv_df[CT+"_I"] = [1]*len(inv_df)
    inv_df = inv_df[[CT+"_I"]]

    all_df = pd.concat([all_df,dyn_df,inv_df],axis=1)
    all_df = all_df.fillna(0)



### Pseudobulk ATAC only --> build cell type annotations

## read in all cell type (B,T,M) ATAC normalized counts matrix
full_ATAC = sc.read_h5ad("/data/srlab/agupta/data/HB_PBMCs/PBMCs_allCT_ATAC.h5ad")
print(full_ATAC)
# subset to cells in the individiual CT matrices & add subtype info
cells_keep = []
subtypes = []
cell_types = []
for CT in ['B','T','myeloid']:
    subATAC = sc.read_h5ad("/data/srlab/agupta/data/HB_PBMCs/"+CT+"_RNA.h5ad")
    cells_keep.extend(list(subATAC.obs_names))
    subtypes.extend(list(subATAC.obs['ct_subtype']))
    cell_types.extend([CT]*len(subATAC))
print(len(cells_keep))
full_ATAC = full_ATAC[cells_keep]
full_ATAC.obs['ct_subtype'] = subtypes
full_ATAC.obs['cell_type'] = cell_types

CT_ATAC_df = full_ATAC.to_df()
CT_ATAC_df['cell_type'] = full_ATAC.obs['cell_type']
CT_ATAC_df = CT_ATAC_df.groupby('cell_type').sum()
CT_ATAC_df['cell_type'] = CT_ATAC_df.index

CT_onehotdf = pd.get_dummies(CT_ATAC_df['cell_type'])
CT_ATAC_df = CT_ATAC_df.drop('cell_type',axis=1)

import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS

#for each cell type:
CT_peak_Ttest = pd.DataFrame()
#get cell type binary vectors from the one hot encoded DF
for CT in CT_onehotdf.columns:
    print(CT)
    X_use = CT_onehotdf[[CT]]
    peak_tvals = []
    #for each peak
    #get the normalized peak accessibility across samples from the CT_ATAC_df
    for peak in CT_ATAC_df.columns:
        Y_use = CT_ATAC_df[[peak]]

        results = OLS(Y_use, X_use).fit()
        r = np.zeros_like(results.params)
        if results.tvalues[0]>15:
            print(peak, np.round(results.tvalues[0],2))
        peak_tvals.append(results.tvalues[0])
    CT_peak_Ttest[CT] = peak_tvals

CT_peak_Ttest.index = CT_ATAC_df.columns

CT_peak_Ttest.replace([np.inf, -np.inf], np.nan, inplace=True)
df_updated = CT_peak_Ttest.dropna(axis=0)
df_updated = df_updated[df_updated.index.str.contains("chrX")==False]

#get top 10% of genes WRT t statistic values for each cell type
top_peaks_perCT = pd.DataFrame()
for CT in df_updated.columns:
    this_CT_peaks = df_updated.nlargest(int(len(df_updated)/10), CT)[[CT]]
    top_peaks_perCT = pd.concat([top_peaks_perCT, this_CT_peaks],axis=1)

top_peaks_perCT = top_peaks_perCT.fillna(0)

#to binarize the scores for LDSC annotations input
top_peaks_perCT[top_peaks_perCT>0]=1
print(np.sum(top_peaks_perCT,axis=0))
top_peaks_perCT