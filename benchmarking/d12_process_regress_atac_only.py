from regression_dependencies import *

## READ IN DATA

CT = 'T'
ATAC = sc.read_h5ad("/data/srlab/agupta/data/d12_matrices/processed_h5ads/peaksXcells_"+CT+"_norm_scaled.h5ad")

## ASSIGN METADATA TO CELLS: DONOR ID AND TOTAL COUNTS ACROSS GENES
ATAC = ATAC[RNA_subsetted.obs_names,:]
ATAC.obs['donor'] = [str.split(cname,"_")[0] for cname in RNA_subsetted.obs_names]
ATAC = ATAC[ATAC.obs['ct_subtype']!='T-19: MT-high (low quality)']
ATAC = ATAC[ATAC.obs['ct_subtype']!='M-13: pDC']
ATAC = ATAC[ATAC.obs['ct_subtype']!='B-4: AICDA+BCL6+ GC-like']
ATAC = ATAC[ATAC.obs['ct_subtype']!='B-6: IgM+ plasma']

## HARMONY

## PRE HARMONY
import harmonypy as hm

numPCs = 20
## donor level info
covar_df = pd.DataFrame(ATAC.obs['donor'])

### prep ATAC PC-based input for harmony
print('calculating PCs')
sc.tl.pca(ATAC, n_comps=numPCs)

print('calculating pre-harmony UMAP')
sc.pp.neighbors(ATAC, n_neighbors=100, n_pcs=numPCs,
               use_rep='X_pca')
sc.tl.umap(ATAC)

## RUN HARMONY
ATAC_PCs_input = ATAC.obsm['X_pca'] #use ATAC-based PCs in our 'ATAC only' world
ho = hm.run_harmony(ATAC_PCs_input, covar_df, ['donor'], theta=1)
# Write the Harmony-adjusted PCs to a new file.
res = pd.DataFrame(ho.Z_corr).T

#add donor and total counts as covariates
res['donor']=list(ATAC.obs['donor'])
res.index = ATAC.obs_names
ATAC.obsm['X_harmonized_pca'] = res.iloc[:,:-3]

new_cols=[]
for col in res.columns:
    if type(col)==int:
        new_cols.append("PC"+str(col+1))
    else:
        new_cols.append(str(col))
res.columns=new_cols

#UMAP
print('calculating POST-harmony UMAP')
sc.tl.umap(ATAC)

## prepare for regressions
print('reading in covariates df')
inp_folder = 'srlab'
threshold = 25

# read in peaks X cells DF
ATAC_filepath = "/data/"+inp_folder+"/agupta/data/d12_matrices/h5ads/peaksXcells_"+CT+".h5ad"

print('reading in ATAC cell type subset matrix')
ATAC_CT = read_ATAC(ATAC_filepath)
ATAC_CT_subsetted = ATAC_CT[res.index,:]
ATAC_CT_subsetted.obs['total_ATAC_counts'] = ATAC_CT[res.index,:].X.sum(axis=1)
print('subsetted peaks X cells matrix shape: ', ATAC_CT_subsetted.shape, res.shape)

# subset to peaks open in at least threshold of all cells in cell type
ATAC_BINARIZED = anndata.AnnData(X = (ATAC_CT_subsetted.X>0)*1)
num_cells_open = list(np.array(ATAC_BINARIZED.X.sum(axis=0)))
peaks_open_enough = [idx for idx, val in enumerate(num_cells_open[0]) if val > threshold] ##ARG
print('number of peaks open in at least', str(threshold), 'cells: ', len(peaks_open_enough))

# MAKE A DF OF COVARIATES:
reg_input_df = pd.DataFrame()
reg_input_df['donor'] = res['donor'].astype('category')
reg_input_df[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']] = res[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']] #use ATAC-based PCs in our 'ATAC only' world
reg_input_df['ATAC_total_counts'] = list(ATAC_CT_subsetted.obs['total_ATAC_counts'])
reg_input_df['log_ATAC_total_counts'] = np.log1p(reg_input_df['ATAC_total_counts'])
reg_input_df_copy = reg_input_df.copy()

def get_glm_fit_predict(formula, df, family, **fit_kwargs):
    # works if family in ('binomial', 'poisson')
    f = Formula(formula)
    with localconverter(ro.default_converter + pandas2ri.converter):
        modelfit = stats.glm(f, df, family=family, **fit_kwargs) #either use stats.glm or lme4.glmer
        prediction = stats.predict(modelfit, df, type="response")

    coef_df = r['as.data.frame'](stats.coef(base.summary(modelfit)))
    coef_df = pandas2ri.rpy2py(coef_df)

    return modelfit, prediction, coef_df

def get_LRT(null_model, full_model):
    model_lrt = stats.anova(null_model, full_model, test='LRT')
    model_lrt_df = pandas2ri.rpy2py(r['as.data.frame'](model_lrt))
    chisq_pval = model_lrt_df['Pr(>Chi)'][1]

    return chisq_pval

## REGRESSIONS
peaks_subset_num = 1
peaks = ATAC_CT_subsetted[:,peaks_open_enough].var_names[::peaks_subset_num]
print("num peaks:", len(peaks))

chisq_pvals = {}
regression_pvals, regression_zscores = {}, {}

for peak in peaks:
    reg_input_df_copy['peak_raw_counts'] = ATAC_CT_subsetted[:,peak].X.A

    # to get overall variance explained
    print('fitting & predicting null and full models')
    #null model
    null_formula = "peak_raw_counts ~  offset(log(ATAC_total_counts)) + donor"
    #full model
    full_formula = "peak_raw_counts ~  offset(log(ATAC_total_counts)) + donor + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"

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
    print('done with peak', peak)

overall_fit_df = pd.DataFrame(chisq_pvals,index=['overall_model_fit_pval']).T

# FDR-BH
overall_fit_qvals_df = overall_fit_df.copy()
overall_fit_qvals_df['overall_model_fit_qval'] = BH(overall_fit_qvals_df['overall_model_fit_pval'])

overall_df = overall_fit_qvals_df.copy()