#load in dependencies
import pandas as pd
import numpy as np
import scipy
import random
import scipy.sparse
from scipy.stats import spearmanr, pearsonr
import scanpy as sc
import anndata
import seaborn
from sklearn import decomposition
from sklearn import random_projection
from sklearn.linear_model import LinearRegression
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pylab as pylab
params = {'legend.fontsize': '18',
          'figure.figsize': (10,10),
         'axes.labelsize': '18',
         'axes.titlesize':'22',
         'xtick.labelsize':'14',
         'ytick.labelsize':'14',
         'axes.linewidth': '0.5',
         'pdf.fonttype': '42',
         'font.sans-serif': 'Helvetica'}
pylab.rcParams.update(params)
plt.style.use('seaborn-white')

from intervaltree import Interval, IntervalTree

from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.decomposition import PCA

from statsmodels.genmod.bayes_mixed_glm import PoissonBayesMixedGLM
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import Poisson
from statsmodels.stats.multitest import multipletests

##calculate lambda (to check for pvalue inflation):
def inflation(ps):
    chisq = scipy.stats.chi2.ppf(1-ps, 1)
    lambda_val = np.median(chisq) / scipy.stats.chi2.ppf(0.5, 1)
    print('inflation: ', lambda_val)

def read_ATAC(filepath):
    peaks = sc.read_h5ad(filepath)
    try:
        peaks.obs_names = [str.split(i,"_")[0]+"_"+str.split(i,"#")[1][:-2] for i in peaks.obs_names]
    except:
        print('peak names already adjusted')
        pass
    return peaks

from pathlib import Path
from tqdm.auto import tqdm
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests
from statsmodels.genmod.bayes_mixed_glm import PoissonBayesMixedGLM
from IPython.display import display
from collections import Counter
import random
import time

# R integration
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, numpy2ri, r, Formula
from rpy2.robjects.vectors import StrVector, FloatVector, ListVector
import rpy2.robjects as ro

lme4 = importr('lme4')
lmer = importr('lmerTest') # overloads lmer function from lme4 package
base = importr('base')
stats = importr('stats')
rcompanion = importr('rcompanion')

## METHOD FOR GLMM ##
def fit_lme(formula, df, family='gaussian', random_effect=True, **fit_kwargs):
    f = Formula(formula)

    with localconverter(ro.default_converter + pandas2ri.converter):
        if family == 'gaussian':
            if random_effect:
                control = lme4.lmerControl(**{'calc.derivs': True,
                                              'check.rankX': 'silent.drop.cols',
                                              'check.conv.singular': r('lme4::.makeCC')(action = "ignore",  tol = 1e-4)})
                fit = lmer.lmer(f, df, control=control, **fit_kwargs)
            else:
                fit = stats.lm(f, df, **fit_kwargs)
        elif family in ('binomial', 'poisson'):
            if random_effect:
                fit = lme4.glmer(f, df, nAGQ=0, family=family, **fit_kwargs)  #nAGQ = 0 --> accepts fit when 'good enough'
            else:
                fit = stats.glm(f, df, family=family, **fit_kwargs)
        else:
            if random_effect:
                control = lme4.glmerControl(**{'optimizer': 'nloptwrap',
                                   'calc.derivs': True,
                                   'check.rankX': 'silent.drop.cols',
                                   'check.conv.singular': r('lme4::.makeCC')(action = "ignore",  tol = 1e-4)})
                fit = r('lme4::glmer.nb')(f, df, **{'nb.control': control}, **fit_kwargs)
            else:
                fit = r('MASS::glm.nb')(f, df, **fit_kwargs)

    coef_df = r['as.data.frame'](stats.coef(base.summary(fit)))
    coef_df = pandas2ri.rpy2py(coef_df)

    return fit, coef_df

### FOR NO RANDOM EFFECT (GLM, not GLMM)
def get_glm_fit_predict(formula, df, family, **fit_kwargs):
    # works if family in ('binomial', 'poisson')
    f = Formula(formula)
    with localconverter(ro.default_converter + pandas2ri.converter):
        modelfit = stats.glm(f, df, family=family, **fit_kwargs) #either use stats.glm or lme4.glmer (depending on whether using random effect for donor --> lme4)
        prediction = stats.predict(modelfit, df, type="response")

    coef_df = r['as.data.frame'](stats.coef(base.summary(modelfit)))
    coef_df = pandas2ri.rpy2py(coef_df)

    return modelfit, prediction, coef_df

#use an ANOVA with a chi-sq/likelihood ratio test (LRT) to determine whether overall model is sig
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova.glm
def get_LRT(null_model, full_model):
    model_lrt = stats.anova(null_model, full_model, test='LRT')
    model_lrt_df = pandas2ri.rpy2py(r['as.data.frame'](model_lrt))
    try: #if no random effect
        chisq_pval = model_lrt_df['Pr(>Chi)'][1]
    except: #if random effect
        chisq_pval = model_lrt_df['Pr(>Chisq)'][1]

    return chisq_pval

def BH(pvals):
            sig, pval, _, _ = multipletests(list(pvals), alpha=0.1, method='fdr_bh')
            return pval

inp_folder = 'srlab'
CT='T'
threshold=25

print('reading in covariates df')
res = pd.read_csv("/data/"+inp_folder+"/agupta/data/d12_matrices/cellPC-loadings/"+CT+"_RNA_PC_loadings_noScaling.csv",index_col=0)
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
peaks_open_enough = [idx for idx, val in enumerate(num_cells_open[0]) if val > threshold]
print('number of peaks open in at least', str(threshold), 'cells: ', len(peaks_open_enough))

# MAKE A DF OF COVARIATES:
reg_input_df = pd.DataFrame()
reg_input_df['donor'] = res['donor'].astype('category')
reg_input_df['ATAC_total_counts'] = list(ATAC_CT_subsetted.obs['total_ATAC_counts'])
reg_input_df['log_ATAC_total_counts'] = np.log1p(reg_input_df['ATAC_total_counts'])
reg_input_df[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']] = res[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']]
reg_input_df_copy = reg_input_df.copy()

## REGRESSIONS
peaks_subset_num=3000
peaks = ATAC_CT_subsetted[:,peaks_open_enough].var_names[::peaks_subset_num]
print("num peaks:", len(peaks))

chisq_pvals = {}
regression_pvals, regression_zscores = {}, {}

for peak in peaks:
    print(peak)
    reg_input_df_copy['peak_raw_counts'] = ATAC_CT_subsetted[:,peak].X.A

    ### TO GET OVERALL VARIANCE EXPLAINED ###
    print('fitting & predicting null and full models')
    ###### IF DONOR IS A RANDOM EFFECT #######
    #null model
#     null_formula = "peak_raw_counts ~  offset(log(ATAC_total_counts)) + (1|donor)"
#     #full model
#     full_formula = "peak_raw_counts ~  offset(log(ATAC_total_counts)) + (1|donor) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"

#     #yhat_null
#     nullmodelfit, coef_df = fit_lme(null_formula, reg_input_df_copy, family='poisson', random_effect=True)
#     #yhat_alt
#     ### run PCs-peak multiple regression; coef_df: to get each covariate's info
#     altmodelfit, coef_df = fit_lme(full_formula, reg_input_df_copy, family='poisson', random_effect=True)

    ##### IF DONOR IS A FIXED EFFECT ######
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

    #get PC-peak pvals and zscores
    PC_peak_pvals = list(coef_df['Pr(>|z|)'][-10:])
    regression_pvals[peak] = PC_peak_pvals

    PC_peak_Zscores = list(coef_df['z value'][-10:])
    regression_zscores[peak] = PC_peak_Zscores
    print('done with peak', peak)

zscores_df = pd.DataFrame(regression_zscores,index=['PC'+str(i)+' zscore' for i in range(1,11)]).T
pvals_df = pd.DataFrame(regression_pvals,index=['PC'+str(i)+' pval' for i in range(1,11)]).T
overall_fit_df = pd.DataFrame(chisq_pvals,index=['model_LRT_pval']).T

# FDR-BH
qvals_df = pvals_df.copy()
for PC in ['PC'+str(i)+' pval' for i in range(1,11)]:
    #defined BH correction in dependencies.py file
    qvals_df[PC] = BH(qvals_df[PC])
qvals_df.columns=['PC'+str(i)+' qval' for i in range(1,11)]

overall_fit_qvals_df = overall_fit_df.copy()
overall_fit_qvals_df['model_LRT_qval'] = BH(overall_fit_qvals_df['model_LRT_pval'])

overall_df = pd.concat([overall_fit_qvals_df,zscores_df,pvals_df,qvals_df],axis=1)



### READ IN PVALS FROM BOTH KINDS OF MODELS:
## 1. DONOR AS FIXED EFFECT (original)
## 2. DONOR AS RANDOM EFFECT (modified)

CT="myeloid"
random_df = pd.read_csv('/data/srlab/agupta/data/peak_gene_scores/test_fixed_v_random/'+CT+'_poisson_donorAsRandom_outputs.csv',index_col=0)
fixed_df = pd.read_csv('/data/srlab/agupta/data/peak_gene_scores/d12_scores/'+CT+'/083122_poisson_regression_10PC_all_outputs.csv',index_col=0)
random_df, fixed_df

plt.figure(figsize=(5,5))
plt.scatter(-np.log10(fixed_df['model_LRT_pval']),
                     -np.log10(random_df['model_LRT_pval']),color='grey')
print(spearmanr(-np.log10(fixed_df['model_LRT_pval']),
                     -np.log10(random_df['model_LRT_pval'])))
top_val = np.max(-np.log10(fixed_df['model_LRT_pval']))
plt.plot([0,205],[0,205],color='black')
plt.xlabel("Donor as Fixed Effect\n-log10(p)")
plt.ylabel("Donor as Random Effect\n-log10(p)")
sns.despine()

