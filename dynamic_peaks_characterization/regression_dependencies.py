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
    return lambda_val

def BH(pvals):
            sig, pval, _, _ = multipletests(list(pvals), alpha=0.1, method='fdr_bh')
            return pval

## FOR GENERALIZED LINEAR MIXED EFFECT MODEL-BASED REGRESSIONS ##

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
                fit = stats.glm(f, df, family=family, **fit_kwargs) #nAGQ=0, 
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
    
    return coef_df

### FOR NO RANDOM EFFECT (GLM, not GLMM)
def get_glm_fit_predict(formula, df, family, **fit_kwargs):
    # works if family in ('binomial', 'poisson')
    f = Formula(formula)
    with localconverter(ro.default_converter + pandas2ri.converter):
        modelfit = stats.glm(f, df, family=family, **fit_kwargs)
        prediction = stats.predict(modelfit, df, type="response")
    
    coef_df = r['as.data.frame'](stats.coef(base.summary(modelfit)))
    coef_df = pandas2ri.rpy2py(coef_df)
        
    return modelfit, prediction, coef_df

#use an ANOVA with a chi-sq/likelihood ratio test (LRT) to determine whether overall model is sig
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova.glm
def get_LRT(null_model, full_model):
    model_lrt = stats.anova(null_model, full_model, test='LRT')
    model_lrt_df = pandas2ri.rpy2py(r['as.data.frame'](model_lrt))
    chisq_pval = model_lrt_df['Pr(>Chi)'][1]
    
    return chisq_pval


def read_ATAC(filepath):
    peaks = sc.read_h5ad(filepath)
    try:
        peaks.obs_names = [str.split(i,"_")[0]+"_"+str.split(i,"#")[1][:-2] for i in peaks.obs_names]
    except:
        print('peak names already adjusted')
        pass
    return peaks

def plot_peak_pval_distribs(peak_PC_df_binarized, filepath, cell_type_res, family, date):
    #how many peaks is each PC a significant predictor for?
    plt.figure(figsize=(8,6))
    plt.bar(np.arange(1,11),
            list(peak_PC_df_binarized.sum(axis=0)),
           color='darkslateblue',alpha=0.7)
    plt.ylabel("number of peaks\nwith a signif. PC beta",fontsize=20)
    plt.title(cell_type_res,fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel("PC",fontsize=20)
    sns.despine()
    plt.savefig(filepath+date+'_'+cell_type_res+'_peak_PCs_bar'+'_'+family) ## UPDATE

    #how many peaks have significant betas for how many PCs?
    plt.figure(figsize=(8,6))
    plt.hist(np.sum(peak_PC_df_binarized,axis=1),
            color='darkslateblue',alpha=0.7,histtype='stepfilled')
    plt.xlabel("# of PCs w significant betas",fontsize=18)
    plt.title(cell_type_res,fontsize=20)
    plt.ylabel("# of peaks",fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    sns.despine()
    plt.savefig(filepath+date+'_'+cell_type_res+'_peak_PCs_hist'+'_'+family) ## UPDATE

def plot_QQ(this_pvals_df, scramble, cell_type_res, family, date):
    uniform_pvals = np.random.uniform(low=0, high=1,size=len(this_pvals_df))
    first_plot, second_plot = uniform_pvals, this_pvals_df
    xmax = np.max(-np.log10(first_plot))

    plt.style.use('default')
    fig, axs = plt.subplots(nrows=2, ncols=5, sharex=True, sharey=True,figsize=(12,8))

    # subset_to = -1
    X = sorted(-np.log10(first_plot))

    ## PC1
    axs[0, 0].scatter(X,sorted(-np.log10(second_plot['PC1 pval'])),s=5,color='black')
    axs[0, 0].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC1 pval']))))
    axs[0, 0].set_title('PC1\ninflation = '+str(np.round(PC_inflation,2)))

    ## PC2
    axs[0, 1].scatter(X,sorted(-np.log10(second_plot['PC2 pval'])),s=5,color='black')
    axs[0, 1].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC2 pval']))))
    axs[0, 1].set_title('PC2\ninflation = '+str(np.round(PC_inflation,2)))

    ## PC3
    axs[0, 2].scatter(X,sorted(-np.log10(second_plot['PC3 pval'])),s=5,color='black')
    axs[0, 2].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC3 pval']))))
    axs[0, 2].set_title('PC3\ninflation = '+str(np.round(PC_inflation,2)))

    ## PC4
    axs[0, 3].scatter(X,sorted(-np.log10(second_plot['PC4 pval'])),s=5,color='black')
    axs[0, 3].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC4 pval']))))
    axs[0,3 ].set_title('PC4\ninflation = '+str(np.round(PC_inflation,2)))

    ## PC5
    axs[0, 4].scatter(X,sorted(-np.log10(second_plot['PC5 pval'])),s=5,color='black')
    axs[0, 4].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC5 pval']))))
    axs[0, 4].set_title('PC5\ninflation = '+str(np.round(PC_inflation,2)))
    
    ## PC6
    axs[1, 0].scatter(X,sorted(-np.log10(second_plot['PC6 pval'])),s=5,color='black')
    axs[1, 0].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC6 pval']))))
    axs[1, 0].set_title('PC6\ninflation = '+str(np.round(PC_inflation,2)))

    ## PC7
    axs[1, 1].scatter(X,sorted(-np.log10(second_plot['PC7 pval'])),s=5,color='black')
    axs[1, 1].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC7 pval']))))
    axs[1, 1].set_title('PC7\ninflation = '+str(np.round(PC_inflation,2)))

    ## PC8
    axs[1, 2].scatter(X,sorted(-np.log10(second_plot['PC8 pval'])),s=5,color='black')
    axs[1, 2].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC8 pval']))))
    axs[1, 2].set_title('PC8\ninflation = '+str(np.round(PC_inflation,2)))

    ## PC9
    axs[1, 3].scatter(X,sorted(-np.log10(second_plot['PC9 pval'])),s=5,color='black')
    axs[1, 3].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC9 pval']))))
    axs[1, 3].set_title('PC9\ninflation = '+str(np.round(PC_inflation,2)))

    ## PC10
    axs[1, 4].scatter(X,sorted(-np.log10(second_plot['PC10 pval'])),s=5,color='black')
    axs[1, 4].plot([0,xmax],[0,xmax], color='grey', linestyle='--')
    PC_inflation = inflation(np.array(sorted(np.array(second_plot['PC10 pval']))))
    axs[1, 4].set_title('PC10\ninflation = '+str(np.round(PC_inflation,2)))

    if scramble==True:
        status='permuted'
    else:
        status='raw'
        
    plt.ylabel("-log10P(Observed)\n"+status)
    plt.xlabel("-log10P(Expected)\nUniform Distribution")
    sns.despine()
    plt.savefig('/data/srlab/agupta/data/d12_regression_plots/'+date+'_QQ_'+cell_type_res+'_'+status+'_'+family)