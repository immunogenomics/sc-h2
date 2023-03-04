from regression_dependencies import *
import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--CT', type=str)
parser.add_argument('--threshold', type=int, default=25)
parser.add_argument('--peaks_subset_num', type=int, default=10000)
parser.add_argument('--family', type=str)
parser.add_argument('--i', type=int)
parser.add_argument('--date', type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

inp_folder = 'srlab'

print('reading in covariates df')
res = pd.read_csv("/data/"+inp_folder+"/agupta/data/d12_matrices/cellPC-loadings/"+args.CT+"_RNA_PC_loadings_noScaling.csv",index_col=0) ##ARG
# read in peaks X cells DF
ATAC_filepath = "/data/"+inp_folder+"/agupta/data/d12_matrices/h5ads/peaksXcells_"+args.CT+".h5ad"
## UPDATE: OUTPUT ^ FOR OTHER CELL TYPES
print('reading in ATAC cell type subset matrix')
ATAC_CT = read_ATAC(ATAC_filepath)
ATAC_CT_subsetted = ATAC_CT[res.index,:]
#ENSURE SAME CELLS, IN SAME ORDER
ATAC_CT_subsetted.obs['total_ATAC_counts'] = ATAC_CT[res.index,:].X.sum(axis=1)
print('subsetted peaks X cells matrix shape: ', ATAC_CT_subsetted.shape, res.shape)

## SUBSET to peaks open in at least threshold of all cells in cell type
ATAC_BINARIZED = anndata.AnnData(X = (ATAC_CT_subsetted.X>0)*1)
num_cells_open = list(np.array(ATAC_BINARIZED.X.sum(axis=0)))
peaks_open_enough = [idx for idx, val in enumerate(num_cells_open[0]) if val > args.threshold]
print('number of peaks open in at least', str(args.threshold), 'cells: ', len(peaks_open_enough))

# MAKE A DF OF COVARIATES:
reg_input_df = pd.DataFrame()
reg_input_df['donor'] = res['donor'].astype('category')
reg_input_df['ATAC_total_counts'] = list(ATAC_CT_subsetted.obs['total_ATAC_counts'])
reg_input_df['log_ATAC_total_counts'] = np.log1p(reg_input_df['ATAC_total_counts'])
reg_input_df[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']] = res[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']]
reg_input_df_copy = reg_input_df.copy()


## REGRESSIONS ##
peaks_subset_num = args.peaks_subset_num #40000 ##ARG
peaks = ATAC_CT_subsetted[:,peaks_open_enough].var_names[::args.peaks_subset_num]
print("num peaks:", len(peaks))

chisq_pvals = {}
regression_pvals, regression_zscores = {}, {}

for scramble in [False,True]:
    ###############################################
    if scramble == True:
        print('scrambling')
        for PC in ['PC'+str(i) for i in range(1,11)]:
            shuffled_PC_cells = []
            for donor in list(set(reg_input_df_copy['donor'])): #for each donor, shuffle cell PC loadings
                dcells = reg_input_df_copy[reg_input_df_copy['donor']==donor]
                shuffled_PC_dcells = dcells[PC][np.random.permutation(len(dcells))]
                shuffled_PC_cells.extend(shuffled_PC_dcells)
            reg_input_df_copy[PC] = shuffled_PC_cells
    ###############################################

    for peak in peaks:
        reg_input_df_copy['peak_raw_counts'] = ATAC_CT_subsetted[:,peak].X.A

        #########################################
        ### TO GET OVERALL MODEL GOODNESS OF FIT ###
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
        
        #get PC-peak pvals and zscores
        PC_peak_pvals = list(coef_df['Pr(>|z|)'][-10:])
        regression_pvals[peak] = PC_peak_pvals

        PC_peak_Zscores = list(coef_df['z value'][-10:])
        regression_zscores[peak] = PC_peak_Zscores
        print('done with peak', peak)
        #########################################

    zscores_df = pd.DataFrame(regression_zscores,index=['PC'+str(i)+' zscore' for i in range(1,11)]).T
    pvals_df = pd.DataFrame(regression_pvals,index=['PC'+str(i)+' pval' for i in range(1,11)]).T
    overall_fit_df = pd.DataFrame(chisq_pvals,index=['model_LRT_pval']).T
    
    # FDR-BH
    qvals_df = pvals_df.copy()
    for PC in ['PC'+str(i)+' pval' for i in range(1,11)]:
        #defined BH correction in regression_dependencies.py file
        qvals_df[PC] = BH(qvals_df[PC])
    qvals_df.columns=['PC'+str(i)+' qval' for i in range(1,11)]

    overall_fit_qvals_df = overall_fit_df.copy()
    overall_fit_qvals_df['model_LRT_qval'] = BH(overall_fit_qvals_df['model_LRT_pval'])

    plot_QQ(pvals_df, scramble, args.CT, 'poisson', '083122')
    
    if scramble==False:
        overall_df = pd.concat([overall_fit_qvals_df,zscores_df,pvals_df,qvals_df],axis=1)
        qvals_df_REAL = qvals_df.copy()
        zscores_df_REAL = zscores_df.copy()