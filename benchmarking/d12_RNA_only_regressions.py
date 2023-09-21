from regression_dependencies import *
import argparse

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--CT', type=str)
parser.add_argument('--threshold', type=int, default=50)
parser.add_argument('--genes_subset_num', type=int, default=10000)
parser.add_argument('--family', type=str)
parser.add_argument('--date', type=str)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

print('reading in covariates df')
res = pd.read_csv("/data/srlab/agupta/data/d12_matrices/cellPC-loadings/"+args.CT+"_RNA_PC_loadings.csv",index_col=0)
# read in peaks X cells DF
RNA_filepath = "/data/srlab/agupta/data/d12_matrices/h5ads/genesXcells_"+args.CT+".h5ad"
print('reading in RNA cell type subset matrix')
RNA_CT = sc.read_h5ad(RNA_filepath)
RNA_CT_subsetted = RNA_CT[res.index,:]
#ensure same cells, in same order
RNA_CT_subsetted.obs['total_RNA_counts'] = RNA_CT[res.index,:].X.sum(axis=1)
print('subsetted genes X cells matrix shape: ', RNA_CT_subsetted.shape, res.shape)

## SUBSET to peaks open in at least threshold of all cells in cell type
RNA_BINARIZED = anndata.AnnData(X = (RNA_CT_subsetted.X>0)*1)
num_cells_open = list(np.array(RNA_BINARIZED.X.sum(axis=0)))
genes_expressed_enough = [idx for idx, val in enumerate(num_cells_open[0]) if val > args.threshold]
print('number of genes expressed in at least', str(args.threshold), 'cells: ', len(genes_expressed_enough))

# MAKE A DF OF COVARIATES:
reg_input_df = pd.DataFrame()
reg_input_df['donor'] = res['donor'].astype('category')
reg_input_df['RNA_total_counts'] = list(RNA_CT_subsetted.obs['total_RNA_counts'])
reg_input_df['log_RNA_total_counts'] = np.log1p(reg_input_df['RNA_total_counts'])
reg_input_df[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']] = res[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']]
reg_input_df_copy = reg_input_df.copy()


## REGRESSIONS ##
genes_subset_num = args.genes_subset_num
genes = RNA_CT_subsetted[:,genes_expressed_enough].var_names[::args.genes_subset_num]
print("num genes:", len(genes))

chisq_pvals = {}
regression_pvals, regression_zscores = {}, {}

for scramble in [False,True]:
    if scramble == True:
        print('scrambling')
        for PC in ['PC'+str(i) for i in range(1,11)]:
            shuffled_PC_cells = []
            for donor in list(set(reg_input_df_copy['donor'])): #for each donor, shuffle cell PC loadings
                dcells = reg_input_df_copy[reg_input_df_copy['donor']==donor]
                shuffled_PC_dcells = dcells[PC][np.random.permutation(len(dcells))]
                shuffled_PC_cells.extend(shuffled_PC_dcells)
            reg_input_df_copy[PC] = shuffled_PC_cells

    for gene in genes: ##iterate over genes since we only have RNA to work with
        reg_input_df_copy['gene_raw_counts'] = RNA_CT_subsetted[:,gene].X.A

        ## to get overall variance explained
        print('fitting & predicting null and full models')
        #null model
        null_formula = "gene_raw_counts ~  offset(log(RNA_total_counts)) + donor"
        #full model
        full_formula = "gene_raw_counts ~  offset(log(RNA_total_counts)) + donor + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"

        #yhat_null
        nullmodelfit, yhat_null, _ = get_glm_fit_predict(null_formula, reg_input_df_copy, 'poisson')
        #yhat_alt
        ### run PCs-peak multiple regression; coef_df: to get each covariate's info
        altmodelfit, yhat_alt, coef_df = get_glm_fit_predict(full_formula, reg_input_df_copy, 'poisson')

        #compare null from full model
        print('comparing null vs full models for model fit')
        # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova.glm
        chisq_pval = get_LRT(nullmodelfit, altmodelfit)
        chisq_pvals[gene] = chisq_pval
        
    overall_fit_df = pd.DataFrame(chisq_pvals,index=['model_LRT_pval']).T
    
    # FDR-BH
    overall_fit_qvals_df = overall_fit_df.copy()
    overall_fit_qvals_df['model_LRT_FDR'] = BH(overall_fit_qvals_df['model_LRT_pval'])

    overall_df = pd.concat([overall_fit_qvals_df,zscores_df,pvals_df,qvals_df],axis=1)