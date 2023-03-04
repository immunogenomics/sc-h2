from regression_dependencies import *

ct='T'

reg_outs = pd.DataFrame(pd.read_csv('/data/srlab/agupta/data/peak_gene_scores/d12_scores/'+ct+'/083122_poisson_regression_10PC_all_outputs.csv',
                        index_col=0))
uniform_pvals = np.random.uniform(low=0, high=1,size=len(reg_outs))
second_plot = reg_outs['model_LRT_pval']
xmax=np.max(-np.log10(uniform_pvals))

# for real
plt.figure(figsize=(5,8))
plt.scatter(sorted(-np.log10(uniform_pvals)),
           sorted([-np.log10(i) for i in second_plot]),color='#F892EB',s=40)

# for null
reg_outs_NULL = pd.DataFrame(pd.read_csv('/data/srlab/agupta/data/peak_gene_scores/d12_scores/T/083122_poisson_regression_TCELL_10PC_all_outputs_shuffled.csv',
                        index_col=0))
uniform_pvals_NULL = np.random.uniform(low=0, high=1,size=len(reg_outs_NULL))
second_plot_NULL = reg_outs_NULL['model_LRT_pval']
plt.scatter(sorted(-np.log10(uniform_pvals_NULL)),
            sorted([-np.log10(i) for i in second_plot_NULL]),color='grey',s=40)
plt.plot([0,xmax],[0,xmax], color='black', linestyle='--')
sns.despine()
plt.ylabel('-log10p(observed)',fontsize=30)
plt.xlabel('-log10p(expected)',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
inflation_val = inflation(np.array(list(second_plot)))
plt.title(ct+' signal = '+str(np.round(inflation_val,1)),fontsize=30)
print('null inflation:', inflation(np.array(list(second_plot_NULL))))

#colors for F M T E B: ['#00BFC4','#619CFF','#F892EB','#00BA38','#F8766D']
