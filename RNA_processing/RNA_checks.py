from regression_dependencies import *

CT='T'
RNA = sc.read_h5ad("/data/srlab/agupta/data/d12_matrices/processed_h5ads/"+CT+"_RNA_adata.h5ad")

## number of cells per donor
donorcolors=['lightgrey']*10
donors, total_cells = [], []
for d in list(sorted(set(RNA.obs['donor']))):
    RNA_DONOR = RNA[RNA.obs['donor']==d]
    donors.append(d)
    total_cells.append(len(RNA_DONOR))
plt.figure(figsize=(10,5))
plt.bar(donors, total_cells, color=donorcolors)
sns.despine()
plt.xticks("", rotation=90)
plt.xlabel("")
plt.ylabel("# cells")
plt.show()

## number of cells per subtype
ct='t'
if ct=='t':
    subtypecolors=['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928', '#B3E2CD', '#FDCDAC', '#CBD5E8', '#F4CAE4', '#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666']
if ct=='m':
    subtypecolors=['#d0b4dc', '#FCCDE5', '#945cb4','#842bd7','yellow','#B38072','#d11141','#FDB462','grey','#1F78B4','#A6CEE3','#66C2A4','#CCECE6','#238B45','#A1D99B']
if ct=='b':
    subtypecolors=['#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#6A3D9A','#CAB2D6','#FDBF6F','brown','#FF7F00']
if ct=='e':
    subtypecolors=['#8DD3C7','#BEBADA','#FB8072','#80B1D3','#FDB462']
if ct=='f':
    subtypecolors=['#A6CEE3','#FDBF6F','#B2DF8A','#33A02C','#FB9A99','#E31A1C','#1F78B4', '#FF7F00', '#6495ED']

sts, total_cells = [], []
for st in list(sorted(set(RNA.obs['ct_subtype']))):
    RNA_ST = RNA[RNA.obs['ct_subtype']==st]
    # print(st, len(RNA_ST))
    sts.append(st)
    total_cells.append(len(RNA_ST))
plt.figure(figsize=(10,5))
plt.bar(sts, total_cells, color=subtypecolors)
sns.despine()
plt.xticks(rotation=90)
plt.ylabel("# cells")
plt.xlabel("cell subtype")

#cell type UMAP colored by cell states
RNA.uns['ct_subtype_colors'] = subtypecolors
sc.pl.umap(RNA,color='ct_subtype',size=300,palette=RNA.uns['ct_subtype_colors'],frameon=False)


## check variance explained by each of the top PCs
numPCs=10
ct_colors = ['#F8766D','#F892EB','#619CFF','#00BFC4','#00BA38']
i=0
for ct in ['B','T','myeloid','fibroblast','endothelial']:
    RNA_forPCs = sc.read_h5ad("/data/srlab/agupta/data/d12_matrices/processed_h5ads/"+ct+"_RNA_adata.h5ad")
    plt.figure(figsize=(5,3))
    plt.scatter(np.arange(0,numPCs), 
                [i*100 for i in RNA_forPCs.uns['pca']['variance_ratio']][:numPCs],
                color=ct_colors[i],s=100)
    plt.ylabel("% variance explained")
    plt.ylim(0,7)
    plt.xticks(np.arange(0,10,1),np.arange(1,11,1))
    plt.xlabel("PC")
    sns.despine()
    plt.show()
    i+=1