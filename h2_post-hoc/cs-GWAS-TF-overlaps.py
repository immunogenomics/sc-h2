from regression_dependencies import *

### RNA and ATAC dataframes

## read in RNA and ATAC
T_cells_RNA_adata = sc.read("/data/srlab/agupta/data/d12_matrices/processed_h5ads/T_RNA_adata.h5ad")

## to get ATAC peak normalized accessibility for each cell
ATAC_T_normed = sc.read("/data/srlab/agupta/data/d12_matrices/processed_h5ads/peaksXcells_T_norm_scaled.h5ad")
ATAC_T_normed = ATAC_T_normed[list(T_cells_RNA_adata.obs_names),]

RNA_all = T_cells_RNA_adata.copy()
ATAC = ATAC_T_normed.copy()

### Checking overlap between enriched TF motifs for interesting cell states and RA

## formatting HOMER output so that enriched TF motifs --> bed files genomic coordinates

## input coords of all peaks such that we can map the enriched motifs to genomic coordinates
inpdir = '/data/srlab/agupta/data/hg19_for_HOMER_reviews/only1s/'
peakcoordsbed = pd.read_csv(inpdir+'TM_foreground_peaks_hg19.bed',sep='\t',header=None)

## read in file that has, for each motif, all of the genomic coordinates (then map to GENOME using peak coords +/- offset)
motifposs = pd.read_csv(inpdir+'output/motifMatches/TM_foreground_peaks_homer_knownMotif_motifMatches.txt',sep='\t')

## info on enriched motifs (TF name, enrichment pval) (sequence, offset from peak, and TF name) for each cell state
all_cs_TFs_df = pd.DataFrame()
for cs in ['T-3','T-7','T-8','T-9']:
    print(cs)
    print('reading in motifs info')
    csmotifs = pd.read_csv(inpdir+'output/'+cs+'_condOnI_enriched_TFs/knownResults.txt',sep='\t')
    #only keep sig ones
    print('subsetting to only those TFs sig enriched in this cell state')
    csmotifs = csmotifs[csmotifs['P-value']<0.0001]
    csmotifs['TF'] = [i.split('(')[0] for i in csmotifs['Motif Name']]
    enriched_TFs = list(csmotifs['Motif Name'])
    enriched_TFs_withMotifPos = motifposs[motifposs['Motif Name'].isin(enriched_TFs)];
    print(len(motifposs), len(enriched_TFs_withMotifPos), len(enriched_TFs_withMotifPos)/len(motifposs));
    enriched_TFs_withMotifPos['PositionID'] = [i-1 for i in enriched_TFs_withMotifPos['PositionID']]
    enriched_TFs_withMotifPos
    print('mapping TFs to genomic coord windows')
    enriched_TFs_withMotifPos['peak_chr'] = list(peakcoordsbed.iloc[list(enriched_TFs_withMotifPos['PositionID']),:][0])
    enriched_TFs_withMotifPos['peak_start'] = list(peakcoordsbed.iloc[list(enriched_TFs_withMotifPos['PositionID']),:][1])
    enriched_TFs_withMotifPos['peak_end'] = list(peakcoordsbed.iloc[list(enriched_TFs_withMotifPos['PositionID']),:][2])
    enriched_TFs_withMotifPos['peak_center'] = enriched_TFs_withMotifPos['peak_end']-600
    enriched_TFs_withMotifPos['motif_start'] = enriched_TFs_withMotifPos['peak_center'] + enriched_TFs_withMotifPos['Offset']
    enriched_TFs_withMotifPos['motif_end'] = enriched_TFs_withMotifPos['motif_start'] + enriched_TFs_withMotifPos['Sequence'].str.len()
    print('creating final DF for overlap with GWAS SNPs with trait of interest')
    enriched_TFs_withMotifPos = enriched_TFs_withMotifPos[['Motif Name','peak_chr','motif_start','motif_end']]
    enriched_TFs_withMotifPos.index = enriched_TFs_withMotifPos['Motif Name']
    enriched_TFs_withMotifPos = enriched_TFs_withMotifPos.drop('Motif Name',axis=1)
    enriched_TFs_withMotifPos.columns = ['chr','start','end']
    enriched_TFs_withMotifPos['state'] = [cs]*len(enriched_TFs_withMotifPos)
    enriched_TFs_withMotifPos.index = [i.split('(')[0] for i in enriched_TFs_withMotifPos.index]
    enriched_TFs_withMotifPos['row'] = np.arange(0,len(enriched_TFs_withMotifPos),1)
    enriched_TFs_withMotifPos['fullmotifrange_CS'] = enriched_TFs_withMotifPos['chr']+":"+enriched_TFs_withMotifPos['start'].astype(str)+"-"+enriched_TFs_withMotifPos['end'].astype(str)+"_"+enriched_TFs_withMotifPos['state']
    all_cs_TFs_df = pd.concat([all_cs_TFs_df,enriched_TFs_withMotifPos],axis=0)

### GWAS SNPs for our AAI diseases df

GWAS_df = pd.DataFrame(pd.read_csv("/data/srlab/agupta/data/ref_files/AAI_GWAS.csv",index_col=0))

GWAS_df.index = ['chr'+str(GWAS_df['CHR_ID'][i]) + ":" + str(GWAS_df['CHR_POS'][i]) for i in range(len(GWAS_df))]
GWAS_df = GWAS_df[['JOURNAL','LINK','DISEASE/TRAIT','INITIAL SAMPLE SIZE','CHR_ID','CHR_POS','REPORTED GENE(S)','P-VALUE']]

GWAS_to_overlap = GWAS_df[['CHR_ID','CHR_POS','P-VALUE','DISEASE/TRAIT','REPORTED GENE(S)','LINK']]
GWAS_to_overlap.columns = ['chr','SNP_position','p-value','Disease/Trait','Reported Gene(s)','Study']
GWAS_to_overlap['chr'] = ['chr'+str(i) for i in GWAS_to_overlap['chr']]
GWAS_to_overlap['SNP'] = GWAS_to_overlap.index
GWAS_to_overlap.index = np.arange(0,len(GWAS_to_overlap),1)

RAtraits = ['Rheumatoid arthritis',
 'Rheumatoid arthritis (ACPA-negative)',
 'Rheumatoid arthritis (ACPA-positive)',
 'Rheumatoid arthritis (rheumatoid factor and/or anti-cyclic citrullinated peptide seropositive)']
GWAS_to_overlap = GWAS_to_overlap[GWAS_to_overlap['Disease/Trait'].isin(RAtraits)] #GUTtraits RAtraits

### get peaks for important cell states

states_use = T_states
folder='/data/srlab/agupta/data/h2/bed/d12/subtypes_hg38/'

impt_states = pd.DataFrame()
for state in states_use:
    state_peaks_df = pd.DataFrame(pd.read_csv(folder+'083122_'+state+'_14kPeaks_binary_hg38.bed',sep='\t',header=None)))
    state_peaks_df['state'] = state
    state_peaks_df = state_peaks_df[state_peaks_df[3]==1]
    impt_states = pd.concat([impt_states, state_peaks_df], axis=0)

impt_states = impt_states[[0,1,2,'state']]
impt_states.columns = ['chr','start','end','state'] #these peaks are already 1.2k in window size
impt_states['start']=[i+500 for i in impt_states['start']]
impt_states['end']=[i-500 for i in impt_states['end']]
impt_states

### get peak set + GWAS SNP overlaps

## if doing TF motif and GWAS overlap
impt_states = all_cs_TFs_df

from intervaltree import IntervalTree, Interval

peak_gene_mapping_dict = dict()
cell_states = []
ALL_assoc_SNPS = []
ALL_assoc_peaks = []
diseases = []
pvals = []
studies = []
all_genes = []
TFs = []
output_df = pd.DataFrame()

for cluster in states_use:
    print(cluster)
    this_peak_set = impt_states[impt_states['state']==cluster]#[::10]
    peak_set = this_peak_set[['chr','start','end']]

    overlapping_genes_set = set()
    for chrnum in np.arange(1,23,1):
        print(chrnum)
        chrname = 'chr'+str(chrnum)
        #subset peaks to this chr
        chr_DYN_peaks_df = peak_set[peak_set['chr']==chrname]
        #subset SNPs to this chr
        chr_TSS_df = GWAS_to_overlap[GWAS_to_overlap['chr']==chrname]
        # get all GWAS SNPs in this chromosome
        SNPs_check = chr_TSS_df['SNP_position']
        indices_for_SNPs_check = chr_TSS_df['SNP_position'].index

        # print('looping through peaks of interest, building IntervalTree')
        peaks_tree = IntervalTree()
        for peak_i in range(len(chr_DYN_peaks_df)):
            try:
                peakstart,peakend = chr_DYN_peaks_df.iloc[peak_i,:]['start'],chr_DYN_peaks_df.iloc[peak_i,:]['end']
                peaks_tree[peakstart:peakend]=chrname+":"+str(peakstart)+"-"+str(peakend)
            except:
                pass

        ## GETTING OVERLAPPING SNPS
        assoc_SNPS = []
        for SNP_i in indices_for_SNPs_check:
            SNP = chr_TSS_df['SNP_position'][SNP_i]
            if peaks_tree.overlaps(SNP):
                pval=chr_TSS_df['p-value'][SNP_i]
                disease=chr_TSS_df['Disease/Trait'][SNP_i]
                genes=chr_TSS_df['Reported Gene(s)'][SNP_i]
                study=chr_TSS_df['Study'][SNP_i]
                SNP_report=chr_TSS_df['SNP'][SNP_i]
                print(SNP_report)
                assoc_SNPS.append(SNP)
                ##add values associated with this SNP/disease/study to eventual DF
                ALL_assoc_SNPS.append(SNP_report)
                diseases.append(disease)
                cell_states.append(cluster)
                pvals.append(pval)
                studies.append(study)
                all_genes.append(genes)

        #get the peaks corresponding to above SNP(s)
        assoc_peaks = []
        for i in peaks_tree:
            for SNP in assoc_SNPS:
                # print(chrnum, SNP)
                if i[0] <= SNP and i[1] >= SNP:
                    assoc_peakname = 'chr'+str(chrnum)+':'+str(i[0])+'-'+str(i[1])
                    #to get row index to map motif window to TF name:
                    peak_cs_name = assoc_peakname+"_"+cs
                    print(assoc_peakname)
                    assoc_peaks.append(assoc_peakname)
                    overlapping_TFs = ', '.join(list(enriched_TFs_withMotifPos[enriched_TFs_withMotifPos['fullmotifrange_CS']==peak_cs_name].index))
                    TFs.append(overlapping_TFs)
        ALL_assoc_peaks.extend(assoc_peaks)

output_df['cell state'] = cell_states
# output_df['TF(s)'] = TFs
output_df['Disease/Trait'] = diseases
output_df['SNP'] = ALL_assoc_SNPS
output_df['p-value'] = pvals
output_df['Reported Gene(s)'] = all_genes
output_df['Study'] = studies
output_df.index = output_df['cell state']
output_df = output_df.drop('cell state',axis=1)