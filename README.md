## Mapping cell state heritability in autoimmunity with single-cell data

Code related to the manuscript "Dynamic regulatory elements in single-cell multimodal data implicate key immune cell states enriched for autoimmune disease heritability" (Accepted in principle at Nature Genetics).

Preprint now available at https://www.medrxiv.org/content/10.1101/2023.02.24.23286364v1

All scripts required to reproduce the analysis in this paper can also be found at: https://doi.org/10.5281/zenodo.8436010

Contact: Anika Gupta anikagupta@g.harvard.edu

To follow along with the analyses described in the manuscript, we would recommend using the scripts in the following directory order:

1. RNA processing (preparing the counts data for the regressions; also includes a subdirectory with the QC and ATAC processing code)
2. regressions (underlying basis for this work and the key novelty)
3. dynamic_peaks_characterization (plots for Figure 2)
4. stat_gen (preparing the regression outputs for heritability analyses)
5. ldsc (running heritability enrichment analyses)
6. h2_post-hoc (making sense of and plotting the heritability outputs; for h2 plots in Figures 3-6)
7. benchmarking (trying out modified versions of this approach OR the approach as-is but on PBMCs)
