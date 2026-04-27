Question 8 guidance
Inputs were validated against the ATAC sample manifest before analysis.
Consensus strategy now matches the pairwise ATAC script: union of replicate peaks per condition, then union across conditions.
Signal values come from summed bigBed peak scores over overlaps, so this clustering is exploratory rather than a substitute for BAM-based peak counting.
Use sample_hclust.pdf and sample_pca.pdf to describe how chromatin accessibility organizes the four cell types.
Use top_variable_peaks_heatmap.pdf to describe accessibility clusters that are cell-type-specific.
Then compare the overall tree structure to the RNA tree from Q4. Similar broad ordering supports coordinated regulatory and transcriptional differentiation.
