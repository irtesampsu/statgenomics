ATAC pairwise analysis notes
Inputs: processed bigBed peak files only
Consensus strategy: union of replicate peaks per condition, then union across conditions
Quantification: summed bigBed peak scores over overlaps per consensus peak
Limitation: DESeq2 is being applied to score-derived pseudo-counts, not BAM-derived read counts.
Interpretation: treat differential accessibility statistics as approximate unless BAM-based counting is added.
