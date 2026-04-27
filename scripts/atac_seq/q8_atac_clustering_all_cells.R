suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(DESeq2)
  library(pheatmap)
})

source("scripts/utils/atac_utils.R")

dir.create("results/clustering/q8_atac_all_cells", recursive = TRUE, showWarnings = FALSE)

sample_info <- get_all_atac()
gr_list <- import_bigbed_replicates(sample_info)
all_peaks <- build_union_consensus_peaks(gr_list)
atac_quant <- build_atac_count_matrix(sample_info, consensus = all_peaks)
count_mat <- atac_quant$count_mat
quant_method <- atac_quant$quant_method

condition <- factor(
  sample_info$condition,
  levels = c("HSC", "CMP", "CFUE", "Erythroblast")
)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = data.frame(row.names = colnames(count_mat), condition),
  design = ~ condition
)

vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

sample_cor <- cor(vsd_mat)
write.csv(sample_cor, "results/clustering/q8_atac_all_cells/sample_correlation.csv")

annotation_df <- data.frame(condition = condition, row.names = colnames(count_mat))

pdf("results/clustering/q8_atac_all_cells/sample_correlation_heatmap.pdf", width = 7, height = 6)
pheatmap(
  sample_cor,
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  main = "ATAC-seq Sample Correlation"
)
dev.off()

sample_dist <- dist(t(vsd_mat))
hc <- hclust(sample_dist)

pdf("results/clustering/q8_atac_all_cells/sample_hclust.pdf", width = 8, height = 6)
plot(hc, main = "Q8 ATAC-seq Hierarchical Clustering")
dev.off()

pdf("results/clustering/q8_atac_all_cells/sample_pca.pdf", width = 7, height = 6)
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

peak_var <- apply(vsd_mat, 1, var)
top_peak_ids <- order(peak_var, decreasing = TRUE)[seq_len(min(100, length(peak_var)))]
top_peak_mat <- vsd_mat[top_peak_ids, , drop = FALSE]

pdf("results/clustering/q8_atac_all_cells/top_variable_peaks_heatmap.pdf", width = 8, height = 10)
pheatmap(
  top_peak_mat,
  scale = "row",
  show_rownames = FALSE,
  annotation_col = annotation_df,
  main = "Top Variable ATAC Peaks"
)
dev.off()

cluster_summary <- c(
  "Question 8 guidance",
  "Inputs were validated against the ATAC sample manifest before analysis.",
  "Consensus strategy now matches the pairwise ATAC script: union of replicate peaks per condition, then union across conditions.",
  if (quant_method == "bam_counts") {
    "Quantification uses BAM read counts per consensus peak, which is preferred for clustering and downstream differential analysis."
  } else {
    "Signal values come from summed bigBed peak scores over overlaps, so this clustering is exploratory rather than a substitute for BAM-based peak counting."
  },
  "Use sample_hclust.pdf and sample_pca.pdf to describe how chromatin accessibility organizes the four cell types.",
  "Use top_variable_peaks_heatmap.pdf to describe accessibility clusters that are cell-type-specific.",
  "Then compare the overall tree structure to the RNA tree from Q4. Similar broad ordering supports coordinated regulatory and transcriptional differentiation."
)
writeLines(cluster_summary, "results/clustering/q8_atac_all_cells/README.txt")
