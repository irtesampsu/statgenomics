suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(DESeq2)
  library(pheatmap)
})

dir.create("results/clustering/q8_atac_all_cells", recursive = TRUE, showWarnings = FALSE)

files <- list(
  HSC = c("data/atac_seq/HSC_rep1_ENCFF662DYG.bigBed", "data/atac_seq/HSC_rep2_ENCFF255IVU.bigBed"),
  CMP = c("data/atac_seq/CMP_rep1_ENCFF832UUS.bigBed", "data/atac_seq/CMP_rep2_ENCFF343PTQ.bigBed"),
  CFUE = c("data/atac_seq/CFUE_rep1_ENCFF796ZSB.bigBed", "data/atac_seq/CFUE_rep2_ENCFF599ZDJ.bigBed"),
  Erythroblast = c("data/atac_seq/Erythroblast_rep1_ENCFF181AMY.bigBed", "data/atac_seq/Erythroblast_rep2_ENCFF616EWK.bigBed")
)

missing_files <- unlist(files)[!file.exists(unlist(files))]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Q8 ATAC clustering is missing required bigBed files:",
      paste(missing_files, collapse = ", ")
    )
  )
}

gr_list <- lapply(files, function(reps) {
  lapply(reps, function(f) {
    gr <- import(f, format = "bigBed")
    seqlevelsStyle(gr) <- "UCSC"
    gr
  })
})

consensus_per_condition <- lapply(gr_list, function(reps) {
  reduce(intersect(reps[[1]], reps[[2]]))
})

all_peaks <- Reduce(c, consensus_per_condition)
all_peaks <- reduce(all_peaks)

get_signal <- function(gr, consensus) {
  hits <- findOverlaps(consensus, gr)
  signal <- numeric(length(consensus))
  if ("score" %in% colnames(mcols(gr))) {
    signal[queryHits(hits)] <- signal[queryHits(hits)] + mcols(gr)$score[subjectHits(hits)]
  } else {
    signal[queryHits(hits)] <- signal[queryHits(hits)] + 1
  }
  signal
}

all_samples <- do.call(c, gr_list)
count_mat <- sapply(all_samples, get_signal, consensus = all_peaks)
count_mat <- round(as.matrix(count_mat))
colnames(count_mat) <- c(
  "HSC_rep1", "HSC_rep2",
  "CMP_rep1", "CMP_rep2",
  "CFUE_rep1", "CFUE_rep2",
  "Erythroblast_rep1", "Erythroblast_rep2"
)

condition <- factor(
  c("HSC", "HSC", "CMP", "CMP", "CFUE", "CFUE", "Erythroblast", "Erythroblast"),
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
  "Use sample_hclust.pdf and sample_pca.pdf to describe how chromatin accessibility organizes the four cell types.",
  "Use top_variable_peaks_heatmap.pdf to describe accessibility clusters that are cell-type-specific.",
  "Then compare the overall tree structure to the RNA tree from Q4. Similar broad ordering supports coordinated regulatory and transcriptional differentiation."
)
writeLines(cluster_summary, "results/clustering/q8_atac_all_cells/README.txt")
