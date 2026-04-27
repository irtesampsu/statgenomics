suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(matrixStats)
  library(RColorBrewer)
})

source("scripts/utils/rna_utils.R")

out_dir <- "results/clustering/pairwise_cmp_vs_ery"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_info <- get_pairwise_rna_scriptseq("CMP", "Erythroblast")
counts <- build_rna_count_matrix(sample_info)
counts <- filter_rna_counts(counts, min_count = 10L, min_samples = 2L)

sample_info$condition <- factor(sample_info$condition, levels = c("CMP", "Erythroblast"))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info[, c("sample", "condition")],
  design = ~ condition
)
rownames(colData(dds)) <- sample_info$sample

vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

gene_vars <- rowVars(vsd_mat)
top_n <- min(1000L, length(gene_vars))
top_ids <- names(sort(gene_vars, decreasing = TRUE))[seq_len(top_n)]
top_mat <- vsd_mat[top_ids, , drop = FALSE]
top_mat_scaled <- t(scale(t(top_mat)))

annotation_col <- data.frame(
  Condition = sample_info$condition,
  row.names = sample_info$sample
)
annotation_colors <- list(
  Condition = c(CMP = "#67A9CF", Erythroblast = "#EF8A62")
)

png(
  filename = file.path(out_dir, "Gene_heatmap_top1000.png"),
  width = 800,
  height = 1000,
  res = 150
)
pheatmap(
  top_mat_scaled,
  show_rownames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  clustering_method = "complete",
  cutree_rows = 4,
  main = "Top 1000 Variable Genes - Hierarchical Clustering\nCMP vs Erythroblast"
)
dev.off()

writeLines(
  c(
    "Top-1000 variable gene heatmap (pairwise CMP vs Erythroblast)",
    sprintf("Samples: %s", paste(sample_info$sample, collapse = ", ")),
    sprintf("Genes plotted: %d", top_n),
    "Matrix: VST-transformed counts with per-gene z-score scaling"
  ),
  con = file.path(out_dir, "Gene_heatmap_top1000_notes.txt")
)
