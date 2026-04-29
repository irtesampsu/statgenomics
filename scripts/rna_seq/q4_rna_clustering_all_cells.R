suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
})

source("scripts/utils/rna_utils.R")

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  cells = "CMP,CFUE,Erythroblast",
  out_dir = "results/clustering/q4_rna_all_cells",
  deseq_path = "results/deseq2/de_genes_deseq2.csv",
  limma_path = "results/limma_voom/de_genes_limma_voom.csv",
  top_n = "1000"
)
if (length(args) %% 2 != 0) stop("Arguments must be passed as --key value pairs.")
if (length(args) > 0) {
  keys <- sub("^--", "", args[seq(1, length(args), by = 2)])
  vals <- args[seq(2, length(args), by = 2)]
  for (i in seq_along(keys)) defaults[[keys[i]]] <- vals[i]
}

cells <- trimws(strsplit(as.character(defaults$cells), ",", fixed = TRUE)[[1]])
out_dir <- as.character(defaults$out_dir)
top_n <- as.integer(defaults$top_n)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sample_info <- get_all_rna_scriptseq()
sample_info <- sample_info[sample_info$condition %in% cells, , drop = FALSE]
counts <- build_rna_count_matrix(sample_info)
counts <- filter_rna_counts(counts, min_count = 10L, min_samples = 2L)

sample_info$condition <- factor(
  sample_info$condition,
  levels = cells
)
sample_info$cell_line <- sample_info$condition

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info[, c("sample", "condition", "cell_line")],
  design = ~ condition
)
rownames(colData(dds)) <- sample_info$sample

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
vsd_mat_full <- assay(vsd)  # normalized log2-scale matrix

read_results <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
clean_ensembl <- function(x) sub("\\..*$", "", x)
select_top_ids <- function(df, fc_col, padj_col, n) {
  keep <- !is.na(df[[fc_col]]) & !is.na(df[[padj_col]])
  d <- df[keep, , drop = FALSE]
  d$gene_clean <- clean_ensembl(d$gene_id)
  d <- d[order(d[[padj_col]], -abs(d[[fc_col]])), , drop = FALSE]
  unique(head(d$gene_clean, n))
}

deseq_df <- read_results(as.character(defaults$deseq_path))
limma_df <- read_results(as.character(defaults$limma_path))
selected_ids <- unique(c(
  select_top_ids(deseq_df, "log2FoldChange", "padj", top_n),
  select_top_ids(limma_df, "logFC", "adj.P.Val", top_n)
))
row_ids <- clean_ensembl(rownames(vsd_mat_full))
vsd_mat <- vsd_mat_full[row_ids %in% selected_ids, , drop = FALSE]
if (nrow(vsd_mat) == 0) stop("No overlapping top differential genes found for clustering matrix.")

sample_cor <- cor(vsd_mat)

sample_group_annot <- data.frame("Cell line" = sample_info$cell_line, row.names = sample_info$sample, check.names = FALSE)

pdf(file.path(out_dir, "sample_correlation_heatmap.pdf"), width = 7, height = 6)
pheatmap(
  sample_cor,
  annotation_col = sample_group_annot,
  annotation_row = sample_group_annot,
  main = "RNA Sample Correlation (top DE genes, log2 scale)"
)
dev.off()

pdf(file.path(out_dir, "sample_pca.pdf"), width = 7, height = 6)
pca_obj <- prcomp(t(vsd_mat), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  PC1 = pca_obj$x[, 1],
  PC2 = pca_obj$x[, 2],
  cell_line = sample_info$cell_line,
  sample = colnames(vsd_mat),
  stringsAsFactors = FALSE
)
plot(
  pca_df$PC1,
  pca_df$PC2,
  col = as.integer(pca_df$cell_line),
  pch = 19,
  xlab = "PC1",
  ylab = "PC2",
  main = "RNA PCA (top DE genes, log2 scale)"
)
text(pca_df$PC1, pca_df$PC2, labels = pca_df$sample, pos = 3, cex = 0.7)
legend("topright", legend = levels(sample_info$cell_line), col = seq_along(levels(sample_info$cell_line)), pch = 19, bty = "n")
dev.off()

pdf(file.path(out_dir, "top_differential_genes_heatmap.pdf"), width = 8, height = 9)
pheatmap(
  t(scale(t(vsd_mat))),
  annotation_col = sample_group_annot,
  show_rownames = FALSE,
  main = "RNA Heatmap (top DE genes, row-scaled log2 values)"
)
dev.off()
