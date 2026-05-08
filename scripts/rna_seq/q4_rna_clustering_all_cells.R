suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

source("scripts/utils/rna_utils.R")

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  cells = "CMP,CFUE,HSC,Erythroblast",
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

pca_obj <- stats::prcomp(t(vsd_mat), center = TRUE, scale. = FALSE)
var_expl <- summary(pca_obj)$importance[2, ] * 100
pca_df <- data.frame(
  PC1 = pca_obj$x[, 1],
  PC2 = pca_obj$x[, 2],
  cell_line = sample_info$cell_line,
  sample = colnames(vsd_mat),
  stringsAsFactors = FALSE
)
pca_df$cell_line <- factor(pca_df$cell_line, levels = cells)
pal <- stats::setNames(
  grDevices::colorRampPalette(c("#2166AC", "#FDAE61", "#B2182B"))(length(cells)),
  cells
)
xr <- range(pca_df$PC1)
yr <- range(pca_df$PC2)
pad_x <- max(diff(xr) * 0.28, 0.5)
pad_y <- max(diff(yr) * 0.28, 0.5)
x_bounds <- c(xr[1] - pad_x, xr[2] + pad_x)
y_bounds <- c(yr[1] - pad_y, yr[2] + pad_y)
xl <- sprintf("PC1 (%.1f%% variance)", var_expl[1])
yl <- sprintf("PC2 (%.1f%% variance)", var_expl[2])

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = cell_line)) +
  geom_point(size = 3.6, alpha = 0.92) +
  geom_text_repel(
    aes(label = sample),
    size = 3,
    segment.size = 0.28,
    segment.color = "grey55",
    box.padding = 0.45,
    point.padding = 0.28,
    min.segment.length = 0,
    max.overlaps = Inf,
    xlim = x_bounds,
    ylim = y_bounds,
    show.legend = FALSE
  ) +
  scale_color_manual(values = pal, drop = FALSE, name = "Cell line") +
  coord_cartesian(
    xlim = x_bounds,
    ylim = y_bounds,
    expand = FALSE,
    clip = "on"
  ) +
  labs(title = "RNA-seq PCA (top DE genes, vst)", x = xl, y = yl) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = margin(14, 18, 20, 18),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4)))

ggsave(file.path(out_dir, "sample_pca.pdf"), p_pca, width = 8, height = 7.2)

pdf(file.path(out_dir, "top_differential_genes_heatmap.pdf"), width = 8, height = 9)
pheatmap(
  t(scale(t(vsd_mat))),
  annotation_col = sample_group_annot,
  show_rownames = FALSE,
  main = "RNA Heatmap (top DE genes, row-scaled log2 values)"
)
dev.off()
