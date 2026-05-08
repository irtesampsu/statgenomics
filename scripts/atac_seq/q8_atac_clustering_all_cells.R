suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

source("scripts/utils/atac_utils.R")

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  cells = "CMP,CFUE,HSC,Erythroblast",
  out_dir = "results/clustering/q8_atac_all_cells",
  dar_deseq2 = "results/atac_seq/pairwise_cmp_vs_ery/DARs_consensus_deseq2.csv",
  dar_limma = "results/atac_seq/pairwise_cmp_vs_ery/DARs_consensus_limma_voom.csv",
  top_n = "2000"
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
sample_info <- get_all_atac()
sample_info <- sample_info[sample_info$condition %in% cells, , drop = FALSE]
gr_list <- import_bigbed_replicates(sample_info)
all_peaks <- build_union_consensus_peaks(gr_list)
atac_quant <- build_atac_count_matrix(sample_info, consensus = all_peaks)
count_mat <- atac_quant$count_mat

condition <- factor(
  sample_info$condition,
  levels = cells
)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = data.frame(row.names = colnames(count_mat), condition),
  design = ~ condition
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
vsd_mat_full <- assay(vsd)  # normalized log2-scale matrix
peak_keys <- paste(as.character(seqnames(all_peaks)), start(all_peaks), end(all_peaks), sep = ":")
rownames(vsd_mat_full) <- peak_keys

read_results <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
select_top_peaks <- function(path, n) {
  if (!file.exists(path)) return(character(0))
  d <- read_results(path)
  if (!all(c("chr", "start", "end", "log2FoldChange", "padj") %in% colnames(d))) return(character(0))
  d <- d[!is.na(d$log2FoldChange) & !is.na(d$padj), , drop = FALSE]
  d <- d[order(d$padj, -abs(d$log2FoldChange)), , drop = FALSE]
  unique(head(paste(d$chr, d$start, d$end, sep = ":"), n))
}
selected_peaks <- unique(c(
  select_top_peaks(as.character(defaults$dar_deseq2), top_n),
  select_top_peaks(as.character(defaults$dar_limma), top_n)
))
vsd_mat <- vsd_mat_full[rownames(vsd_mat_full) %in% selected_peaks, , drop = FALSE]
if (nrow(vsd_mat) == 0) stop("No overlapping top differential peaks found for clustering matrix.")

sample_cor <- cor(vsd_mat)

annotation_df <- data.frame(condition = condition, row.names = colnames(count_mat))

pdf(file.path(out_dir, "sample_correlation_heatmap.pdf"), width = 7, height = 6)
pheatmap(
  sample_cor,
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  main = "ATAC Sample Correlation (top DAR peaks, log2 scale)"
)
dev.off()

pca_obj <- stats::prcomp(t(vsd_mat), center = TRUE, scale. = FALSE)
var_expl <- summary(pca_obj)$importance[2, ] * 100
pca_df <- data.frame(
  PC1 = pca_obj$x[, 1],
  PC2 = pca_obj$x[, 2],
  condition = condition,
  sample = colnames(vsd_mat),
  stringsAsFactors = FALSE
)
pca_df$condition <- factor(pca_df$condition, levels = cells)
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

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
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
  labs(title = "ATAC-seq PCA (top DAR peaks, vst)", x = xl, y = yl) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    plot.margin = margin(14, 18, 20, 18),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4)))

ggsave(file.path(out_dir, "sample_pca.pdf"), p_pca, width = 8, height = 7.2)

pdf(file.path(out_dir, "top_differential_peaks_heatmap.pdf"), width = 8, height = 9)
pheatmap(
  t(scale(t(vsd_mat))),
  annotation_col = annotation_df,
  show_rownames = FALSE,
  main = "ATAC Heatmap (top DAR peaks, row-scaled log2 values)"
)
dev.off()
