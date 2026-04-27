suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

source("scripts/utils/rna_utils.R")

dir.create("results/clustering/q4_rna_all_cells", recursive = TRUE, showWarnings = FALSE)

sample_info <- get_all_rna_scriptseq()
counts <- build_rna_count_matrix(sample_info)
counts <- filter_rna_counts(counts, min_count = 10L, min_samples = 2L)

sample_info$condition <- factor(
  sample_info$condition,
  levels = c("HSC", "CMP", "CFUE", "Erythroblast")
)
sample_info$cell_line <- sample_info$condition

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info[, c("sample", "condition", "cell_line")],
  design = ~ condition
)
rownames(colData(dds)) <- sample_info$sample

vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

sample_cor <- cor(vsd_mat)
write.csv(sample_cor, "results/clustering/q4_rna_all_cells/sample_correlation.csv")

sample_group_annot <- data.frame("Cell line" = sample_info$cell_line, row.names = sample_info$sample, check.names = FALSE)

pdf("results/clustering/q4_rna_all_cells/sample_correlation_heatmap.pdf", width = 7, height = 6)
pheatmap(
  sample_cor,
  annotation_col = sample_group_annot,
  annotation_row = sample_group_annot,
  main = "RNA-seq Sample Correlation"
)
dev.off()

sample_dist <- dist(t(vsd_mat))
hc <- hclust(sample_dist)

pdf("results/clustering/q4_rna_all_cells/sample_hclust.pdf", width = 8, height = 6)
plot(hc, main = "Q4 RNA-seq Hierarchical Clustering")
dev.off()

pdf("results/clustering/q4_rna_all_cells/sample_pca.pdf", width = 7, height = 6)
print(plotPCA(vsd, intgroup = "cell_line"))
dev.off()

gene_var <- apply(vsd_mat, 1, var)
top_gene_ids <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(50, length(gene_var)))]
top_gene_mat <- vsd_mat[top_gene_ids, , drop = FALSE]
top_gene_scaled <- t(scale(t(top_gene_mat)))
top_gene_scaled[is.na(top_gene_scaled)] <- 0

ensembl_ids <- sub("\\..*$", "", rownames(top_gene_scaled))
gene_symbols <- suppressMessages(
  AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
)
gene_labels <- ifelse(is.na(gene_symbols) | !nzchar(gene_symbols), ensembl_ids, unname(gene_symbols))
rownames(top_gene_scaled) <- make.unique(gene_labels)

cell_line_levels <- levels(sample_info$condition)
cell_line_colors <- c(
  HSC = "#4DAF4A",
  CMP = "#E78AC3",
  CFUE = "#4DBBD5",
  Erythroblast = "#8A6BBE"
)
cell_line_colors <- cell_line_colors[cell_line_levels]

ha_top <- HeatmapAnnotation(
  `Cell line` = sample_info$condition,
  col = list(`Cell line` = cell_line_colors),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  show_annotation_name = TRUE
)

heat_cols <- colorRamp2(c(-2, 0, 2), c("#2C7BB6", "white", "#D7191C"))

pdf("results/clustering/q4_rna_all_cells/top_variable_genes_heatmap.pdf", width = 8, height = 10)
ht <- Heatmap(
  top_gene_scaled,
  name = "Expression\n(z-score)",
  col = heat_cols,
  top_annotation = ha_top,
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 6),
  show_column_names = TRUE,
  column_title = "Top Variable RNA-seq Genes",
  cluster_rows = TRUE,
  cluster_columns = TRUE
)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "top")
dev.off()

cluster_summary <- c(
  "Question 4 guidance",
  "Inputs were validated against the ScriptSeq sample manifest before analysis.",
  "Heatmap interpretation note: the top annotation bar is sample group identity (cell line), while the heatmap body is row-scaled expression intensity (relative high/low per gene, not raw counts).",
  "Use sample_hclust.pdf and sample_pca.pdf to describe how the cell types relate.",
  "Use top_variable_genes_heatmap.pdf to describe gene clusters that rise or fall across the lineage.",
  "Biological expectation: HSC should be the most distinct stem-like group, CMP should be intermediate, and CFUE/Erythroblast should be closer to each other than to HSC."
)
writeLines(cluster_summary, "results/clustering/q4_rna_all_cells/README.txt")
