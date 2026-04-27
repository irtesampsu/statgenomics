suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

source("scripts/utils/rna_utils.R")

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

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

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Erythroblast", "CMP"))
res <- lfcShrink(dds, coef = "condition_Erythroblast_vs_CMP", res = res, type = "apeglm")

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, "results/deseq2/all_genes_deseq2.csv", row.names = FALSE)

deg_df <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg_df, "results/deseq2/de_genes_deseq2.csv", row.names = FALSE)

sig_gene_ids <- map_ensembl_to_entrez(deg_df$gene_id, org.Mm.eg.db)
bg_entrez_ids <- map_ensembl_to_entrez(rownames(counts), org.Mm.eg.db)

if (length(sig_gene_ids) == 0) {
  stop("GO enrichment failed: no significant DE genes could be mapped from ENSEMBL IDs to mouse ENTREZID values.")
}

ego <- enrichGO(
  gene = sig_gene_ids,
  universe = bg_entrez_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH"
)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  stop("GO enrichment failed: enrichGO returned no terms after mapping valid mouse ENTREZID values.")
}

write.csv(as.data.frame(ego), "results/deseq2/go_enrichment_deseq2.csv", row.names = FALSE)
writeLines(
  c(
    "RNA pairwise analysis inputs",
    "Assay: ScriptSeq only",
    paste("Samples:", paste(sample_info$sample, collapse = ", ")),
    sprintf("Genes tested after filtering: %d", nrow(counts)),
    sprintf("Significant DE genes: %d", nrow(deg_df))
  ),
  "results/deseq2/analysis_notes.txt"
)

vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

pdf("results/deseq2/pca_deseq2.pdf")
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

sample_cor <- cor(vsd_mat)
write.csv(sample_cor, "results/deseq2/sample_correlation_deseq2.csv")

pdf("results/deseq2/sample_heatmap_deseq2.pdf")
pheatmap(sample_cor, main = "Sample Correlation")
dev.off()

sample_dist <- dist(t(vsd_mat))
pdf("results/deseq2/hclust_deseq2.pdf")
plot(hclust(sample_dist), main = "Hierarchical Clustering of Samples")
dev.off()

plot_top_variable_genes_heatmap(
  expr_mat = vsd_mat,
  sample_conditions = sample_info$condition,
  sample_names = sample_info$sample,
  out_file = "results/deseq2/top1000_variable_genes_heatmap_deseq2.pdf",
  title_text = "Top 1000 Variable Genes - Hierarchical Clustering\nCMP vs Erythroblast",
  top_n = 1000L,
  annotation_palette = c(CMP = "#67A9CF", Erythroblast = "#EF8A62")
)

top_n <- min(50, nrow(deg_df))
if (top_n > 1) {
  top_gene_ids <- deg_df$gene_id[seq_len(top_n)]
  heatmap_mat <- vsd_mat[rownames(vsd_mat) %in% top_gene_ids, , drop = FALSE]
  heatmap_mat <- heatmap_mat[match(top_gene_ids, rownames(heatmap_mat)), , drop = FALSE]

  pdf("results/deseq2/top_de_genes_heatmap_deseq2.pdf", width = 7, height = 9)
  pheatmap(
    heatmap_mat,
    scale = "row",
    annotation_col = data.frame(condition = sample_info$condition, row.names = sample_info$sample),
    show_rownames = FALSE,
    main = "Top Differentially Expressed Genes"
  )
  dev.off()
}

volcano_df <- res_df[!is.na(res_df$padj), ]
volcano_df$significant <- volcano_df$padj < 0.05 & abs(volcano_df$log2FoldChange) > 1

p <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(size = 1) +
  theme_minimal() +
  labs(
    title = "DESeq2 Volcano Plot",
    x = "Log2 Fold Change (Erythroblast vs CMP)",
    y = "-log10 adjusted p-value"
  )

ggsave("results/deseq2/volcano_deseq2.pdf", plot = p, width = 6, height = 5)
