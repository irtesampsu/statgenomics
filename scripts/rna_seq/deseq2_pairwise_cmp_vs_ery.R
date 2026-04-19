suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)

sample_info <- data.frame(
  sample = c("CMP_rep1", "CMP_rep2", "Erythroblast_rep1", "Erythroblast_rep2"),
  condition = c("CMP", "CMP", "Erythroblast", "Erythroblast"),
  file = c(
    "data/rna_seq/CMP_rep1_ENCFF623OLU.tsv",
    "data/rna_seq/CMP_rep2_ENCFF691MHW.tsv",
    "data/rna_seq/Erythroblast_rep1_ENCFF342WUL.tsv",
    "data/rna_seq/Erythroblast_rep2_ENCFF858JHF.tsv"
  ),
  stringsAsFactors = FALSE
)

read_expected_counts <- function(file_path, sample_name) {
  x <- read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  x <- x[, c("gene_id", "expected_count")]
  colnames(x)[2] <- sample_name
  x
}

count_tables <- mapply(
  read_expected_counts,
  sample_info$file,
  sample_info$sample,
  SIMPLIFY = FALSE
)

count_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_tables)
count_df[is.na(count_df)] <- 0

counts <- as.matrix(count_df[, -1])
rownames(counts) <- count_df$gene_id
storage.mode(counts) <- "integer"

keep <- rowSums(counts >= 10) >= 2
counts <- counts[keep, ]

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

sig_gene_ids <- unique(deg_df$gene_id)
sig_gene_ids <- as.character(sig_gene_ids[!is.na(sig_gene_ids)])
sig_gene_ids <- sub("\\..*$", "", sig_gene_ids)

id_map <- bitr(
  sig_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

sig_gene_ids <- unique(id_map$ENTREZID)

if (length(sig_gene_ids) == 0) {
  stop("GO enrichment failed: no significant DE genes could be mapped from ENSEMBL IDs to mouse ENTREZID values.")
}

ego <- enrichGO(
  gene = sig_gene_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH"
)

if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  stop("GO enrichment failed: enrichGO returned no terms after mapping valid mouse ENTREZID values.")
}

write.csv(as.data.frame(ego), "results/deseq2/go_enrichment_deseq2.csv", row.names = FALSE)

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
