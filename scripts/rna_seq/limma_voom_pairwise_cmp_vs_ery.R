suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

dir.create("results/limma_voom", recursive = TRUE, showWarnings = FALSE)

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

group <- factor(sample_info$condition, levels = c("CMP", "Erythroblast"))

y <- DGEList(counts = counts, group = group)
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~ group)
v <- voom(y, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

res_df <- topTable(
  fit,
  coef = "groupErythroblast",
  number = Inf,
  sort.by = "P"
)
res_df$gene_id <- rownames(res_df)
write.csv(res_df, "results/limma_voom/all_genes_limma_voom.csv", row.names = FALSE)

deg_df <- subset(res_df, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(deg_df, "results/limma_voom/de_genes_limma_voom.csv", row.names = FALSE)

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

write.csv(as.data.frame(ego), "results/limma_voom/go_enrichment_limma_voom.csv", row.names = FALSE)

log_cpm <- cpm(y, log = TRUE, prior.count = 2)
sample_cor <- cor(log_cpm)
write.csv(sample_cor, "results/limma_voom/sample_correlation_limma_voom.csv")

pdf("results/limma_voom/sample_heatmap_limma_voom.pdf")
pheatmap(sample_cor, main = "Sample Correlation")
dev.off()

sample_dist <- dist(t(log_cpm))
pdf("results/limma_voom/hclust_limma_voom.pdf")
plot(hclust(sample_dist), main = "Hierarchical Clustering of Samples")
dev.off()

pca_df <- as.data.frame(prcomp(t(log_cpm), scale. = FALSE)$x[, 1:2])
pca_df$sample <- rownames(pca_df)
pca_df$condition <- group

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7) +
  theme_minimal() +
  labs(title = "limma-voom PCA")

ggsave("results/limma_voom/pca_limma_voom.pdf", plot = p, width = 6, height = 5)

top_n <- min(50, nrow(deg_df))
if (top_n > 1) {
  top_gene_ids <- deg_df$gene_id[seq_len(top_n)]
  heatmap_mat <- log_cpm[rownames(log_cpm) %in% top_gene_ids, , drop = FALSE]
  heatmap_mat <- heatmap_mat[match(top_gene_ids, rownames(heatmap_mat)), , drop = FALSE]

  pdf("results/limma_voom/top_de_genes_heatmap_limma_voom.pdf", width = 7, height = 9)
  pheatmap(
    heatmap_mat,
    scale = "row",
    annotation_col = data.frame(condition = group, row.names = sample_info$sample),
    show_rownames = FALSE,
    main = "Top Differentially Expressed Genes"
  )
  dev.off()
}

volcano_df <- res_df[!is.na(res_df$adj.P.Val), ]
volcano_df$significant <- volcano_df$adj.P.Val < 0.05 & abs(volcano_df$logFC) > 1

p <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(size = 1) +
  theme_minimal() +
  labs(
    title = "limma-voom Volcano Plot",
    x = "Log2 Fold Change (Erythroblast vs CMP)",
    y = "-log10 adjusted p-value"
  )

ggsave("results/limma_voom/volcano_limma_voom.pdf", plot = p, width = 6, height = 5)
