suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
})

dir.create("results/clustering/q4_rna_all_cells", recursive = TRUE, showWarnings = FALSE)

sample_info <- data.frame(
  sample = c(
    "HSC_rep1", "HSC_rep2",
    "CMP_rep1", "CMP_rep2",
    "CFUE_rep1", "CFUE_rep2",
    "Erythroblast_rep1", "Erythroblast_rep2"
  ),
  condition = c(
    "HSC", "HSC",
    "CMP", "CMP",
    "CFUE", "CFUE",
    "Erythroblast", "Erythroblast"
  ),
  file = c(
    "data/rna_seq/HSC_rep1_ENCFF247FEJ.tsv",
    "data/rna_seq/HSC_rep2_ENCFF064MKY.tsv",
    "data/rna_seq/CMP_rep1_ENCFF623OLU.tsv",
    "data/rna_seq/CMP_rep2_ENCFF691MHW.tsv",
    "data/rna_seq/CFUE_rep1_ENCFF667IDY.tsv",
    "data/rna_seq/CFUE_rep2_ENCFF655LMK.tsv",
    "data/rna_seq/Erythroblast_rep1_ENCFF342WUL.tsv",
    "data/rna_seq/Erythroblast_rep2_ENCFF858JHF.tsv"
  ),
  stringsAsFactors = FALSE
)

missing_files <- sample_info$file[!file.exists(sample_info$file)]
if (length(missing_files) > 0) {
  stop(
    paste(
      "Q4 RNA clustering is missing required ScriptSeq files:",
      paste(missing_files, collapse = ", ")
    )
  )
}

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

sample_info$condition <- factor(
  sample_info$condition,
  levels = c("HSC", "CMP", "CFUE", "Erythroblast")
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info[, c("sample", "condition")],
  design = ~ condition
)
rownames(colData(dds)) <- sample_info$sample

vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

sample_cor <- cor(vsd_mat)
write.csv(sample_cor, "results/clustering/q4_rna_all_cells/sample_correlation.csv")

pdf("results/clustering/q4_rna_all_cells/sample_correlation_heatmap.pdf", width = 7, height = 6)
pheatmap(
  sample_cor,
  annotation_col = data.frame(condition = sample_info$condition, row.names = sample_info$sample),
  annotation_row = data.frame(condition = sample_info$condition, row.names = sample_info$sample),
  main = "RNA-seq Sample Correlation"
)
dev.off()

sample_dist <- dist(t(vsd_mat))
hc <- hclust(sample_dist)

pdf("results/clustering/q4_rna_all_cells/sample_hclust.pdf", width = 8, height = 6)
plot(hc, main = "Q4 RNA-seq Hierarchical Clustering")
dev.off()

pdf("results/clustering/q4_rna_all_cells/sample_pca.pdf", width = 7, height = 6)
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

gene_var <- apply(vsd_mat, 1, var)
top_gene_ids <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(100, length(gene_var)))]
top_gene_mat <- vsd_mat[top_gene_ids, , drop = FALSE]

pdf("results/clustering/q4_rna_all_cells/top_variable_genes_heatmap.pdf", width = 8, height = 10)
pheatmap(
  top_gene_mat,
  scale = "row",
  show_rownames = FALSE,
  annotation_col = data.frame(condition = sample_info$condition, row.names = sample_info$sample),
  main = "Top Variable RNA-seq Genes"
)
dev.off()

cluster_summary <- c(
  "Question 4 guidance",
  "Use sample_hclust.pdf and sample_pca.pdf to describe how the cell types relate.",
  "Use top_variable_genes_heatmap.pdf to describe gene clusters that rise or fall across the lineage.",
  "Biological expectation: HSC should be the most distinct stem-like group, CMP should be intermediate, and CFUE/Erythroblast should be closer to each other than to HSC."
)
writeLines(cluster_summary, "results/clustering/q4_rna_all_cells/README.txt")
