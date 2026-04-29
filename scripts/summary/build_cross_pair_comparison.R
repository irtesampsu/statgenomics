suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(pheatmap)
})

dir.create("results/summary", recursive = TRUE, showWarnings = FALSE)

clean_ensembl <- function(ids) sub("\\..*$", "", ids)

parse_entrez <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  out <- sub("^\\s*([0-9]+).*$", "\\1", x)
  out[!grepl("^[0-9]+$", out)] <- NA_character_
  out
}

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(data.frame())
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

jaccard <- function(a, b) {
  a <- unique(a[!is.na(a) & nzchar(a)])
  b <- unique(b[!is.na(b) & nzchar(b)])
  if (length(union(a, b)) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

run_go_from_ensembl <- function(ensembl_ids) {
  ensembl_ids <- unique(clean_ensembl(ensembl_ids))
  ensembl_ids <- ensembl_ids[!is.na(ensembl_ids) & nzchar(ensembl_ids)]
  if (length(ensembl_ids) < 10) return(data.frame())
  entrez_ids <- suppressMessages(mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"))
  entrez_ids <- unique(unname(entrez_ids))
  entrez_ids <- entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)]
  if (length(entrez_ids) < 10) return(data.frame())
  as.data.frame(enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH"))
}

pair_specs <- list(
  cmp_vs_ery = list(
    deseq = "results/deseq2/de_genes_deseq2.csv",
    limma = "results/limma_voom/de_genes_limma_voom.csv",
    atac = "results/atac_seq/pairwise_cmp_vs_ery/DARs_annotated.csv"
  ),
  cfue_vs_ery = list(
    deseq = "results/deseq2/pairwise_cfue_vs_ery/de_genes_deseq2.csv",
    limma = "results/limma_voom/pairwise_cfue_vs_ery/de_genes_limma_voom.csv",
    atac = "results/atac_seq/pairwise_cfue_vs_ery/DARs_annotated.csv"
  )
)

de_sets <- list()
dar_gene_sets <- list()
go_by_pair <- list()
comparison_rows <- list()

for (pair_name in names(pair_specs)) {
  spec <- pair_specs[[pair_name]]
  deseq <- read_csv_safe(spec$deseq)
  limma <- read_csv_safe(spec$limma)
  atac <- read_csv_safe(spec$atac)

  deseq_sig <- clean_ensembl(deseq$gene_id[!is.na(deseq$padj) & deseq$padj < 0.05 & abs(deseq$log2FoldChange) > 1])
  limma_sig <- clean_ensembl(limma$gene_id[!is.na(limma$adj.P.Val) & limma$adj.P.Val < 0.05 & abs(limma$logFC) > 1])
  shared_de <- intersect(deseq_sig, limma_sig)

  atac_entrez <- parse_entrez(atac$geneId[!is.na(atac$padj) & atac$padj < 0.05 & abs(atac$log2FoldChange) > 1])
  atac_entrez <- unique(atac_entrez[!is.na(atac_entrez)])
  de_sets[[pair_name]] <- shared_de
  dar_gene_sets[[pair_name]] <- atac_entrez

  go_df <- run_go_from_ensembl(shared_de)
  go_by_pair[[pair_name]] <- go_df
  write.csv(go_df, file.path("results/summary", paste0("go_shared_rna_", pair_name, ".csv")), row.names = FALSE)

  comparison_rows[[pair_name]] <- data.frame(
    pair = pair_name,
    deseq_sig = length(unique(deseq_sig)),
    limma_sig = length(unique(limma_sig)),
    shared_de = length(unique(shared_de)),
    atac_dar_genes = length(unique(atac_entrez)),
    stringsAsFactors = FALSE
  )
}

comparison_df <- do.call(rbind, comparison_rows)
write.csv(comparison_df, "results/summary/cross_pair_counts_summary.csv", row.names = FALSE)

pair_names <- names(pair_specs)
de_j <- matrix(NA_real_, nrow = length(pair_names), ncol = length(pair_names), dimnames = list(pair_names, pair_names))
dar_j <- de_j
for (i in seq_along(pair_names)) {
  for (j in seq_along(pair_names)) {
    de_j[i, j] <- jaccard(de_sets[[pair_names[i]]], de_sets[[pair_names[j]]])
    dar_j[i, j] <- jaccard(dar_gene_sets[[pair_names[i]]], dar_gene_sets[[pair_names[j]]])
  }
}
write.csv(as.data.frame(de_j), "results/summary/jaccard_shared_de_genes.csv")
write.csv(as.data.frame(dar_j), "results/summary/jaccard_dar_nearby_genes.csv")

pdf("results/summary/jaccard_shared_de_genes_heatmap.pdf", width = 5.5, height = 4.5)
pheatmap(de_j, cluster_rows = FALSE, cluster_cols = FALSE, main = "Jaccard: shared DE genes")
dev.off()

pdf("results/summary/jaccard_dar_nearby_genes_heatmap.pdf", width = 5.5, height = 4.5)
pheatmap(dar_j, cluster_rows = FALSE, cluster_cols = FALSE, main = "Jaccard: DAR-nearby genes")
dev.off()

top_terms <- unique(unlist(lapply(go_by_pair, function(df) head(df$Description, 15))))
top_terms <- top_terms[!is.na(top_terms) & nzchar(top_terms)]
if (length(top_terms) > 0) {
  go_mat <- matrix(0, nrow = length(top_terms), ncol = length(pair_names), dimnames = list(top_terms, pair_names))
  for (pair_name in pair_names) {
    df <- go_by_pair[[pair_name]]
    if (nrow(df) == 0) next
    idx <- match(df$Description, top_terms)
    valid <- which(!is.na(idx))
    go_mat[cbind(idx[valid], rep(which(pair_names == pair_name), length(valid)))] <- -log10(pmax(df$p.adjust[valid], 1e-300))
  }
  write.csv(as.data.frame(go_mat), "results/summary/cross_pair_go_neglog10fdr_matrix.csv")
  pdf("results/summary/cross_pair_go_heatmap.pdf", width = 8, height = 10)
  pheatmap(go_mat, main = "Cross-pair RNA GO enrichment (-log10 FDR)")
  dev.off()
}
