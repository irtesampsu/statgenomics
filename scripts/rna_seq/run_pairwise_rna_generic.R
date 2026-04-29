suppressPackageStartupMessages({
  library(clusterProfiler)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(org.Mm.eg.db)
})

source("scripts/utils/rna_utils.R")

write_rna_results <- function(res, out_dir, method) {
  if (method == "limma_voom") {
    all_file <- "all_genes_limma_voom.csv"
    de_file <- "de_genes_limma_voom.csv"
    fc_col <- "logFC"
    padj_col <- "adj.P.Val"
  } else {
    all_file <- "all_genes_deseq2.csv"
    de_file <- "de_genes_deseq2.csv"
    fc_col <- "log2FoldChange"
    padj_col <- "padj"
  }

  write.csv(res, file.path(out_dir, all_file), row.names = FALSE)
  keep <- !is.na(res[[padj_col]]) & res[[padj_col]] < 0.05 & abs(res[[fc_col]]) > 1
  write.csv(res[keep, , drop = FALSE], file.path(out_dir, de_file), row.names = FALSE)
}

write_go_results <- function(res, background_gene_ids, out_dir, method) {
  if (method == "limma_voom") {
    fc_col <- "logFC"
    padj_col <- "adj.P.Val"
    go_file <- "go_enrichment_limma_voom.csv"
  } else {
    fc_col <- "log2FoldChange"
    padj_col <- "padj"
    go_file <- "go_enrichment_deseq2.csv"
  }

  keep <- !is.na(res[[padj_col]]) & res[[padj_col]] < 0.05 & abs(res[[fc_col]]) > 1
  sig_gene_ids <- map_ensembl_to_entrez(res$gene_id[keep], org.Mm.eg.db)
  bg_entrez_ids <- map_ensembl_to_entrez(background_gene_ids, org.Mm.eg.db)

  if (length(sig_gene_ids) < 10 || length(bg_entrez_ids) < 10) {
    write.csv(data.frame(), file.path(out_dir, go_file), row.names = FALSE)
    return(invisible(NULL))
  }

  ego <- enrichGO(
    gene = sig_gene_ids,
    universe = bg_entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH"
  )
  write.csv(as.data.frame(ego), file.path(out_dir, go_file), row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  cell_a = "CMP",
  cell_b = "Erythroblast",
  method = "limma_voom",
  out_dir = NULL,
  min_count = 10L,
  min_samples = 2L
)

if (length(args) %% 2 != 0) {
  stop("Arguments must be passed as --key value pairs.")
}
if (length(args) > 0) {
  keys <- sub("^--", "", args[seq(1, length(args), by = 2)])
  vals <- args[seq(2, length(args), by = 2)]
  for (i in seq_along(keys)) {
    k <- keys[i]
    if (!k %in% names(defaults)) {
      stop(sprintf("Unknown argument: --%s", k))
    }
    defaults[[k]] <- vals[i]
  }
}

opt <- list(
  cell_a = as.character(defaults$cell_a),
  cell_b = as.character(defaults$cell_b),
  method = as.character(defaults$method),
  out_dir = if (is.null(defaults$out_dir) || identical(defaults$out_dir, "NULL")) NULL else as.character(defaults$out_dir),
  min_count = as.integer(defaults$min_count),
  min_samples = as.integer(defaults$min_samples)
)

method <- tolower(opt$method)
if (!method %in% c("limma_voom", "deseq2")) {
  stop("--method must be one of: limma_voom, deseq2")
}

if (is.null(opt$out_dir)) {
  opt$out_dir <- sprintf(
    "results/%s/pairwise_%s_vs_%s",
    method,
    tolower(opt$cell_a),
    tolower(opt$cell_b)
  )
}
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

sample_info <- get_pairwise_rna_scriptseq(opt$cell_a, opt$cell_b)
counts <- build_rna_count_matrix(sample_info)
counts <- filter_rna_counts(counts, min_count = opt$min_count, min_samples = opt$min_samples)
group <- factor(sample_info$condition, levels = c(opt$cell_a, opt$cell_b))
design <- model.matrix(~ group)

if (method == "limma_voom") {
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$gene_id <- rownames(res)
  res$comparison <- paste(opt$cell_b, "vs", opt$cell_a)
  write_rna_results(res, opt$out_dir, method)
  write_go_results(res, rownames(dge), opt$out_dir, method)
} else {
  col_data <- data.frame(condition = group)
  rownames(col_data) <- colnames(counts)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds, contrast = c("condition", opt$cell_b, opt$cell_a)))
  res$gene_id <- rownames(res)
  res$comparison <- paste(opt$cell_b, "vs", opt$cell_a)
  res <- res[order(res$padj), , drop = FALSE]
  write_rna_results(res, opt$out_dir, method)
  write_go_results(res, rownames(counts), opt$out_dir, method)
}
