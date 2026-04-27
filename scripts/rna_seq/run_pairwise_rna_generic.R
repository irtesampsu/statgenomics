suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
})

source("scripts/utils/rna_utils.R")

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
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res$gene_id <- rownames(res)
  res$comparison <- paste(opt$cell_b, "vs", opt$cell_a)
  write.csv(res, file.path(opt$out_dir, "de_genes_limma_voom.csv"), row.names = FALSE)
} else {
  col_data <- data.frame(condition = group)
  rownames(col_data) <- colnames(counts)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds, contrast = c("condition", opt$cell_b, opt$cell_a)))
  res$gene_id <- rownames(res)
  res$comparison <- paste(opt$cell_b, "vs", opt$cell_a)
  write.csv(res, file.path(opt$out_dir, "de_genes_deseq2.csv"), row.names = FALSE)
}

notes <- c(
  "Pairwise RNA generic run notes",
  "==============================",
  sprintf("Cell A: %s", opt$cell_a),
  sprintf("Cell B: %s", opt$cell_b),
  sprintf("Method: %s", method),
  sprintf("Output directory: %s", opt$out_dir),
  sprintf("Samples used: %s", paste(sample_info$sample, collapse = ", ")),
  sprintf("Genes tested after filtering: %d", nrow(counts)),
  sprintf("Filter: min_count >= %d in at least %d samples", opt$min_count, opt$min_samples),
  "",
  "This script is for rapid expansion to additional assigned pairs.",
  "Use limma_voom by default based on first-pair method comparison unless a rerun decision changes."
)
writeLines(notes, file.path(opt$out_dir, "analysis_notes.txt"))
