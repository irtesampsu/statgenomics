source("scripts/utils/project_metadata.R")

is_valid_ensembl_gene_id <- function(ids) {
  grepl("^ENSMUSG[0-9]+(?:\\.[0-9]+)?$", ids)
}

clean_ensembl_ids <- function(ids) {
  sub("\\..*$", "", ids)
}

read_expected_counts <- function(file_path, sample_name) {
  x <- read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  x <- x[, c("gene_id", "expected_count")]
  x <- x[is_valid_ensembl_gene_id(x$gene_id), , drop = FALSE]
  colnames(x)[2] <- sample_name
  x
}

build_rna_count_matrix <- function(sample_info) {
  count_tables <- mapply(
    read_expected_counts,
    sample_info$file,
    sample_info$sample,
    SIMPLIFY = FALSE
  )

  count_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_tables)
  count_df[is.na(count_df)] <- 0

  counts <- as.matrix(count_df[, -1, drop = FALSE])
  rownames(counts) <- count_df$gene_id
  storage.mode(counts) <- "integer"
  counts
}

filter_rna_counts <- function(counts, min_count = 10L, min_samples = 2L) {
  keep <- rowSums(counts >= min_count) >= min_samples
  counts[keep, , drop = FALSE]
}

map_ensembl_to_entrez <- function(ids, org_db) {
  ids <- unique(clean_ensembl_ids(as.character(ids)))
  ids <- ids[is_valid_ensembl_gene_id(ids)]
  if (length(ids) == 0) {
    return(character(0))
  }

  mapped <- withCallingHandlers(
    suppressMessages(
      AnnotationDbi::mapIds(
        org_db,
        keys = ids,
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    ),
    warning = function(w) {
      # Benign warning from AnnotationDbi when one ENSEMBL ID has multiple
      # possible mappings; we intentionally keep the first mapping.
      if (grepl("1:many mapping between keys and columns", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  unique(unname(mapped[!is.na(mapped)]))
}

plot_top_variable_genes_heatmap <- function(
  expr_mat,
  sample_conditions,
  sample_names,
  out_file,
  title_text,
  top_n = 1000L,
  annotation_palette = NULL
) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("pheatmap is required for plot_top_variable_genes_heatmap().")
  }

  if (is.null(rownames(expr_mat)) || is.null(colnames(expr_mat))) {
    stop("expr_mat must have rownames (genes) and colnames (samples).")
  }
  if (length(sample_conditions) != ncol(expr_mat)) {
    stop("sample_conditions length must match number of columns in expr_mat.")
  }
  if (length(sample_names) != ncol(expr_mat)) {
    stop("sample_names length must match number of columns in expr_mat.")
  }

  gene_var <- apply(expr_mat, 1, var)
  top_n <- min(as.integer(top_n), length(gene_var))
  top_gene_ids <- names(sort(gene_var, decreasing = TRUE))[seq_len(top_n)]
  top_gene_mat <- expr_mat[top_gene_ids, , drop = FALSE]
  top_gene_scaled <- t(scale(t(top_gene_mat)))

  annotation_col <- data.frame(
    condition = sample_conditions,
    row.names = sample_names
  )

  ann_colors <- NULL
  if (!is.null(annotation_palette)) {
    ann_colors <- list(condition = annotation_palette)
  }

  grDevices::pdf(out_file, width = 8, height = 10)
  pheatmap::pheatmap(
    top_gene_scaled,
    scale = "none",
    show_rownames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    clustering_method = "complete",
    color = grDevices::colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
    main = title_text
  )
  grDevices::dev.off()
}

