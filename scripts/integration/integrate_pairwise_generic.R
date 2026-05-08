suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  cell_a = "CMP",
  cell_b = "Erythroblast",
  deseq_path = "results/deseq2/de_genes_deseq2.csv",
  limma_path = "results/limma_voom/de_genes_limma_voom.csv",
  atac_annot_path = "results/atac_seq/pairwise_cmp_vs_ery/DARs_annotated.csv",
  out_dir = "results/integration/pairwise_cmp_vs_ery"
)
if (length(args) %% 2 != 0) {
  stop("Arguments must be passed as --key value pairs.")
}
if (length(args) > 0) {
  keys <- sub("^--", "", args[seq(1, length(args), by = 2)])
  vals <- args[seq(2, length(args), by = 2)]
  for (i in seq_along(keys)) defaults[[keys[i]]] <- vals[i]
}

clean_ensembl <- function(ids) sub("\\..*$", "", ids)
read_results <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

map_ensembl_to_symbol <- function(ids) {
  ids <- unique(clean_ensembl(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) return(character(0))
  unname(suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))[ids])
}
map_entrez_to_ensembl <- function(ids) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) return(character(0))
  unname(suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = ids, column = "ENSEMBL", keytype = "ENTREZID", multiVals = "first"))[ids])
}

prepare_rna_method <- function(df, fc_col, padj_col) {
  out <- df[, c("gene_id", fc_col, padj_col)]
  colnames(out) <- c("gene_id", "fold_change", "padj")
  out$gene_id <- clean_ensembl(out$gene_id)
  out$gene_symbol <- map_ensembl_to_symbol(out$gene_id)
  out
}

prepare_atac_gene_table <- function(atac_df) {
  raw_gene_id <- as.character(atac_df$geneId); raw_gene_id[is.na(raw_gene_id)] <- ""
  parsed_gene_id <- sub("^\\s*([0-9]+).*$", "\\1", raw_gene_id)
  parsed_gene_id[!grepl("^[0-9]+$", parsed_gene_id)] <- NA_character_
  atac_df$entrez_id <- parsed_gene_id
  atac_df <- atac_df[!is.na(atac_df$entrez_id), , drop = FALSE]
  if (nrow(atac_df) == 0) return(data.frame())
  split_by_gene <- split(atac_df, atac_df$entrez_id)
  out <- do.call(rbind, lapply(names(split_by_gene), function(gene_id) {
    sub_df <- split_by_gene[[gene_id]]
    bp <- suppressWarnings(min(sub_df$padj, na.rm = TRUE)); if (!is.finite(bp)) bp <- NA_real_
    data.frame(entrez_id = gene_id, mean_log2FoldChange = mean(sub_df$log2FoldChange, na.rm = TRUE), best_padj = bp, representative_annotation = as.character(sub_df$annotation[1]), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out$ensembl_id <- map_entrez_to_ensembl(out$entrez_id)
  out
}

dir.create(defaults$out_dir, recursive = TRUE, showWarnings = FALSE)
deseq_tbl <- prepare_rna_method(read_results(defaults$deseq_path), "log2FoldChange", "padj")
limma_tbl <- prepare_rna_method(read_results(defaults$limma_path), "logFC", "adj.P.Val")

rna_shared <- merge(
  deseq_tbl[, c("gene_id", "fold_change")],
  limma_tbl[, c("gene_id", "fold_change")],
  by = "gene_id",
  suffixes = c("_deseq2", "_limma")
)
rna_shared$mean_rna_log2fc <- rowMeans(rna_shared[, c("fold_change_deseq2", "fold_change_limma")], na.rm = TRUE)

make_concordance_plot <- function(atac_annot_path, out_file, title_suffix = NULL) {
  if (!file.exists(atac_annot_path)) return(invisible(NULL))
  atac_tbl <- prepare_atac_gene_table(read_results(atac_annot_path))
  if (nrow(atac_tbl) == 0) return(invisible(NULL))
  atac_tbl$ensembl_id <- clean_ensembl(atac_tbl$ensembl_id)

  quad <- merge(
    rna_shared[, c("gene_id", "mean_rna_log2fc")],
    atac_tbl[, c("ensembl_id", "mean_log2FoldChange", "representative_annotation")],
    by.x = "gene_id",
    by.y = "ensembl_id"
  )
  if (nrow(quad) == 0) return(invisible(NULL))

  quad$concordance <- ifelse(sign(quad$mean_rna_log2fc) == sign(quad$mean_log2FoldChange), "Concordant", "Discordant")
  quad$quadrant <- ifelse(
    quad$mean_rna_log2fc >= 0 & quad$mean_log2FoldChange >= 0, "Q1: RNA up, ATAC up",
    ifelse(
      quad$mean_rna_log2fc < 0 & quad$mean_log2FoldChange >= 0, "Q2: RNA down, ATAC up",
      ifelse(
        quad$mean_rna_log2fc < 0 & quad$mean_log2FoldChange < 0, "Q3: RNA down, ATAC down",
        "Q4: RNA up, ATAC down"
      )
    )
  )

  top_gene_ids <- head(quad$gene_id[order(-(abs(quad$mean_rna_log2fc) + abs(quad$mean_log2FoldChange)))], 20)
  top_gene_ids <- unique(top_gene_ids[!is.na(top_gene_ids) & nzchar(top_gene_ids)])
  top_symbols <- map_ensembl_to_symbol(top_gene_ids)
  symbol_map <- setNames(top_symbols, top_gene_ids)
  quad$gene_label <- ifelse(
    quad$gene_id %in% names(symbol_map) & !is.na(symbol_map[quad$gene_id]) & nzchar(symbol_map[quad$gene_id]),
    unname(symbol_map[quad$gene_id]),
    quad$gene_id
  )
  quad$label_flag <- quad$gene_id %in% top_gene_ids

  q_counts <- as.data.frame(table(quad$quadrant), stringsAsFactors = FALSE)
  colnames(q_counts) <- c("quadrant", "n")
  q_labels <- data.frame(
    quadrant = c("Q1: RNA up, ATAC up", "Q2: RNA down, ATAC up", "Q3: RNA down, ATAC down", "Q4: RNA up, ATAC down"),
    x = c(Inf, -Inf, -Inf, Inf),
    y = c(Inf, Inf, -Inf, -Inf),
    hjust = c(1.05, -0.05, -0.05, 1.05),
    vjust = c(1.4, 1.4, -0.4, -0.4),
    stringsAsFactors = FALSE
  )
  q_labels <- merge(q_labels, q_counts, by = "quadrant", all.x = TRUE)
  q_labels$n[is.na(q_labels$n)] <- 0

  title_text <- sprintf("RNA-ATAC Direction Concordance: %s vs %s", defaults$cell_b, defaults$cell_a)
  if (!is.null(title_suffix) && nzchar(title_suffix)) {
    title_text <- sprintf("%s (%s)", title_text, title_suffix)
  }

  p <- ggplot(quad, aes(mean_rna_log2fc, mean_log2FoldChange, color = concordance)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey60") +
    geom_point(alpha = 0.6, size = 1) +
    geom_text(
      data = quad[quad$label_flag, , drop = FALSE],
      aes(label = gene_label),
      size = 2.4,
      vjust = -0.5,
      show.legend = FALSE
    ) +
    geom_text(
      data = q_labels,
      aes(x = x, y = y, label = paste0(quadrant, "\nN=", n)),
      inherit.aes = FALSE,
      size = 3.0,
      hjust = q_labels$hjust,
      vjust = q_labels$vjust,
      color = "grey20"
    ) +
    theme_minimal(base_size = 11) +
    labs(title = title_text, x = "Mean RNA log2FC", y = "ATAC mean log2FC")

  ggsave(file.path(defaults$out_dir, out_file), plot = p, width = 8, height = 6)
}

# Consensus (existing output path retained)
make_concordance_plot(
  defaults$atac_annot_path,
  "rna_atac_quadrant_plot.pdf",
  "consensus ATAC DARs"
)

# Optional method-specific concordance plots if ATAC per-method annotations exist.
atac_dir <- dirname(defaults$atac_annot_path)
make_concordance_plot(
  file.path(atac_dir, "DARs_annotated_deseq2.csv"),
  "rna_atac_quadrant_plot_deseq2.pdf",
  "ATAC DESeq2 DARs"
)
make_concordance_plot(
  file.path(atac_dir, "DARs_annotated_limma_voom.csv"),
  "rna_atac_quadrant_plot_limma_voom.pdf",
  "ATAC limma-voom DARs"
)
