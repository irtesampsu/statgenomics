suppressPackageStartupMessages({
  library(ChIPseeker)
  library(clusterProfiler)
  library(GenomicRanges)
  library(ggplot2)
  library(org.Mm.eg.db)
  library(scales)
  library(stringr)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
})

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  pair_dir = "results/atac_seq/pairwise_cmp_vs_ery",
  cell_a = NULL,
  cell_b = NULL
)
if (length(args) %% 2 != 0) {
  stop("Arguments must be passed as --key value pairs.")
}
if (length(args) > 0) {
  keys <- sub("^--", "", args[seq(1, length(args), by = 2)])
  vals <- args[seq(2, length(args), by = 2)]
  for (i in seq_along(keys)) {
    defaults[[keys[i]]] <- vals[i]
  }
}

out_dir <- as.character(defaults$pair_dir)
in_path <- file.path(out_dir, "DARs_annotated.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (is.null(defaults$cell_a) || is.null(defaults$cell_b) || !nzchar(defaults$cell_a) || !nzchar(defaults$cell_b)) {
  pair_base <- basename(out_dir)
  parts <- strsplit(pair_base, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    defaults$cell_a <- sub("^pairwise_", "", parts[1])
    defaults$cell_b <- parts[2]
  } else {
    defaults$cell_a <- "first cell"
    defaults$cell_b <- "second cell"
  }
}
cell_a <- as.character(defaults$cell_a)
cell_b <- as.character(defaults$cell_b)

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  lines <- readLines(path, warn = FALSE)
  # Guard against blank/whitespace-only CSV files.
  if (length(lines) == 0 || all(!nzchar(trimws(lines)))) return(data.frame())
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) data.frame()
  )
}

dar_deseq2 <- safe_read_csv(file.path(out_dir, "DARs_annotated_deseq2.csv"))
dar_limma <- safe_read_csv(file.path(out_dir, "DARs_annotated_limma_voom.csv"))
dar <- safe_read_csv(in_path)
if (nrow(dar_deseq2) == 0 && nrow(dar_limma) == 0 && nrow(dar) == 0) quit(save = "no", status = 0)
if (nrow(dar_deseq2) == 0 && nrow(dar_limma) == 0 && nrow(dar) > 0) {
  if ("source_method" %in% colnames(dar)) {
    dar_deseq2 <- dar[tolower(as.character(dar$source_method)) == "deseq2", , drop = FALSE]
    dar_limma <- dar[tolower(as.character(dar$source_method)) == "limma_voom", , drop = FALSE]
  }
}

parse_entrez <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  out <- sub("^\\s*([0-9]+).*$", "\\1", x)
  out[!grepl("^[0-9]+$", out)] <- NA_character_
  out
}

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
    as.numeric(parts[1]) / as.numeric(parts[2])
  }, numeric(1))
}

annotate_peak_table <- function(peak_df) {
  if (nrow(peak_df) == 0 || !all(c("chr", "start", "end") %in% colnames(peak_df))) {
    return(data.frame())
  }
  gr <- GRanges(
    seqnames = peak_df$chr,
    ranges = IRanges(start = peak_df$start, end = peak_df$end)
  )
  anno <- as.data.frame(annotatePeak(gr, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene))
  anno$peak_key <- paste(anno$seqnames, anno$start, anno$end, sep = ":")
  anno
}

is_proximal_annotation <- function(annotation_vec) {
  ann <- tolower(annotation_vec)
  !grepl("distal intergenic", ann, fixed = TRUE)
}

run_directional_go <- function(entrez_ids, universe_ids, label, method_label) {
  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)])
  universe_ids <- unique(universe_ids[!is.na(universe_ids) & nzchar(universe_ids)])
  if (length(entrez_ids) < 10 || length(universe_ids) < 10) {
    return(data.frame())
  }
  ego <- enrichGO(gene = entrez_ids, universe = universe_ids, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH")
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(data.frame())
  ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  out <- as.data.frame(ego)
  if (nrow(out) == 0) return(out)
  out$direction <- label
  out$method <- method_label
  out$gene_ratio <- parse_ratio(out$GeneRatio)
  out$bg_ratio <- parse_ratio(out$BgRatio)
  out$fold_enrichment <- out$gene_ratio / out$bg_ratio
  out
}

all_peak_anno_path <- file.path(out_dir, "all_tested_peaks_annotated.csv")
if (file.exists(all_peak_anno_path)) {
  all_peak_anno <- safe_read_csv(all_peak_anno_path)
} else {
  all_peak_path <- file.path(out_dir, "all_peaks_deseq2.csv")
  all_peak_anno <- if (file.exists(all_peak_path)) annotate_peak_table(safe_read_csv(all_peak_path)) else data.frame()
  write.csv(all_peak_anno, all_peak_anno_path, row.names = FALSE)
}

universe_entrez <- if (nrow(all_peak_anno) > 0 && "geneId" %in% colnames(all_peak_anno)) {
  parse_entrez(all_peak_anno$geneId)
} else {
  parse_entrez(c(dar_deseq2$geneId, dar_limma$geneId, dar$geneId))
}

run_go_for_method <- function(dar_df, method_label) {
  if (nrow(dar_df) == 0) return(data.frame())
  open_in_second <- parse_entrez(dar_df$geneId[dar_df$log2FoldChange > 0 & dar_df$padj < 0.05])
  open_in_first <- parse_entrez(dar_df$geneId[dar_df$log2FoldChange < 0 & dar_df$padj < 0.05])
  rbind(
    run_directional_go(open_in_first, universe_entrez, paste("More open in", cell_a), method_label),
    run_directional_go(open_in_second, universe_entrez, paste("More open in", cell_b), method_label)
  )
}

go_combined <- rbind(
  run_go_for_method(dar_deseq2, "DESeq2"),
  run_go_for_method(dar_limma, "limma-voom")
)
write.csv(go_combined, file.path(out_dir, "atac_go_simplified_combined.csv"), row.names = FALSE)

plot_go_bar <- function(go_df, title_text, out_path) {
  if (nrow(go_df) == 0) return(invisible(NULL))
  go_df <- go_df[is.finite(go_df$p.adjust) & is.finite(go_df$fold_enrichment), , drop = FALSE]
  if (nrow(go_df) == 0) return(invisible(NULL))
  go_df <- do.call(rbind, lapply(split(go_df, list(go_df$method, go_df$direction), drop = TRUE), function(x) {
    top <- head(x[order(x$p.adjust), , drop = FALSE], 10)
    top
  }))
  go_df$signed_fold <- ifelse(go_df$direction == paste("More open in", cell_a), -abs(go_df$fold_enrichment), abs(go_df$fold_enrichment))
  go_df <- go_df[order(go_df$method, go_df$signed_fold), , drop = FALSE]
  go_df$term_label <- str_wrap(go_df$Description, width = 45)
  go_df$term_key <- paste0(go_df$term_label, "___", seq_len(nrow(go_df)))
  go_df$term_key <- factor(go_df$term_key, levels = unique(go_df$term_key))
  p <- ggplot(go_df, aes(x = signed_fold, y = term_key, fill = direction)) +
    geom_col(width = 0.75, alpha = 0.9) +
    facet_wrap(~ method, scales = "free_y") +
    scale_fill_manual(values = c(setNames("#2C7BB6", paste("More open in", cell_a)), setNames("#D7191C", paste("More open in", cell_b)))) +
    scale_y_discrete(labels = function(x) sub("___[0-9]+$", "", x)) +
    labs(
      title = title_text,
      x = "Signed fold enrichment (absolute magnitude)",
      y = "Biological process",
      fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )
  ggsave(out_path, p, width = 12, height = 7)
}

plot_go_bar(go_combined, "ATAC GO Enrichment (top directional terms)", file.path(out_dir, "atac_go_directional_barplot.pdf"))
