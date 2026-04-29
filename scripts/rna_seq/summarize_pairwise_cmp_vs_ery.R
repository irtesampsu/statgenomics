suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(limma)
  library(org.Mm.eg.db)
  library(scales)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  deseq_path = "results/deseq2/de_genes_deseq2.csv",
  limma_path = "results/limma_voom/de_genes_limma_voom.csv",
  deseq_all_path = "results/deseq2/all_genes_deseq2.csv",
  limma_all_path = "results/limma_voom/all_genes_limma_voom.csv",
  deseq_go = "results/deseq2/go_enrichment_deseq2.csv",
  limma_go = "results/limma_voom/go_enrichment_limma_voom.csv",
  out_dir = "results/summary",
  label = "RNA",
  cell_a = "CMP",
  cell_b = "Erythroblast"
)
if (length(args) %% 2 != 0) stop("Arguments must be passed as --key value pairs.")
if (length(args) > 0) {
  keys <- sub("^--", "", args[seq(1, length(args), by = 2)])
  vals <- args[seq(2, length(args), by = 2)]
  for (i in seq_along(keys)) defaults[[keys[i]]] <- vals[i]
}

dir.create(as.character(defaults$out_dir), recursive = TRUE, showWarnings = FALSE)
cell_a <- as.character(defaults$cell_a)
cell_b <- as.character(defaults$cell_b)

clean_ensembl <- function(ids) sub("\\..*$", "", ids)
read_results <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

map_symbols <- function(ids) {
  clean_ids <- clean_ensembl(ids)
  mapped <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = clean_ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
  out <- unname(mapped)
  ifelse(is.na(out) | !nzchar(out), clean_ids, out)
}

prepare_method <- function(df, fc_col, padj_col, method_label) {
  keep <- !is.na(df[[fc_col]]) & !is.na(df[[padj_col]])
  out <- data.frame(
    gene_id = clean_ensembl(df$gene_id[keep]),
    log2fc = df[[fc_col]][keep],
    padj = df[[padj_col]][keep],
    method = method_label,
    stringsAsFactors = FALSE
  )
  out$gene <- map_symbols(out$gene_id)
  out
}

label_top_genes <- function(df, n = 10) {
  sig <- df[df$padj < 0.05 & abs(df$log2fc) >= 1, , drop = FALSE]
  up <- head(sig[sig$log2fc > 0, ][order(sig[sig$log2fc > 0, ]$padj), ], n)
  down <- head(sig[sig$log2fc < 0, ][order(sig[sig$log2fc < 0, ]$padj), ], n)
  rbind(up, down)
}

significant_genes <- function(df) {
  unique(clean_ensembl(df$gene_id[!is.na(df$padj) & df$padj < 0.05 & abs(df$log2fc) > 1]))
}

map_ensembl_to_entrez <- function(ids) {
  ids <- unique(clean_ensembl(as.character(ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) return(character(0))
  mapped <- suppressMessages(AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = ids,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  ))
  unique(unname(mapped[!is.na(mapped)]))
}

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
    as.numeric(parts[1]) / as.numeric(parts[2])
  }, numeric(1))
}

simplified_go_from_rna <- function(all_df, fc_col, padj_col, method_label, cell_a_label, cell_b_label) {
  sig <- all_df[!is.na(all_df[[padj_col]]) & all_df[[padj_col]] < 0.05 & abs(all_df[[fc_col]]) > 1, , drop = FALSE]
  universe <- map_ensembl_to_entrez(all_df$gene_id)
  run_dir <- function(gene_ids, direction_label) {
    genes <- map_ensembl_to_entrez(gene_ids)
    if (length(genes) < 10 || length(universe) < 10) return(data.frame())
    ego <- enrichGO(
      gene = genes,
      universe = universe,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH"
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(data.frame())
    ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
    out <- as.data.frame(ego)
    if (nrow(out) == 0) return(out)
    out$method <- method_label
    out$direction <- direction_label
    out$gene_ratio <- parse_ratio(out$GeneRatio)
    out$bg_ratio <- parse_ratio(out$BgRatio)
    out$fold_enrichment <- out$gene_ratio / out$bg_ratio
    out
  }
  up_b <- run_dir(sig$gene_id[sig[[fc_col]] > 0], paste("Up in", cell_b_label))
  up_a <- run_dir(sig$gene_id[sig[[fc_col]] < 0], paste("Up in", cell_a_label))
  rbind(up_a, up_b)
}

deseq <- read_results(as.character(defaults$deseq_path))
limma_df <- read_results(as.character(defaults$limma_path))
deseq_all <- read_results(as.character(defaults$deseq_all_path))
limma_all <- read_results(as.character(defaults$limma_all_path))

deseq_tbl <- prepare_method(deseq_all, "log2FoldChange", "padj", "DESeq2")
limma_tbl <- prepare_method(limma_all, "logFC", "adj.P.Val", "limma-voom")
volcano_tbl <- rbind(deseq_tbl, limma_tbl)
volcano_tbl$significant <- volcano_tbl$padj < 0.05 & abs(volcano_tbl$log2fc) >= 1
volcano_tbl$direction <- ifelse(
  volcano_tbl$significant & volcano_tbl$log2fc > 0,
  paste("Up in", cell_b),
  ifelse(volcano_tbl$significant & volcano_tbl$log2fc < 0, paste("Up in", cell_a), "Not significant")
)
volcano_tbl$neg_log10_padj <- -log10(pmax(volcano_tbl$padj, 1e-300))
top_tbl <- do.call(rbind, lapply(split(volcano_tbl, volcano_tbl$method), label_top_genes))
if (nrow(top_tbl) > 0) {
  top_tbl$direction <- ifelse(top_tbl$log2fc > 0, paste("Up in", cell_b), paste("Up in", cell_a))
  top_tbl$neg_log10_padj <- -log10(pmax(top_tbl$padj, 1e-300))
}

p_volcano <- ggplot(volcano_tbl, aes(x = log2fc, y = neg_log10_padj, color = direction)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.35, color = "grey45") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35, color = "grey45") +
  geom_text_repel(
    data = top_tbl,
    aes(label = gene),
    size = 2.4,
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE,
    segment.color = "grey55",
    segment.size = 0.25
  ) +
  facet_wrap(~ method, scales = "free_y") +
  scale_color_manual(values = c(setNames("#2C7BB6", paste("Up in", cell_a)), "Not significant" = "grey78", setNames("#D7191C", paste("Up in", cell_b)))) +
  labs(
    title = "RNA-seq Volcano Plot (Top 10 per direction labeled)",
    x = "log2 fold change",
    y = expression(-log[10](adjusted~p)),
    color = NULL
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(as.character(defaults$out_dir), "rna_volcano_labeled.pdf"), p_volcano, width = 10, height = 5.5)

plot_combined_go <- function(go_df, out_path, title_text) {
  go_df <- go_df[is.finite(go_df$p.adjust) & is.finite(go_df$fold_enrichment), , drop = FALSE]
  if (nrow(go_df) == 0) return(invisible(NULL))
  go_df <- do.call(rbind, lapply(split(go_df, list(go_df$method, go_df$direction), drop = TRUE), function(x) {
    top <- head(x[order(x$p.adjust), , drop = FALSE], 10)
    top
  }))
  go_df$signed_fold <- ifelse(go_df$direction == paste("Up in", cell_a), -abs(go_df$fold_enrichment), abs(go_df$fold_enrichment))
  go_df <- go_df[order(go_df$method, go_df$signed_fold), , drop = FALSE]
  go_df$term_label <- stringr::str_wrap(go_df$Description, width = 45)
  go_df$term_key <- paste0(go_df$term_label, "___", seq_len(nrow(go_df)))
  go_df$term_key <- factor(go_df$term_key, levels = unique(go_df$term_key))
  p <- ggplot(go_df, aes(x = signed_fold, y = term_key, fill = direction)) +
    geom_col(width = 0.75, alpha = 0.9) +
    facet_wrap(~ method, scales = "free_y") +
    scale_fill_manual(values = c(setNames("#2C7BB6", paste("Up in", cell_a)), setNames("#D7191C", paste("Up in", cell_b)))) +
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

rna_go <- rbind(
  simplified_go_from_rna(deseq_all, "log2FoldChange", "padj", "DESeq2", cell_a, cell_b),
  simplified_go_from_rna(limma_all, "logFC", "adj.P.Val", "limma-voom", cell_a, cell_b)
)
write.csv(rna_go, file.path(as.character(defaults$out_dir), "rna_go_simplified_combined.csv"), row.names = FALSE)
plot_combined_go(rna_go, file.path(as.character(defaults$out_dir), "rna_go_directional_barplot.pdf"), "RNA GO Enrichment (top directional terms)")

deseq_sig_ids <- significant_genes(deseq_tbl)
limma_sig_ids <- significant_genes(limma_tbl)
all_ids <- unique(c(deseq_sig_ids, limma_sig_ids))
venn_input <- cbind(DESeq2 = all_ids %in% deseq_sig_ids, limma_voom = all_ids %in% limma_sig_ids)
pdf(file.path(as.character(defaults$out_dir), "rna_method_overlap_venn.pdf"), width = 6, height = 6)
vennDiagram(vennCounts(venn_input), main = "DE Gene Overlap: DESeq2 vs limma-voom")
grid.text(
  sprintf(
    "DESeq2: %d   limma-voom: %d   overlap: %d",
    length(deseq_sig_ids),
    length(limma_sig_ids),
    length(intersect(deseq_sig_ids, limma_sig_ids))
  ),
  y = unit(0.06, "npc")
)
dev.off()

shared <- merge(
  deseq_tbl[, c("gene_id", "log2fc")],
  limma_tbl[, c("gene_id", "log2fc")],
  by = "gene_id",
  suffixes = c("_deseq2", "_limma")
)
if (nrow(shared) > 0) {
  shared$support_class <- ifelse(
    shared$gene_id %in% deseq_sig_ids & shared$gene_id %in% limma_sig_ids,
    "Called by both",
    ifelse(
      shared$gene_id %in% deseq_sig_ids,
      "DESeq2 only",
      ifelse(shared$gene_id %in% limma_sig_ids, "limma-voom only", "Not DE in either")
    )
  )
  r <- suppressWarnings(cor(shared$log2fc_deseq2, shared$log2fc_limma, method = "pearson"))
  p <- ggplot(shared, aes(x = log2fc_deseq2, y = log2fc_limma, color = support_class)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey60") +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(
      values = c(
        "Called by both" = "#1b9e77",
        "DESeq2 only" = "#d95f02",
        "limma-voom only" = "#7570b3",
        "Not DE in either" = "grey75"
      )
    ) +
    labs(
      title = sprintf("RNA DESeq2 vs limma-voom Concordance (r = %.3f)", r),
      x = "DESeq2 log2FC",
      y = "limma-voom log2FC",
      color = NULL
    ) +
    theme_minimal(base_size = 11)
  ggsave(file.path(as.character(defaults$out_dir), "rna_method_pearson_concordance.pdf"), p, width = 6.5, height = 5.5)
}
