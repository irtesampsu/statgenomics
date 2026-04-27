suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(ggplot2)
})

dir.create("results/summary", recursive = TRUE, showWarnings = FALSE)

read_results <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

clean_ensembl <- function(ids) {
  sub("\\..*$", "", ids)
}

map_gene_symbols <- function(ids) {
  ids <- clean_ensembl(ids)
  mapped <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
  unname(mapped)
}

prepare_gene_table <- function(df, gene_col, fc_col, padj_col, method_name) {
  out <- df[, c(gene_col, fc_col, padj_col)]
  colnames(out) <- c("gene_id", "fold_change", "padj")
  out$gene_symbol <- map_gene_symbols(out$gene_id)
  out$gene_label <- ifelse(is.na(out$gene_symbol), clean_ensembl(out$gene_id), out$gene_symbol)
  out$method <- method_name
  out
}

top_genes_text <- function(df, direction = c("up", "down"), n = 10) {
  direction <- match.arg(direction)
  if (direction == "up") {
    x <- df[df$fold_change > 0, ]
    x <- x[order(-x$fold_change, x$padj), ]
  } else {
    x <- df[df$fold_change < 0, ]
    x <- x[order(x$fold_change, x$padj), ]
  }
  x <- head(x, n)
  if (nrow(x) == 0) {
    return("None")
  }
  paste(sprintf("%s (%.2f)", x$gene_label, x$fold_change), collapse = ", ")
}

top_go_text <- function(df, n = 10) {
  if (nrow(df) == 0) {
    return("None")
  }
  x <- df[order(df$p.adjust, -df$Count), c("Description", "p.adjust", "Count")]
  x <- head(x, n)
  paste(sprintf("%s (FDR %.2e, n=%s)", x$Description, x$p.adjust, x$Count), collapse = "\n")
}

extract_themes <- function(go_df, n = 8) {
  if (nrow(go_df) == 0) {
    return(character(0))
  }
  head(go_df[order(go_df$p.adjust), "Description"], n)
}

safe_read_results <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read_results(path)
}

integration_line <- function(path, prefix) {
  if (!file.exists(path)) {
    return(paste(prefix, "Not generated yet."))
  }
  x <- readLines(path, warn = FALSE)
  key <- x[grepl(prefix, x, fixed = TRUE)]
  if (length(key) == 0) {
    return(paste(prefix, "Not available."))
  }
  key[1]
}

deseq_de <- read_results("results/deseq2/de_genes_deseq2.csv")
limma_de <- read_results("results/limma_voom/de_genes_limma_voom.csv")
deseq_go <- read_results("results/deseq2/go_enrichment_deseq2.csv")
limma_go <- read_results("results/limma_voom/go_enrichment_limma_voom.csv")
integration_summary_path <- "results/integration/pairwise_cmp_vs_ery/integration_summary.txt"
integration_summary <- if (file.exists(integration_summary_path)) readLines(integration_summary_path, warn = FALSE) else character(0)

deseq_tbl <- prepare_gene_table(deseq_de, "gene_id", "log2FoldChange", "padj", "DESeq2")
limma_tbl <- prepare_gene_table(limma_de, "gene_id", "logFC", "adj.P.Val", "limma-voom")

deseq_tbl$gene_clean <- clean_ensembl(deseq_tbl$gene_id)
limma_tbl$gene_clean <- clean_ensembl(limma_tbl$gene_id)

deseq_ids <- unique(deseq_tbl$gene_clean)
limma_ids <- unique(limma_tbl$gene_clean)
shared_ids <- intersect(deseq_ids, limma_ids)

comparison_df <- merge(
  unique(deseq_tbl[, c("gene_clean", "gene_label", "fold_change", "padj")]),
  unique(limma_tbl[, c("gene_clean", "gene_label", "fold_change", "padj")]),
  by = "gene_clean",
  suffixes = c("_deseq2", "_limma")
)

comparison_df$same_direction <- sign(comparison_df$fold_change_deseq2) == sign(comparison_df$fold_change_limma)
direction_agreement <- if (nrow(comparison_df) > 0) {
  mean(comparison_df$same_direction) * 100
} else {
  NA_real_
}

summary_lines <- c(
  "RNA-seq Summary for CMP vs Erythroblast",
  "=======================================",
  "",
  "Question 1. Differentially expressed genes",
  sprintf("DESeq2 identified %d significant genes (padj < 0.05 and |log2FC| > 1).", nrow(deseq_tbl)),
  sprintf("limma-voom identified %d significant genes (adj.P.Val < 0.05 and |logFC| > 1).", nrow(limma_tbl)),
  sprintf("The two methods shared %d DE genes after matching Ensembl IDs without version suffixes.", length(shared_ids)),
  "",
  "Top genes higher in Erythroblast",
  paste("DESeq2:", top_genes_text(deseq_tbl, "up")),
  paste("limma-voom:", top_genes_text(limma_tbl, "up")),
  "",
  "Top genes higher in CMP",
  paste("DESeq2:", top_genes_text(deseq_tbl, "down")),
  paste("limma-voom:", top_genes_text(limma_tbl, "down")),
  "",
  "Question 2. Functional enrichment of DE genes",
  "Top DESeq2 GO terms:",
  top_go_text(deseq_go),
  "",
  "Top limma-voom GO terms:",
  top_go_text(limma_go),
  "",
  "Question 3. Consistency between DESeq2 and limma-voom",
  sprintf("Shared DE genes: %d", length(shared_ids)),
  if (!is.na(direction_agreement)) sprintf("Direction agreement among shared DE genes: %.1f%%", direction_agreement) else "Direction agreement could not be calculated.",
  paste("Top DESeq2 GO themes:", paste(extract_themes(deseq_go), collapse = "; ")),
  paste("Top limma-voom GO themes:", paste(extract_themes(limma_go), collapse = "; ")),
  "Interpretation: if the main biological themes are similar in both methods, the RNA-seq conclusions are robust even if gene rankings differ.",
  "",
  "Question 4. Sample relationships and clustering",
  "Use the PCA plots, sample correlation heatmaps, hierarchical clustering trees, and top-gene heatmaps from each method.",
  "Interpretation: replicates should cluster together, and CMP samples should separate from Erythroblast samples if the biology is strong and the data quality is good.",
  "",
  "Question 6. Relationship between differential chromatin and expression",
  integration_line(integration_summary_path, "ATAC quantification mode:"),
  integration_line(integration_summary_path, "Integrated genes up/open in Erythroblast:"),
  integration_line(integration_summary_path, "Integrated genes up/open in CMP:"),
  if (length(integration_summary) > 0) {
    integration_summary[grepl("^Interpretation guidance$|^Integration results|^Genes appearing", integration_summary)]
  } else {
    "Integration summary has not been generated yet."
  },
  "",
  "Suggested short report summary",
  sprintf("Both DESeq2 and limma-voom detected strong transcriptional differences between CMP and Erythroblast, with %d and %d significant genes respectively.", nrow(deseq_tbl), nrow(limma_tbl)),
  sprintf("The methods shared %d DE genes, and %.1f%% of shared genes had the same direction of change.", length(shared_ids), direction_agreement),
  "GO enrichment from both methods highlighted erythrocyte differentiation, erythrocyte homeostasis, myeloid cell homeostasis, and ribosome/rRNA processing pathways.",
  "These results are biologically consistent with erythroblast maturation and increased biosynthetic activity during differentiation.",
  if (length(integration_summary) > 0) {
    "RNA-ATAC integration files are available in results/integration/pairwise_cmp_vs_ery and should be used to support regulatory interpretation."
  } else {
    "RNA-ATAC integration has not been summarized yet."
  }
)

writeLines(summary_lines, "results/summary/rna_report_summary.txt")

shared_gene_table <- comparison_df[, c(
  "gene_clean", "gene_label_deseq2", "fold_change_deseq2", "padj_deseq2",
  "fold_change_limma", "padj_limma", "same_direction"
)]
colnames(shared_gene_table) <- c(
  "gene_id", "gene_symbol", "log2FC_deseq2", "padj_deseq2",
  "logFC_limma_voom", "padj_limma_voom", "same_direction"
)
write.csv(shared_gene_table, "results/summary/shared_de_genes_comparison.csv", row.names = FALSE)

top_deseq2_up <- head(deseq_tbl[deseq_tbl$fold_change > 0, c("gene_id", "gene_label", "fold_change", "padj")][order(-deseq_tbl[deseq_tbl$fold_change > 0, "fold_change"]), ], 20)
top_deseq2_down <- head(deseq_tbl[deseq_tbl$fold_change < 0, c("gene_id", "gene_label", "fold_change", "padj")][order(deseq_tbl[deseq_tbl$fold_change < 0, "fold_change"]), ], 20)
top_limma_up <- head(limma_tbl[limma_tbl$fold_change > 0, c("gene_id", "gene_label", "fold_change", "padj")][order(-limma_tbl[limma_tbl$fold_change > 0, "fold_change"]), ], 20)
top_limma_down <- head(limma_tbl[limma_tbl$fold_change < 0, c("gene_id", "gene_label", "fold_change", "padj")][order(limma_tbl[limma_tbl$fold_change < 0, "fold_change"]), ], 20)

write.csv(top_deseq2_up, "results/summary/top20_erythroblast_up_deseq2.csv", row.names = FALSE)
write.csv(top_deseq2_down, "results/summary/top20_cmp_up_deseq2.csv", row.names = FALSE)
write.csv(top_limma_up, "results/summary/top20_erythroblast_up_limma_voom.csv", row.names = FALSE)
write.csv(top_limma_down, "results/summary/top20_cmp_up_limma_voom.csv", row.names = FALSE)

build_directional_top_genes <- function(df, method_name, top_n = 10) {
  up_ery <- df[df$fold_change > 0, c("gene_label", "fold_change", "padj"), drop = FALSE]
  up_ery <- up_ery[order(-up_ery$fold_change, up_ery$padj), , drop = FALSE]
  up_ery <- head(up_ery, top_n)
  up_ery$direction <- "Up in Erythroblast"

  up_cmp <- df[df$fold_change < 0, c("gene_label", "fold_change", "padj"), drop = FALSE]
  up_cmp <- up_cmp[order(up_cmp$fold_change, up_cmp$padj), , drop = FALSE]
  up_cmp <- head(up_cmp, top_n)
  up_cmp$fold_change <- abs(up_cmp$fold_change)
  up_cmp$direction <- "Up in CMP"

  out <- rbind(up_ery, up_cmp)
  if (nrow(out) == 0) {
    return(out)
  }
  out$method <- method_name
  out$gene_label <- as.character(out$gene_label)
  out
}

top_cmp_ery_deseq2 <- build_directional_top_genes(deseq_tbl, "DESeq2", top_n = 10)
top_cmp_ery_limma <- build_directional_top_genes(limma_tbl, "limma-voom", top_n = 10)
top_cmp_ery <- rbind(top_cmp_ery_deseq2, top_cmp_ery_limma)

if (nrow(top_cmp_ery) > 0) {
  top_cmp_ery$direction <- factor(top_cmp_ery$direction, levels = c("Up in CMP", "Up in Erythroblast"))
  top_cmp_ery$method <- factor(top_cmp_ery$method, levels = c("DESeq2", "limma-voom"))
  top_cmp_ery$gene_label <- factor(top_cmp_ery$gene_label, levels = rev(unique(top_cmp_ery$gene_label)))

  write.csv(
    top_cmp_ery,
    "results/summary/top10_upregulated_cmp_vs_ery_comparison.csv",
    row.names = FALSE
  )

  p <- ggplot(top_cmp_ery, aes(x = fold_change, y = gene_label, fill = direction)) +
    geom_col(width = 0.8) +
    facet_wrap(~ method, scales = "free_y") +
    scale_fill_manual(values = c("Up in CMP" = "#67A9CF", "Up in Erythroblast" = "#EF8A62")) +
    labs(
      title = "Top Upregulated Genes in CMP vs Erythroblast",
      x = "Absolute fold change",
      y = "Gene",
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "top"
    )

  ggsave(
    "results/summary/top10_upregulated_cmp_vs_ery_comparison_barplot.pdf",
    plot = p,
    width = 10,
    height = 8
  )
}
