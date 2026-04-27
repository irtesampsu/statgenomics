suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

dir.create("results/integration/pairwise_cmp_vs_ery", recursive = TRUE, showWarnings = FALSE)

clean_ensembl <- function(ids) {
  sub("\\..*$", "", ids)
}

map_ensembl_to_symbol <- function(ids) {
  ids <- unique(clean_ensembl(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) {
    return(character(0))
  }
  mapped <- tryCatch(
    suppressMessages(
      AnnotationDbi::mapIds(
        org.Mm.eg.db,
        keys = ids,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    ),
    error = function(e) setNames(rep(NA_character_, length(ids)), ids)
  )
  unname(mapped[ids])
}

map_entrez_to_symbol <- function(ids) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) {
    return(character(0))
  }
  mapped <- tryCatch(
    suppressMessages(
      AnnotationDbi::mapIds(
        org.Mm.eg.db,
        keys = ids,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
      )
    ),
    error = function(e) setNames(rep(NA_character_, length(ids)), ids)
  )
  unname(mapped[ids])
}

map_entrez_to_ensembl <- function(ids) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) {
    return(character(0))
  }
  mapped <- tryCatch(
    suppressMessages(
      AnnotationDbi::mapIds(
        org.Mm.eg.db,
        keys = ids,
        column = "ENSEMBL",
        keytype = "ENTREZID",
        multiVals = "first"
      )
    ),
    error = function(e) setNames(rep(NA_character_, length(ids)), ids)
  )
  unname(mapped[ids])
}

read_results <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

read_atac_mode <- function(path) {
  if (!file.exists(path)) {
    return("unknown")
  }
  lines <- readLines(path, warn = FALSE)
  if (any(grepl("BAM read counts", lines, fixed = TRUE))) {
    return("bam_counts")
  }
  if (any(grepl("bigBed", lines, fixed = TRUE))) {
    return("bigBed_score_fallback")
  }
  "unknown"
}

split_direction <- function(df, fc_col, id_col, positive_label, negative_label) {
  up_pos <- unique(df[df[[fc_col]] > 0, id_col])
  up_neg <- unique(df[df[[fc_col]] < 0, id_col])
  list(
    positive = up_pos,
    negative = up_neg,
    positive_label = positive_label,
    negative_label = negative_label
  )
}

prepare_rna_method <- function(df, fc_col, padj_col, method_name) {
  out <- df[, c("gene_id", fc_col, padj_col)]
  colnames(out) <- c("gene_id", "fold_change", "padj")
  out$gene_id <- clean_ensembl(out$gene_id)
  out$gene_symbol <- map_ensembl_to_symbol(out$gene_id)
  out$gene_label <- ifelse(is.na(out$gene_symbol), out$gene_id, out$gene_symbol)
  out$method <- method_name
  out
}

prepare_atac_gene_table <- function(atac_df) {
  raw_gene_id <- as.character(atac_df$geneId)
  raw_gene_id[is.na(raw_gene_id)] <- ""
  # ChIPseeker geneId may be single IDs or delimited lists. Keep first numeric ENTREZ token.
  parsed_gene_id <- sub("^\\s*([0-9]+).*$", "\\1", raw_gene_id)
  parsed_gene_id[!grepl("^[0-9]+$", parsed_gene_id)] <- NA_character_
  atac_df$entrez_id <- parsed_gene_id
  atac_df <- atac_df[!is.na(atac_df$entrez_id), , drop = FALSE]
  if (nrow(atac_df) == 0) {
    return(data.frame(
      entrez_id = character(0),
      ensembl_id = character(0),
      gene_symbol = character(0),
      mean_log2FoldChange = numeric(0),
      best_padj = numeric(0),
      n_dars = integer(0),
      representative_annotation = character(0),
      gene_label = character(0),
      stringsAsFactors = FALSE
    ))
  }

  split_by_gene <- split(atac_df, atac_df$entrez_id)
  out <- do.call(rbind, lapply(names(split_by_gene), function(gene_id) {
    sub_df <- split_by_gene[[gene_id]]
    best_padj <- suppressWarnings(min(sub_df$padj, na.rm = TRUE))
    if (!is.finite(best_padj)) best_padj <- NA_real_
    data.frame(
      entrez_id = gene_id,
      mean_log2FoldChange = mean(sub_df$log2FoldChange, na.rm = TRUE),
      best_padj = best_padj,
      n_dars = nrow(sub_df),
      representative_annotation = as.character(sub_df$annotation[1]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL

  entrez_ids <- out$entrez_id
  ensembl_mapped <- map_entrez_to_ensembl(entrez_ids)
  symbol_mapped <- map_entrez_to_symbol(entrez_ids)
  out$ensembl_id <- ensembl_mapped
  out$gene_symbol <- symbol_mapped
  out$gene_label <- ifelse(is.na(out$gene_symbol), out$entrez_id, out$gene_symbol)
  out
}

run_go <- function(entrez_ids, label) {
  entrez_ids <- unique(as.character(entrez_ids))
  entrez_ids <- entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)]

  if (length(entrez_ids) < 10) {
    write.csv(data.frame(), file.path("results/integration/pairwise_cmp_vs_ery", paste0(label, "_go.csv")), row.names = FALSE)
    return(data.frame())
  }

  ego <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH"
  )

  ego_df <- as.data.frame(ego)
  write.csv(ego_df, file.path("results/integration/pairwise_cmp_vs_ery", paste0(label, "_go.csv")), row.names = FALSE)
  ego_df
}

deseq_de <- read_results("results/deseq2/de_genes_deseq2.csv")
limma_de <- read_results("results/limma_voom/de_genes_limma_voom.csv")
atac_annot <- read_results("results/atac_seq/pairwise_cmp_vs_ery/DARs_annotated.csv")

deseq_tbl <- prepare_rna_method(deseq_de, "log2FoldChange", "padj", "DESeq2")
limma_tbl <- prepare_rna_method(limma_de, "logFC", "adj.P.Val", "limma-voom")
atac_tbl <- prepare_atac_gene_table(atac_annot)

deseq_dir <- split_direction(deseq_tbl, "fold_change", "gene_id", "Up in Erythroblast", "Up in CMP")
limma_dir <- split_direction(limma_tbl, "fold_change", "gene_id", "Up in Erythroblast", "Up in CMP")
shared_rna_up_ery <- intersect(deseq_dir$positive, limma_dir$positive)
shared_rna_up_cmp <- intersect(deseq_dir$negative, limma_dir$negative)

atac_ery_ensembl <- clean_ensembl(na.omit(atac_tbl$ensembl_id[atac_tbl$mean_log2FoldChange > 0]))
atac_cmp_ensembl <- clean_ensembl(na.omit(atac_tbl$ensembl_id[atac_tbl$mean_log2FoldChange < 0]))

integrated_shared_ery <- intersect(shared_rna_up_ery, atac_ery_ensembl)
integrated_shared_cmp <- intersect(shared_rna_up_cmp, atac_cmp_ensembl)

build_integrated_table <- function(rna_ids, direction_label) {
  if (length(rna_ids) == 0) {
    return(data.frame())
  }

  rna_sub <- merge(
    deseq_tbl[, c("gene_id", "gene_label", "fold_change", "padj")],
    limma_tbl[, c("gene_id", "fold_change", "padj", "gene_symbol")],
    by = "gene_id",
    suffixes = c("_deseq2", "_limma")
  )
  colnames(rna_sub) <- c(
    "gene_id", "gene_label", "log2FC_deseq2", "padj_deseq2",
    "logFC_limma", "padj_limma", "gene_symbol_limma"
  )

  atac_sub <- atac_tbl[, c("ensembl_id", "entrez_id", "gene_symbol", "mean_log2FoldChange", "best_padj", "n_dars", "representative_annotation")]
  colnames(atac_sub)[1] <- "gene_id"
  atac_sub$gene_id <- clean_ensembl(atac_sub$gene_id)

  out <- merge(
    rna_sub,
    atac_sub,
    by = "gene_id"
  )
  out <- out[out$gene_id %in% rna_ids, , drop = FALSE]
  out$direction <- direction_label
  out <- out[order(-abs(out$log2FC_deseq2), out$best_padj), ]
  out
}

integrated_ery_tbl <- build_integrated_table(integrated_shared_ery, "Up in Erythroblast and more open in Erythroblast")
integrated_cmp_tbl <- build_integrated_table(integrated_shared_cmp, "Up in CMP and more open in CMP")
integrated_tbl <- rbind(integrated_ery_tbl, integrated_cmp_tbl)

write.csv(integrated_ery_tbl, "results/integration/pairwise_cmp_vs_ery/integrated_shared_up_in_ery.csv", row.names = FALSE)
write.csv(integrated_cmp_tbl, "results/integration/pairwise_cmp_vs_ery/integrated_shared_up_in_cmp.csv", row.names = FALSE)
write.csv(integrated_tbl, "results/integration/pairwise_cmp_vs_ery/integrated_shared_all.csv", row.names = FALSE)

go_ery <- run_go(unique(integrated_ery_tbl$entrez_id), "integrated_up_in_ery")
go_cmp <- run_go(unique(integrated_cmp_tbl$entrez_id), "integrated_up_in_cmp")

atac_mode <- read_atac_mode("results/atac_seq/pairwise_cmp_vs_ery/README.txt")

summary_lines <- c(
  "RNA-ATAC Integration Summary for CMP vs Erythroblast",
  "====================================================",
  "",
  sprintf("ATAC quantification mode: %s", atac_mode),
  sprintf("DESeq2 significant DE genes: %d", nrow(deseq_tbl)),
  sprintf("limma-voom significant DE genes: %d", nrow(limma_tbl)),
  sprintf("Genes significant in both RNA methods and up in Erythroblast: %d", length(shared_rna_up_ery)),
  sprintf("Genes significant in both RNA methods and up in CMP: %d", length(shared_rna_up_cmp)),
  sprintf("ATAC-linked genes more open in Erythroblast: %d", length(unique(atac_ery_ensembl))),
  sprintf("ATAC-linked genes more open in CMP: %d", length(unique(atac_cmp_ensembl))),
  sprintf("Integrated genes up/open in Erythroblast: %d", nrow(integrated_ery_tbl)),
  sprintf("Integrated genes up/open in CMP: %d", nrow(integrated_cmp_tbl)),
  "",
  "Interpretation guidance",
  if (atac_mode == "bam_counts") {
    "Integration results are based on BAM-derived ATAC counts and are more suitable for stronger biological claims."
  } else {
    "Integration results use the current bigBed-score ATAC fallback and should be described as exploratory."
  },
  "Genes appearing in the integrated tables have matched directionality between RNA expression and nearby chromatin accessibility.",
  "These genes are stronger candidates for regulatory interpretation than RNA-only or ATAC-only findings.",
  "",
  "Suggested reporting angle",
  "Focus on whether erythroid genes tend to be both upregulated and more accessible in Erythroblast, and whether progenitor-associated genes tend to be both upregulated and more accessible in CMP.",
  "If BAM mode is not active yet, avoid treating the ATAC integration counts as definitive."
)

writeLines(summary_lines, "results/integration/pairwise_cmp_vs_ery/integration_summary.txt")
