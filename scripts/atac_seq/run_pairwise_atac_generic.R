suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(DESeq2)
  library(edgeR)
  library(ggrepel)
  library(limma)
  library(ggplot2)
  library(scales)
  library(ChIPseeker)
  library(clusterProfiler)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
})

source("scripts/utils/atac_utils.R")

args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  cell_a = "CMP",
  cell_b = "Erythroblast",
  out_dir = NULL,
  method = "both",
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

cell_a <- as.character(defaults$cell_a)
cell_b <- as.character(defaults$cell_b)
if (is.null(defaults$out_dir) || identical(defaults$out_dir, "NULL")) {
  out_dir <- sprintf("results/atac_seq/pairwise_%s_vs_%s", tolower(cell_a), tolower(cell_b))
} else {
  out_dir <- as.character(defaults$out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_info <- get_pairwise_atac(cell_a, cell_b)
gr_list <- import_bigbed_replicates(sample_info)
all_peaks <- build_union_consensus_peaks(gr_list)
atac_quant <- build_atac_count_matrix(sample_info, consensus = all_peaks)
count_mat <- atac_quant$count_mat
condition <- factor(sample_info$condition, levels = c(cell_a, cell_b))
requested_method <- tolower(as.character(defaults$method))
if (!requested_method %in% c("both", "deseq2", "limma_voom")) {
  stop("--method must be one of: both, deseq2, limma_voom")
}
min_count <- as.integer(defaults$min_count)
min_samples <- as.integer(defaults$min_samples)
keep_accessible <- rowSums(count_mat >= min_count) >= min_samples
if (sum(keep_accessible) == 0) {
  stop("No ATAC peaks passed the accessibility filter.")
}
count_mat <- count_mat[keep_accessible, , drop = FALSE]
all_peaks <- all_peaks[keep_accessible]
peak_key <- paste(as.character(seqnames(all_peaks)), start(all_peaks), end(all_peaks), sep = ":")
rownames(count_mat) <- peak_key
peak_metadata <- data.frame(
  peak_key = peak_key,
  chr = as.character(seqnames(all_peaks)),
  start = start(all_peaks),
  end = end(all_peaks),
  stringsAsFactors = FALSE
)

log_counts <- log2(count_mat + 1)
pdf(file.path(out_dir, "replicate_correlation.pdf"))
plot(log_counts[, 1], log_counts[, 2], main = paste(cell_a, "replicates"))
plot(log_counts[, 3], log_counts[, 4], main = paste(cell_b, "replicates"))
dev.off()

dds_for_pca <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = data.frame(row.names = colnames(count_mat), condition),
  design = ~ condition
)
dds_for_pca <- estimateSizeFactors(dds_for_pca)
vsd <- vst(dds_for_pca, blind = TRUE)
pdf(file.path(out_dir, "PCA_plot.pdf"))
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

run_deseq2 <- function() {
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = data.frame(row.names = colnames(count_mat), condition),
    design = ~ condition
  )
  dds <- DESeq(dds, fitType = "local")
  res <- as.data.frame(results(dds, contrast = c("condition", cell_b, cell_a)))
  res$log2FoldChange <- res$log2FoldChange
  res$padj <- res$padj
  res$pvalue <- res$pvalue
  res
}

run_limma_voom <- function() {
  dge <- DGEList(counts = round(count_mat))
  keep <- filterByExpr(dge, group = condition)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ condition)
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  res <- data.frame(
    log2FoldChange = tt$logFC,
    pvalue = tt$P.Value,
    padj = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
  rownames(res) <- rownames(tt)
  res
}

annotate_dars <- function(dar_df) {
  if (nrow(dar_df) == 0) return(data.frame())
  dar_peaks <- GRanges(
    seqnames = dar_df$chr,
    ranges = IRanges(start = dar_df$start, end = dar_df$end)
  )
  peak_anno <- annotatePeak(dar_peaks, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)
  anno_df <- as.data.frame(peak_anno)
  anno_df$log2FoldChange <- dar_df$log2FoldChange
  anno_df$padj <- dar_df$padj
  anno_df$pvalue <- dar_df$pvalue
  anno_df
}

plot_volcano <- function(res_df, method_label, out_file) {
  res_df_clean <- res_df[!is.na(res_df$padj), , drop = FALSE]
  res_df_clean$significant <- res_df_clean$padj < 0.05 & abs(res_df_clean$log2FoldChange) > 1
  res_df_clean$direction <- ifelse(
    res_df_clean$significant & res_df_clean$log2FoldChange > 0,
    paste("More open in", cell_b),
    ifelse(res_df_clean$significant & res_df_clean$log2FoldChange < 0, paste("More open in", cell_a), "Not significant")
  )
  res_df_clean$neg_log10_padj <- -log10(pmax(res_df_clean$padj, 1e-300))
  sig_df <- res_df_clean[res_df_clean$significant, , drop = FALSE]
  top_up <- head(sig_df[sig_df$log2FoldChange > 0, ][order(sig_df[sig_df$log2FoldChange > 0, ]$padj), ], 10)
  top_down <- head(sig_df[sig_df$log2FoldChange < 0, ][order(sig_df[sig_df$log2FoldChange < 0, ]$padj), ], 10)
  top_labeled <- rbind(top_up, top_down)
  if (nrow(top_labeled) > 0) {
    top_peaks <- GRanges(
      seqnames = top_labeled$chr,
      ranges = IRanges(start = top_labeled$start, end = top_labeled$end)
    )
    top_anno <- as.data.frame(annotatePeak(top_peaks, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene))
    peak_gene <- as.character(top_anno$SYMBOL)
    if (!"SYMBOL" %in% colnames(top_anno) || all(is.na(peak_gene) | !nzchar(peak_gene))) {
      entrez_ids <- sub("^\\s*([0-9]+).*$", "\\1", as.character(top_anno$geneId))
      entrez_ids[!grepl("^[0-9]+$", entrez_ids)] <- NA_character_
      mapped <- suppressMessages(AnnotationDbi::mapIds(
        org.Mm.eg.db,
        keys = entrez_ids,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
      ))
      peak_gene <- unname(mapped)
    }
    fallback <- paste(top_labeled$chr, paste0(top_labeled$start, "-", top_labeled$end), sep = ":")
    top_labeled$peak_label <- ifelse(is.na(peak_gene) | !nzchar(peak_gene), fallback, peak_gene)
  }

  p <- ggplot(res_df_clean, aes(x = log2FoldChange, y = neg_log10_padj, color = direction)) +
    geom_point(size = 0.75, alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.35, color = "grey45") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.35, color = "grey45") +
    geom_text_repel(
      data = top_labeled,
      aes(label = peak_label),
      size = 2.1,
      box.padding = 0.35,
      point.padding = 0.2,
      min.segment.length = 0,
      max.overlaps = Inf,
      show.legend = FALSE,
      segment.color = "grey55",
      segment.size = 0.25
    ) +
    scale_color_manual(
      values = setNames(
        c("#2C7BB6", "grey78", "#D7191C"),
        c(paste("More open in", cell_a), "Not significant", paste("More open in", cell_b))
      )
    ) +
    labs(
      title = paste("ATAC volcano -", method_label),
      x = paste("log2 fold change (", cell_b, " vs ", cell_a, ")", sep = ""),
      y = expression(-log[10](adjusted~p)),
      color = NULL
    ) +
    theme_minimal(base_size = 11)
  ggsave(out_file, plot = p, width = 9, height = 6.5)
}

method_results <- list()
if (requested_method %in% c("both", "deseq2")) method_results$deseq2 <- run_deseq2()
if (requested_method %in% c("both", "limma_voom")) method_results$limma_voom <- run_limma_voom()

dar_tables <- list()
anno_tables <- list()
for (m in names(method_results)) {
  res_df <- method_results[[m]]
  res_df$peak_key <- rownames(res_df)
  peak_idx <- match(res_df$peak_key, peak_metadata$peak_key)
  res_df$chr <- peak_metadata$chr[peak_idx]
  res_df$start <- peak_metadata$start[peak_idx]
  res_df$end <- peak_metadata$end[peak_idx]
  write.csv(res_df, file.path(out_dir, paste0("all_peaks_", m, ".csv")), row.names = FALSE)
  dar_idx <- which(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1)
  dar_df <- res_df[dar_idx, , drop = FALSE]
  dar_df$source_method <- rep(m, nrow(dar_df))
  dar_tables[[m]] <- dar_df
  write.csv(dar_df, file.path(out_dir, paste0("DARs_consensus_", m, ".csv")), row.names = FALSE)
  plot_volcano(res_df, m, file.path(out_dir, paste0("volcano_plot_", m, ".pdf")))

  anno_df <- annotate_dars(dar_df)
  if (nrow(anno_df) > 0) {
    anno_df$peak_key <- paste(anno_df$seqnames, anno_df$start, anno_df$end, sep = ":")
    anno_df$source_method <- m
  }
  anno_tables[[m]] <- anno_df
  write.csv(anno_df, file.path(out_dir, paste0("DARs_annotated_", m, ".csv")), row.names = FALSE)
}

if (length(dar_tables) == 2) {
  shared_keys <- intersect(dar_tables$deseq2$peak_key, dar_tables$limma_voom$peak_key)
  shared_dar <- dar_tables$limma_voom[dar_tables$limma_voom$peak_key %in% shared_keys, , drop = FALSE]
  shared_anno <- anno_tables$limma_voom[anno_tables$limma_voom$peak_key %in% shared_keys, , drop = FALSE]
  if (nrow(shared_dar) > 0) {
    write.csv(shared_dar, file.path(out_dir, "DARs_consensus.csv"), row.names = FALSE)
    write.csv(shared_anno, file.path(out_dir, "DARs_annotated.csv"), row.names = FALSE)
    anno_final <- shared_anno
  } else {
    nonempty_dars <- Filter(function(x) nrow(x) > 0, dar_tables)
    nonempty_annos <- Filter(function(x) nrow(x) > 0, anno_tables)
    combined_dar <- if (length(nonempty_dars) > 0) do.call(rbind, nonempty_dars) else data.frame()
    combined_anno <- if (length(nonempty_annos) > 0) do.call(rbind, nonempty_annos) else data.frame()
    write.csv(combined_dar, file.path(out_dir, "DARs_consensus.csv"), row.names = FALSE)
    write.csv(combined_anno, file.path(out_dir, "DARs_annotated.csv"), row.names = FALSE)
    anno_final <- combined_anno
  }
} else {
  only_method <- names(dar_tables)[1]
  write.csv(dar_tables[[only_method]], file.path(out_dir, "DARs_consensus.csv"), row.names = FALSE)
  write.csv(anno_tables[[only_method]], file.path(out_dir, "DARs_annotated.csv"), row.names = FALSE)
  anno_final <- anno_tables[[only_method]]
}

if (length(method_results) == 2) {
  common_peaks <- intersect(rownames(method_results$deseq2), rownames(method_results$limma_voom))
  concord_df <- data.frame(
    peak_id = common_peaks,
    log2fc_deseq2 = method_results$deseq2[common_peaks, "log2FoldChange"],
    log2fc_limma = method_results$limma_voom[common_peaks, "log2FoldChange"],
    stringsAsFactors = FALSE
  )
  concord_df <- concord_df[is.finite(concord_df$log2fc_deseq2) & is.finite(concord_df$log2fc_limma), , drop = FALSE]
  if (nrow(concord_df) > 0) {
    deseq2_dars <- unique(dar_tables$deseq2$peak_key)
    limma_dars <- unique(dar_tables$limma_voom$peak_key)
    concord_df$support_class <- ifelse(
      concord_df$peak_id %in% deseq2_dars & concord_df$peak_id %in% limma_dars,
      "Called by both",
      ifelse(
        concord_df$peak_id %in% deseq2_dars,
        "DESeq2 only",
        ifelse(concord_df$peak_id %in% limma_dars, "limma-voom only", "Not DAR in either")
      )
    )
    pearson_r <- suppressWarnings(cor(concord_df$log2fc_deseq2, concord_df$log2fc_limma, method = "pearson"))
    p <- ggplot(concord_df, aes(x = log2fc_deseq2, y = log2fc_limma, color = support_class)) +
      geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
      geom_vline(xintercept = 0, linewidth = 0.25, color = "grey60") +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(
        values = c(
          "Called by both" = "#1b9e77",
          "DESeq2 only" = "#d95f02",
          "limma-voom only" = "#7570b3",
          "Not DAR in either" = "grey75"
        )
      ) +
      labs(
        title = sprintf("ATAC DESeq2 vs limma-voom Concordance (r = %.3f)", pearson_r),
        x = "DESeq2 log2FC",
        y = "limma-voom log2FC",
        color = NULL
      ) +
      theme_minimal(base_size = 11)
    ggsave(file.path(out_dir, "atac_method_pearson_concordance.pdf"), p, width = 6.5, height = 5.5)
  }

  d1 <- dar_tables$deseq2$peak_key
  d2 <- dar_tables$limma_voom$peak_key
  all_ids <- unique(c(d1, d2))
  if (length(all_ids) > 0) {
    venn_input <- cbind(DESeq2 = all_ids %in% d1, limma_voom = all_ids %in% d2)
    pdf(file.path(out_dir, "atac_method_overlap_venn.pdf"), width = 7.2, height = 6.8)
    grid::grid.newpage()
    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
      stop("Package 'VennDiagram' is required for colored overlap Venn plots. Install with: install.packages('VennDiagram')")
    }
    VennDiagram::draw.pairwise.venn(
      area1 = length(unique(d1)),
      area2 = length(unique(d2)),
      cross.area = length(intersect(unique(d1), unique(d2))),
      category = c("DESeq2", "limma-voom"),
      fill = c("#2C7BB6", "#D7191C"),
      alpha = c(0.45, 0.45), # overlap appears as a distinct blended color
      col = c("#2C7BB6", "#D7191C"),
      lwd = 2,
      cex = 1.45,
      cat.cex = 1.1,
      cat.col = c("#2C7BB6", "#D7191C"),
      cat.pos = c(180, 0),
      cat.dist = c(0.06, 0.06),
      scaled = FALSE
    )
    grid::grid.text(
      "ATAC DAR overlap: DESeq2 vs limma-voom",
      y = grid::unit(0.95, "npc"),
      gp = grid::gpar(cex = 1.15, col = "grey15", fontface = "bold")
    )
    grid::grid.text(
      sprintf(
        "%s vs %s — DAR overlap\nDESeq2: %d   limma-voom: %d   intersection: %d",
        cell_b,
        cell_a,
        length(unique(d1)),
        length(unique(d2)),
        length(intersect(unique(d1), unique(d2)))
      ),
      y = grid::unit(0.05, "npc"),
      gp = grid::gpar(cex = 0.92, col = "grey25")
    )
    dev.off()
  }
}

plot_annotation_pie <- function(anno_df, out_file, title_text) {
  if (nrow(anno_df) == 0 || !"annotation" %in% colnames(anno_df)) return(invisible(NULL))
  cls <- tolower(anno_df$annotation)
  cls <- ifelse(grepl("^promoter", cls), "Promoter/TSS",
    ifelse(grepl("tes|downstream", cls), "TES/Downstream",
      ifelse(grepl("^exon", cls), "Exon",
        ifelse(grepl("^intron", cls), "Intron", "Intergenic/Other")
      )
    )
  )
  pie_counts <- table(cls)
  pdf(out_file, width = 7, height = 6)
  pie(pie_counts, main = title_text)
  dev.off()
}

if (!exists("anno_final")) {
  anno_final <- data.frame()
}

# Keep original consensus pie output for downstream compatibility.
plot_annotation_pie(
  anno_final,
  file.path(out_dir, "atac_annotation_pie.pdf"),
  "ATAC DAR Annotation Composition (Consensus)"
)

# Also emit method-specific annotation pies when method-specific outputs exist.
if ("deseq2" %in% names(anno_tables)) {
  plot_annotation_pie(
    anno_tables$deseq2,
    file.path(out_dir, "atac_annotation_pie_deseq2.pdf"),
    "ATAC DAR Annotation Composition (DESeq2)"
  )
}
if ("limma_voom" %in% names(anno_tables)) {
  plot_annotation_pie(
    anno_tables$limma_voom,
    file.path(out_dir, "atac_annotation_pie_limma_voom.pdf"),
    "ATAC DAR Annotation Composition (limma-voom)"
  )
}
