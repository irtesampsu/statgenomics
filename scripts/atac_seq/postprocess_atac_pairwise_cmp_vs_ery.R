suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
  library(org.Mm.eg.db)
})

in_path <- "results/atac_seq/pairwise_cmp_vs_ery/DARs_annotated.csv"
out_dir <- "results/atac_seq/pairwise_cmp_vs_ery"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dar <- read.csv(in_path, stringsAsFactors = FALSE, check.names = FALSE)

parse_entrez <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  out <- sub("^\\s*([0-9]+).*$", "\\1", x)
  out[!grepl("^[0-9]+$", out)] <- NA_character_
  out
}

is_proximal_annotation <- function(annotation_vec) {
  ann <- tolower(annotation_vec)
  !grepl("distal intergenic", ann, fixed = TRUE)
}

coarse_annotation <- function(annotation_vec) {
  ann <- tolower(annotation_vec)
  out <- ifelse(grepl("^promoter", ann), "Promoter",
    ifelse(grepl("distal intergenic", ann, fixed = TRUE), "Distal intergenic",
      ifelse(grepl("^intron", ann), "Intron",
        ifelse(grepl("^exon", ann), "Exon",
          ifelse(grepl("utr", ann, fixed = TRUE), "UTR",
            ifelse(grepl("^downstream", ann), "Downstream", "Other")
          )
        )
      )
    )
  )
  factor(out, levels = c("Promoter", "Intron", "Exon", "UTR", "Downstream", "Distal intergenic", "Other"))
}

dar$entrez_id <- parse_entrez(dar$geneId)
dar$is_proximal <- is_proximal_annotation(dar$annotation)
dar$direction <- ifelse(dar$log2FoldChange > 0, "More open in Erythroblast", "More open in CMP")
dar$annotation_group <- coarse_annotation(dar$annotation)
dar$abs_distance_to_tss <- abs(dar$distanceToTSS)
dar$log10_abs_distance_to_tss <- log10(dar$abs_distance_to_tss + 1)

# Stringent subset for defensible reporting under bigBed-score fallback.
dar_stringent <- dar[
  !is.na(dar$padj) &
    !is.na(dar$log2FoldChange) &
    dar$padj < 0.01 &
    abs(dar$log2FoldChange) >= 2 &
    dar$is_proximal,
  ,
  drop = FALSE
]

dar_stringent <- dar_stringent[order(dar_stringent$padj, -abs(dar_stringent$log2FoldChange)), ]
write.csv(dar_stringent, file.path(out_dir, "DARs_stringent.csv"), row.names = FALSE)

if (nrow(dar) > 0) {
  p <- ggplot(
    dar[!is.na(dar$log10_abs_distance_to_tss), ],
    aes(x = log10_abs_distance_to_tss, fill = direction)
  ) +
    geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
    facet_wrap(~ direction, ncol = 1) +
    scale_fill_manual(values = c("More open in CMP" = "#67A9CF", "More open in Erythroblast" = "#EF8A62")) +
    labs(
      title = "ATAC DAR Distance to Nearest TSS",
      x = "log10(|distance to TSS| + 1)",
      y = "DAR count",
      fill = NULL
    ) +
    theme_minimal(base_size = 11)

  ggsave(file.path(out_dir, "distance_to_tss_histogram.pdf"), plot = p, width = 8, height = 7)
}

if (nrow(dar_stringent) > 0) {
  p <- ggplot(
    dar_stringent[!is.na(dar_stringent$log10_abs_distance_to_tss), ],
    aes(x = direction, y = log10_abs_distance_to_tss, fill = direction)
  ) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.8) +
    scale_fill_manual(values = c("More open in CMP" = "#67A9CF", "More open in Erythroblast" = "#EF8A62")) +
    labs(
      title = "Stringent ATAC DAR Distance to TSS by Direction",
      x = NULL,
      y = "log10(|distance to TSS| + 1)",
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")

  ggsave(file.path(out_dir, "distance_to_tss_violin_stringent.pdf"), plot = p, width = 7, height = 5)
}

annotation_summary <- as.data.frame(table(
  annotation_group = dar$annotation_group,
  direction = dar$direction
))
annotation_summary <- annotation_summary[annotation_summary$Freq > 0, , drop = FALSE]

if (nrow(annotation_summary) > 0) {
  p <- ggplot(annotation_summary, aes(x = annotation_group, y = Freq, fill = direction)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("More open in CMP" = "#67A9CF", "More open in Erythroblast" = "#EF8A62")) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
    labs(
      title = "ATAC DAR Annotation Composition by Direction",
      x = NULL,
      y = "Proportion of DARs",
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(out_dir, "annotation_composition_by_direction.pdf"), plot = p, width = 8, height = 5)
}

direction_counts <- data.frame(
  direction = c("More open in CMP", "More open in Erythroblast"),
  baseline = c(sum(dar$log2FoldChange < 0, na.rm = TRUE), sum(dar$log2FoldChange > 0, na.rm = TRUE)),
  stringent = c(sum(dar_stringent$log2FoldChange < 0, na.rm = TRUE), sum(dar_stringent$log2FoldChange > 0, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

direction_counts_long <- rbind(
  data.frame(direction = direction_counts$direction, subset = "Baseline", count = direction_counts$baseline),
  data.frame(direction = direction_counts$direction, subset = "Stringent", count = direction_counts$stringent)
)

p <- ggplot(direction_counts_long, aes(x = subset, y = count, fill = direction)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("More open in CMP" = "#67A9CF", "More open in Erythroblast" = "#EF8A62")) +
  labs(
    title = "ATAC DAR Counts by Direction",
    x = NULL,
    y = "DAR count",
    fill = NULL
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(out_dir, "directional_dar_counts.pdf"), plot = p, width = 7, height = 5)

run_directional_go <- function(entrez_ids, universe_ids, label) {
  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids) & nzchar(entrez_ids)])
  universe_ids <- unique(universe_ids[!is.na(universe_ids) & nzchar(universe_ids)])

  out_path <- file.path(out_dir, paste0("GO_enrichment_", label, ".csv"))
  if (length(entrez_ids) < 10 || length(universe_ids) < 10) {
    write.csv(data.frame(), out_path, row.names = FALSE)
    return(invisible(NULL))
  }

  ego <- enrichGO(
    gene = entrez_ids,
    universe = universe_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH"
  )
  write.csv(as.data.frame(ego), out_path, row.names = FALSE)
}

universe_entrez <- parse_entrez(dar$geneId)
open_in_ery <- parse_entrez(dar_stringent$geneId[dar_stringent$log2FoldChange > 0])
open_in_cmp <- parse_entrez(dar_stringent$geneId[dar_stringent$log2FoldChange < 0])

run_directional_go(open_in_ery, universe_entrez, "stringent_open_in_ery")
run_directional_go(open_in_cmp, universe_entrez, "stringent_open_in_cmp")

summary_lines <- c(
  "ATAC stringent post-processing summary",
  "=====================================",
  "",
  "Base file: results/atac_seq/pairwise_cmp_vs_ery/DARs_annotated.csv",
  "Stringent filter definition: padj < 0.01, |log2FC| >= 2, non-distal annotation",
  sprintf("Total DAR rows (annotated): %d", nrow(dar)),
  sprintf("Stringent DAR rows: %d", nrow(dar_stringent)),
  sprintf("Stringent DARs more open in Erythroblast: %d", sum(dar_stringent$log2FoldChange > 0, na.rm = TRUE)),
  sprintf("Stringent DARs more open in CMP: %d", sum(dar_stringent$log2FoldChange < 0, na.rm = TRUE)),
  sprintf("Unique stringent genes (open in Ery): %d", length(unique(open_in_ery[!is.na(open_in_ery)]))),
  sprintf("Unique stringent genes (open in CMP): %d", length(unique(open_in_cmp[!is.na(open_in_cmp)]))),
  "",
  "Additional exploratory figures generated:",
  "- results/atac_seq/pairwise_cmp_vs_ery/distance_to_tss_histogram.pdf",
  "- results/atac_seq/pairwise_cmp_vs_ery/distance_to_tss_violin_stringent.pdf",
  "- results/atac_seq/pairwise_cmp_vs_ery/annotation_composition_by_direction.pdf",
  "- results/atac_seq/pairwise_cmp_vs_ery/directional_dar_counts.pdf",
  "",
  "Outputs generated:",
  "- results/atac_seq/pairwise_cmp_vs_ery/DARs_stringent.csv",
  "- results/atac_seq/pairwise_cmp_vs_ery/GO_enrichment_stringent_open_in_ery.csv",
  "- results/atac_seq/pairwise_cmp_vs_ery/GO_enrichment_stringent_open_in_cmp.csv",
  "",
  "Interpretation note:",
  "These outputs help visualize promoter/distal balance and TSS proximity, but they do not resolve the core ATAC quantification limitation.",
  "Under bigBed score-based fallback, even the stringent subset remains exploratory rather than a robust differential-accessibility call set."
)

writeLines(summary_lines, file.path(out_dir, "stringent_analysis_notes.txt"))
