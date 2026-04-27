# =========================
# ATAC-seq with replicate consensus (FIXED + VALIDATION)
# =========================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(DESeq2)
  library(ggplot2)          # NEW
  library(pheatmap)         # NEW
  library(ChIPseeker)       # NEW
  library(clusterProfiler)  # NEW
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
})

source("scripts/utils/atac_utils.R")

dir.create("results/atac_seq/pairwise_cmp_vs_ery", recursive = TRUE, showWarnings = FALSE)

# ---- input ----
sample_info <- get_pairwise_atac("CMP", "Erythroblast")

# ---- import peaks ----
gr_list <- import_bigbed_replicates(sample_info)

# ---- condition-level peak sets ----
# Use union of replicate peaks per condition; strict interval intersection can
# be overly conservative when replicate peak boundaries are slightly shifted.
all_peaks <- build_union_consensus_peaks(gr_list)
atac_quant <- build_atac_count_matrix(sample_info, consensus = all_peaks)
count_mat <- atac_quant$count_mat
quant_method <- atac_quant$quant_method
condition <- factor(sample_info$condition, levels = c("CMP", "Erythroblast"))

# ---- DESeq2 ----
dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = data.frame(row.names = colnames(count_mat), condition),
  design = ~ condition
)

dds <- DESeq(dds, fitType = "local")   # FIX: better for ATAC
res <- results(dds, contrast = c("condition","Erythroblast","CMP"))

# ---- attach coordinates ----
res_df <- as.data.frame(res)
res_df$chr   <- as.character(seqnames(all_peaks))
res_df$start <- start(all_peaks)
res_df$end   <- end(all_peaks)

# ---- filter DARs ----
dar_idx <- which(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1)
dar <- res_df[dar_idx, ]
write.csv(dar, "results/atac_seq/pairwise_cmp_vs_ery/DARs_consensus.csv", row.names = FALSE)
writeLines(
  c(
    "ATAC pairwise analysis notes",
    if (quant_method == "bam_counts") {
      "Inputs: processed bigBed peaks for consensus definition plus BAM files for read counting"
    } else {
      "Inputs: processed bigBed peak files only"
    },
    "Consensus strategy: union of replicate peaks per condition, then union across conditions",
    if (quant_method == "bam_counts") {
      "Quantification: BAM read counts per consensus peak using summarizeOverlaps"
    } else {
      "Quantification: summed bigBed peak scores over overlaps per consensus peak"
    },
    if (quant_method == "bam_counts") {
      "Interpretation: BAM-based peak counts are preferred and should be more defensible than bigBed score-derived pseudo-counts."
    } else {
      "Limitation: DESeq2 is being applied to score-derived pseudo-counts, not BAM-derived read counts."
    },
    if (quant_method == "bam_counts") {
      "BAM mode is active because all expected ATAC BAM files for this comparison were found."
    } else {
      "Interpretation: treat differential accessibility statistics as approximate unless BAM-based counting is added."
    }
  ),
  "results/atac_seq/pairwise_cmp_vs_ery/README.txt"
)

# =========================
# VALIDATION SECTION
# =========================

# ---- 1. Replicate correlation ----
log_counts <- log2(count_mat + 1)

pdf("results/atac_seq/pairwise_cmp_vs_ery/replicate_correlation.pdf")
plot(log_counts[,1], log_counts[,2], main="CMP replicates")
plot(log_counts[,3], log_counts[,4], main="Erythroblast replicates")
dev.off()

cor_cmp <- cor(log_counts[,1], log_counts[,2])
cor_ery <- cor(log_counts[,3], log_counts[,4])
print(cor_cmp)
print(cor_ery)

# ---- 2. PCA ----
vsd <- vst(dds, blind=TRUE)

pdf("results/atac_seq/pairwise_cmp_vs_ery/PCA_plot.pdf")
print(plotPCA(vsd, intgroup="condition"))
dev.off()

# ---- 3. p-value distribution ----
pdf("results/atac_seq/pairwise_cmp_vs_ery/pvalue_hist.pdf")
hist(res$pvalue, breaks=50, main="P-value distribution")
dev.off()

# ---- 4. Volcano plot ----
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

res_df_clean <- res_df[!is.na(res_df$padj), ]

p <- ggplot(res_df_clean, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(size=1) +
  theme_minimal()

print(p)  # shows in R

ggsave("results/atac_seq/pairwise_cmp_vs_ery/volcano_plot.pdf", plot = p, width = 6, height = 5)


# ---- 5. Peak annotation on DARs only ----
if (length(dar_idx) > 0) {
  dar_peaks <- all_peaks[dar_idx]
  peakAnno_dar <- annotatePeak(
    dar_peaks,
    TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene
  )
  dar_anno_df <- as.data.frame(peakAnno_dar)
  dar_anno_df$log2FoldChange <- dar$log2FoldChange
  dar_anno_df$padj <- dar$padj
  write.csv(dar_anno_df, "results/atac_seq/pairwise_cmp_vs_ery/DARs_annotated.csv", row.names = FALSE)

  pdf("results/atac_seq/pairwise_cmp_vs_ery/peak_annotation_pie.pdf")
  plotAnnoPie(peakAnno_dar)
  dev.off()

  # ---- 6. GO enrichment for genes near DARs ----
  genes <- unique(dar_anno_df$geneId)
  genes <- genes[!is.na(genes)]

  if (length(genes) > 0) {
    ego <- enrichGO(
      gene = genes,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP"
    )
    write.csv(as.data.frame(ego), "results/atac_seq/pairwise_cmp_vs_ery/GO_enrichment.csv", row.names = FALSE)
  } else {
    write.csv(data.frame(), "results/atac_seq/pairwise_cmp_vs_ery/GO_enrichment.csv", row.names = FALSE)
  }
} else {
  write.csv(data.frame(), "results/atac_seq/pairwise_cmp_vs_ery/DARs_annotated.csv", row.names = FALSE)
  write.csv(data.frame(), "results/atac_seq/pairwise_cmp_vs_ery/GO_enrichment.csv", row.names = FALSE)
}

# =========================
# DONE
# =========================
