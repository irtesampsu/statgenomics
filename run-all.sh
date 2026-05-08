# 0) from project root
cd "/Users/irtesam/Courses/STAT 555/statgenomics-project"

# 1) RNA pairwise: CMP vs ERY (both methods)
Rscript scripts/rna_seq/run_pairwise_rna_generic.R \
  --cell_a CMP \
  --cell_b Erythroblast \
  --method deseq2 \
  --out_dir results/deseq2
Rscript scripts/rna_seq/run_pairwise_rna_generic.R \
  --cell_a CMP \
  --cell_b Erythroblast \
  --method limma_voom \
  --out_dir results/limma_voom

# 2) RNA pairwise: CFUE vs ERY (both methods)
Rscript scripts/rna_seq/run_pairwise_rna_generic.R \
  --cell_a CFUE \
  --cell_b Erythroblast \
  --method deseq2 \
  --out_dir results/deseq2/pairwise_cfue_vs_ery
Rscript scripts/rna_seq/run_pairwise_rna_generic.R \
  --cell_a CFUE \
  --cell_b Erythroblast \
  --method limma_voom \
  --out_dir results/limma_voom/pairwise_cfue_vs_ery

# 2b) RNA pairwise: HSC vs ERY (both methods)
Rscript scripts/rna_seq/run_pairwise_rna_generic.R \
  --cell_a HSC \
  --cell_b Erythroblast \
  --method deseq2 \
  --out_dir results/deseq2/pairwise_hsc_vs_ery
Rscript scripts/rna_seq/run_pairwise_rna_generic.R \
  --cell_a HSC \
  --cell_b Erythroblast \
  --method limma_voom \
  --out_dir results/limma_voom/pairwise_hsc_vs_ery

# 3) RNA summary plots (method-overlap venn, directional GO barplot, volcano, concordance)
Rscript scripts/rna_seq/summarize_pairwise_cmp_vs_ery.R \
  --cell_a CMP \
  --cell_b Erythroblast \
  --deseq_path results/deseq2/de_genes_deseq2.csv \
  --limma_path results/limma_voom/de_genes_limma_voom.csv \
  --deseq_all_path results/deseq2/all_genes_deseq2.csv \
  --limma_all_path results/limma_voom/all_genes_limma_voom.csv \
  --deseq_go results/deseq2/go_enrichment_deseq2.csv \
  --limma_go results/limma_voom/go_enrichment_limma_voom.csv \
  --out_dir results/summary/pairwise_cmp_vs_ery
Rscript scripts/rna_seq/summarize_pairwise_cmp_vs_ery.R \
  --cell_a CFUE \
  --cell_b Erythroblast \
  --deseq_path results/deseq2/pairwise_cfue_vs_ery/de_genes_deseq2.csv \
  --limma_path results/limma_voom/pairwise_cfue_vs_ery/de_genes_limma_voom.csv \
  --deseq_all_path results/deseq2/pairwise_cfue_vs_ery/all_genes_deseq2.csv \
  --limma_all_path results/limma_voom/pairwise_cfue_vs_ery/all_genes_limma_voom.csv \
  --deseq_go results/deseq2/pairwise_cfue_vs_ery/go_enrichment_deseq2.csv \
  --limma_go results/limma_voom/pairwise_cfue_vs_ery/go_enrichment_limma_voom.csv \
  --out_dir results/summary/pairwise_cfue_vs_ery
Rscript scripts/rna_seq/summarize_pairwise_cmp_vs_ery.R \
  --cell_a HSC \
  --cell_b Erythroblast \
  --deseq_path results/deseq2/pairwise_hsc_vs_ery/de_genes_deseq2.csv \
  --limma_path results/limma_voom/pairwise_hsc_vs_ery/de_genes_limma_voom.csv \
  --deseq_all_path results/deseq2/pairwise_hsc_vs_ery/all_genes_deseq2.csv \
  --limma_all_path results/limma_voom/pairwise_hsc_vs_ery/all_genes_limma_voom.csv \
  --deseq_go results/deseq2/pairwise_hsc_vs_ery/go_enrichment_deseq2.csv \
  --limma_go results/limma_voom/pairwise_hsc_vs_ery/go_enrichment_limma_voom.csv \
  --out_dir results/summary/pairwise_hsc_vs_ery

# 4) ATAC pairwise: CMP vs ERY (both methods; includes colored venn now)
Rscript scripts/atac_seq/run_pairwise_atac_generic.R \
  --cell_a CMP \
  --cell_b Erythroblast \
  --method both \
  --out_dir results/atac_seq/pairwise_cmp_vs_ery

# 5) ATAC pairwise: CFUE vs ERY (both methods; includes colored venn now)
Rscript scripts/atac_seq/run_pairwise_atac_generic.R \
  --cell_a CFUE \
  --cell_b Erythroblast \
  --method both \
  --out_dir results/atac_seq/pairwise_cfue_vs_ery

# 5b) ATAC pairwise: HSC vs ERY (both methods; includes colored venn now)
Rscript scripts/atac_seq/run_pairwise_atac_generic.R \
  --cell_a HSC \
  --cell_b Erythroblast \
  --method both \
  --out_dir results/atac_seq/pairwise_hsc_vs_ery

# 6) ATAC postprocess for all pairs (directional GO bars, annotation summaries)
Rscript scripts/atac_seq/postprocess_atac_pairwise_generic.R \
  --pair_dir results/atac_seq/pairwise_cmp_vs_ery \
  --cell_a CMP \
  --cell_b Erythroblast
Rscript scripts/atac_seq/postprocess_atac_pairwise_generic.R \
  --pair_dir results/atac_seq/pairwise_cfue_vs_ery \
  --cell_a CFUE \
  --cell_b Erythroblast
Rscript scripts/atac_seq/postprocess_atac_pairwise_generic.R \
  --pair_dir results/atac_seq/pairwise_hsc_vs_ery \
  --cell_a HSC \
  --cell_b Erythroblast

# 7) RNA-ATAC integration for all pairs
Rscript scripts/integration/integrate_pairwise_generic.R \
  --cell_a CMP \
  --cell_b Erythroblast \
  --deseq_path results/deseq2/de_genes_deseq2.csv \
  --limma_path results/limma_voom/de_genes_limma_voom.csv \
  --atac_annot_path results/atac_seq/pairwise_cmp_vs_ery/DARs_annotated.csv \
  --out_dir results/integration/pairwise_cmp_vs_ery
Rscript scripts/integration/integrate_pairwise_generic.R \
  --cell_a CFUE \
  --cell_b Erythroblast \
  --deseq_path results/deseq2/pairwise_cfue_vs_ery/de_genes_deseq2.csv \
  --limma_path results/limma_voom/pairwise_cfue_vs_ery/de_genes_limma_voom.csv \
  --atac_annot_path results/atac_seq/pairwise_cfue_vs_ery/DARs_annotated.csv \
  --out_dir results/integration/pairwise_cfue_vs_ery
Rscript scripts/integration/integrate_pairwise_generic.R \
  --cell_a HSC \
  --cell_b Erythroblast \
  --deseq_path results/deseq2/pairwise_hsc_vs_ery/de_genes_deseq2.csv \
  --limma_path results/limma_voom/pairwise_hsc_vs_ery/de_genes_limma_voom.csv \
  --atac_annot_path results/atac_seq/pairwise_hsc_vs_ery/DARs_annotated.csv \
  --out_dir results/integration/pairwise_hsc_vs_ery

# 8) all-cell clustering snapshots
Rscript scripts/rna_seq/q4_rna_clustering_all_cells.R
Rscript scripts/atac_seq/q8_atac_clustering_all_cells.R

# 9) cross-pair final comparison summary
# Rscript scripts/summary/build_cross_pair_comparison.R