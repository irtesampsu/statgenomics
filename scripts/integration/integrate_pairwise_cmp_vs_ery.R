status <- system2(
  "Rscript",
  c(
    "scripts/integration/integrate_pairwise_generic.R",
    "--cell_a", "CMP",
    "--cell_b", "Erythroblast",
    "--deseq_path", "results/deseq2/de_genes_deseq2.csv",
    "--limma_path", "results/limma_voom/de_genes_limma_voom.csv",
    "--atac_annot_path", "results/atac_seq/pairwise_cmp_vs_ery/DARs_annotated.csv",
    "--out_dir", "results/integration/pairwise_cmp_vs_ery"
  )
)
if (!identical(status, 0L)) {
  stop("Integration CMP vs Erythroblast run failed.")
}
