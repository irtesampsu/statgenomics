dir.create("results/atac_seq/pairwise_cmp_vs_ery", recursive = TRUE, showWarnings = FALSE)

status <- system2(
  "Rscript",
  c(
    "scripts/atac_seq/run_pairwise_atac_generic.R",
    "--cell_a", "CMP",
    "--cell_b", "Erythroblast",
    "--out_dir", "results/atac_seq/pairwise_cmp_vs_ery",
    "--method", "both"
  )
)

if (!identical(status, 0L)) {
  stop("ATAC pairwise CMP vs Erythroblast run failed.")
}
