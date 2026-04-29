status <- system2(
  "Rscript",
  c(
    "scripts/atac_seq/postprocess_atac_pairwise_generic.R",
    "--pair_dir", "results/atac_seq/pairwise_cmp_vs_ery"
  )
)
if (!identical(status, 0L)) {
  stop("ATAC postprocess CMP vs Erythroblast run failed.")
}
