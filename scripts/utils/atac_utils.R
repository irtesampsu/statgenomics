source("scripts/utils/project_metadata.R")

import_bigbed_replicates <- function(sample_info) {
  condition_levels <- unique(as.character(sample_info$condition))
  condition_factor <- factor(sample_info$condition, levels = condition_levels)
  files_by_condition <- split(sample_info$file, condition_factor)
  lapply(files_by_condition, function(reps) {
    lapply(reps, function(f) {
      gr <- rtracklayer::import(f, format = "bigBed")
      GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
      gr
    })
  })
}

build_union_consensus_peaks <- function(gr_list) {
  consensus_per_condition <- lapply(gr_list, function(reps) {
    reduce(do.call(c, reps))
  })
  all_peaks <- Reduce(c, consensus_per_condition)
  reduce(all_peaks)
}

get_bigbed_signal <- function(gr, consensus) {
  hits <- findOverlaps(consensus, gr)
  signal <- numeric(length(consensus))

  if ("score" %in% colnames(mcols(gr))) {
    score_sums <- rowsum(
      x = mcols(gr)$score[subjectHits(hits)],
      group = queryHits(hits),
      reorder = FALSE
    )
    signal[as.integer(rownames(score_sums))] <- as.numeric(score_sums[, 1])
  } else {
    warning("BigBed file does not contain score column, using 1 for all overlaps")
    signal[unique(queryHits(hits))] <- 1
  }

  signal
}

build_atac_score_matrix <- function(gr_list, consensus, sample_names) {
  all_samples <- unlist(gr_list, recursive = FALSE, use.names = FALSE)
  if (length(all_samples) != length(sample_names)) {
    stop("ATAC replicate count does not match sample metadata.")
  }
  score_mat <- sapply(all_samples, get_bigbed_signal, consensus = consensus)
  score_mat <- round(as.matrix(score_mat))
  colnames(score_mat) <- sample_names
  score_mat
}

build_atac_count_matrix <- function(sample_info, consensus) {
  gr_list <- import_bigbed_replicates(sample_info)
  list(
    count_mat = build_atac_score_matrix(gr_list, consensus = consensus, sample_names = sample_info$sample),
    quant_method = "bigBed_score_fallback",
    gr_list = gr_list
  )
}
