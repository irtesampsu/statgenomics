suppressPackageStartupMessages({
  library(tools)
})

project_rna_scriptseq_samples <- data.frame(
  sample = c(
    "HSC_rep1", "HSC_rep2",
    "CMP_rep1", "CMP_rep2",
    "CFUE_rep1", "CFUE_rep2",
    "Erythroblast_rep1", "Erythroblast_rep2"
  ),
  condition = c(
    "HSC", "HSC",
    "CMP", "CMP",
    "CFUE", "CFUE",
    "Erythroblast", "Erythroblast"
  ),
  assay = "ScriptSeq",
  file = c(
    "data/rna_seq/HSC_rep1_ENCFF247FEJ.tsv",
    "data/rna_seq/HSC_rep2_ENCFF064MKY.tsv",
    "data/rna_seq/CMP_rep1_ENCFF623OLU.tsv",
    "data/rna_seq/CMP_rep2_ENCFF691MHW.tsv",
    "data/rna_seq/CFUE_rep1_ENCFF667IDY.tsv",
    "data/rna_seq/CFUE_rep2_ENCFF655LMK.tsv",
    "data/rna_seq/Erythroblast_rep1_ENCFF342WUL.tsv",
    "data/rna_seq/Erythroblast_rep2_ENCFF858JHF.tsv"
  ),
  stringsAsFactors = FALSE
)

project_atac_samples <- data.frame(
  sample = c(
    "HSC_rep1", "HSC_rep2",
    "CMP_rep1", "CMP_rep2",
    "CFUE_rep1", "CFUE_rep2",
    "Erythroblast_rep1", "Erythroblast_rep2"
  ),
  condition = c(
    "HSC", "HSC",
    "CMP", "CMP",
    "CFUE", "CFUE",
    "Erythroblast", "Erythroblast"
  ),
  bam_file = c(
    "data/bam/atac_seq/HSC_rep1_ENCFF250YAL.bam",
    "data/bam/atac_seq/HSC_rep2_ENCFF958EPJ.bam",
    "data/bam/atac_seq/CMP_rep1_ENCFF711QAL.bam",
    "data/bam/atac_seq/CMP_rep2_ENCFF620WGW.bam",
    "data/bam/atac_seq/CFUE_rep1_ENCFF909QFQ.bam",
    "data/bam/atac_seq/CFUE_rep2_ENCFF780SSI.bam",
    "data/bam/atac_seq/Erythroblast_rep1_ENCFF199ZJX.bam",
    "data/bam/atac_seq/Erythroblast_rep2_ENCFF535OJU.bam"
  ),
  file = c(
    "data/atac_seq/HSC_rep1_ENCFF662DYG.bigBed",
    "data/atac_seq/HSC_rep2_ENCFF255IVU.bigBed",
    "data/atac_seq/CMP_rep1_ENCFF832UUS.bigBed",
    "data/atac_seq/CMP_rep2_ENCFF343PTQ.bigBed",
    "data/atac_seq/CFUE_rep1_ENCFF796ZSB.bigBed",
    "data/atac_seq/CFUE_rep2_ENCFF599ZDJ.bigBed",
    "data/atac_seq/Erythroblast_rep1_ENCFF181AMY.bigBed",
    "data/atac_seq/Erythroblast_rep2_ENCFF616EWK.bigBed"
  ),
  stringsAsFactors = FALSE
)

validate_project_files <- function(paths, label) {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop(
      paste(
        label,
        "is missing required files:",
        paste(missing_paths, collapse = ", ")
      )
    )
  }
}

get_pairwise_rna_scriptseq <- function(condition_a, condition_b) {
  keep <- project_rna_scriptseq_samples$condition %in% c(condition_a, condition_b)
  sample_info <- project_rna_scriptseq_samples[keep, , drop = FALSE]
  sample_info <- sample_info[order(match(sample_info$condition, c(condition_a, condition_b)), sample_info$sample), ]
  validate_project_files(sample_info$file, "RNA pairwise analysis")
  sample_info
}

get_all_rna_scriptseq <- function() {
  validate_project_files(project_rna_scriptseq_samples$file, "RNA all-cell analysis")
  project_rna_scriptseq_samples
}

get_pairwise_atac <- function(condition_a, condition_b) {
  keep <- project_atac_samples$condition %in% c(condition_a, condition_b)
  sample_info <- project_atac_samples[keep, , drop = FALSE]
  sample_info <- sample_info[order(match(sample_info$condition, c(condition_a, condition_b)), sample_info$sample), ]
  validate_project_files(sample_info$file, "ATAC pairwise analysis")
  sample_info
}

get_all_atac <- function() {
  validate_project_files(project_atac_samples$file, "ATAC all-cell analysis")
  project_atac_samples
}
