suppressPackageStartupMessages({
  library(utils)
})

id_file <- "encode ids.txt"
base_url <- "https://www.encodeproject.org/files"
download_enabled <- toupper(Sys.getenv("DOWNLOAD_FILES", "TRUE")) != "FALSE"

dir.create("data/rna_seq", recursive = TRUE, showWarnings = FALSE)
dir.create("data/atac_seq", recursive = TRUE, showWarnings = FALSE)
dir.create("data/manifests", recursive = TRUE, showWarnings = FALSE)

normalize_cell_name <- function(x) {
  x <- trimws(x)
  x <- gsub("CFU-E", "CFUE", x, fixed = TRUE)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

parse_header <- function(line) {
  m <- regexec("^(.+) \\((.+)\\)$", line)
  parts <- regmatches(line, m)[[1]]
  if (length(parts) != 3) {
    stop(sprintf("Could not parse section header: %s", line))
  }
  list(
    cell_type = normalize_cell_name(parts[2]),
    assay = trimws(parts[3])
  )
}

infer_extension <- function(assay) {
  if (assay == "ScriptSeq") {
    return("tsv")
  }
  if (assay == "ATAC-seq") {
    return("bigBed")
  }
  stop(sprintf("Unsupported assay for download: %s", assay))
}

is_supported_assay <- function(assay) {
  assay %in% c("ScriptSeq", "ATAC-seq")
}

infer_output_dir <- function(assay) {
  if (assay == "ScriptSeq") {
    return("data/rna_seq")
  }
  if (assay == "ATAC-seq") {
    return("data/atac_seq")
  }
  stop(sprintf("Unsupported assay for output directory: %s", assay))
}

lines <- readLines(id_file, warn = FALSE)
lines <- trimws(lines)
lines <- lines[nzchar(lines)]

records <- list()
current_section <- NULL
replicate_index <- 0

for (line in lines) {
  if (grepl("^ENCFF", line)) {
    if (is.null(current_section)) {
      stop(sprintf("Found file IDs before any section header: %s", line))
    }
    if (!is_supported_assay(current_section$assay)) {
      next
    }

    ids <- strsplit(line, "\\s+")[[1]]
    if (length(ids) < 2) {
      stop(sprintf("Expected two ENCODE IDs in row: %s", line))
    }

    replicate_index <- replicate_index + 1
    processed_id <- ids[2]
    extension <- infer_extension(current_section$assay)
    output_dir <- infer_output_dir(current_section$assay)
    output_file <- file.path(
      output_dir,
      sprintf(
        "%s_rep%d_%s.%s",
        current_section$cell_type,
        replicate_index,
        processed_id,
        extension
      )
    )

    records[[length(records) + 1]] <- data.frame(
      cell_type = current_section$cell_type,
      assay = current_section$assay,
      replicate = replicate_index,
      file_id = processed_id,
      extension = extension,
      destination = output_file,
      stringsAsFactors = FALSE
    )
  } else {
    current_section <- parse_header(line)
    replicate_index <- 0
  }
}

manifest <- do.call(rbind, records)
manifest$url <- sprintf(
  "%s/%s/@@download/%s.%s",
  base_url,
  manifest$file_id,
  manifest$file_id,
  manifest$extension
)

write.csv(
  manifest,
  "data/manifests/encode_download_manifest.csv",
  row.names = FALSE
)

if (!download_enabled) {
  message("Manifest written to data/manifests/encode_download_manifest.csv")
  message("Dry run enabled; skipping downloads.")
  quit(save = "no", status = 0)
}

for (i in seq_len(nrow(manifest))) {
  url <- manifest$url[i]
  dest <- manifest$destination[i]

  if (file.exists(dest)) {
    message("Skipping existing file: ", dest)
    next
  }

  message("Downloading ", url)
  download.file(url, destfile = dest, mode = "wb")
}
