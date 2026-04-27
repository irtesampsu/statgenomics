suppressPackageStartupMessages({
  library(utils)
})

id_file <- "encode ids.txt"
base_url <- "https://www.encodeproject.org/files"
download_enabled <- toupper(Sys.getenv("DOWNLOAD_FILES", "TRUE")) != "FALSE"

dir.create("data/rna_seq", recursive = TRUE, showWarnings = FALSE)
dir.create("data/atac_seq", recursive = TRUE, showWarnings = FALSE)
dir.create("data/bam/rna_seq", recursive = TRUE, showWarnings = FALSE)
dir.create("data/bam/atac_seq", recursive = TRUE, showWarnings = FALSE)
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

infer_bam_output_dir <- function(assay) {
  if (assay == "ATAC-seq") {
    return("data/bam/atac_seq")
  }
  stop(sprintf("Unsupported assay for BAM output directory: %s", assay))
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
    bam_id <- ids[1]
    processed_id <- ids[2]

    if (current_section$assay == "ATAC-seq") {
      bam_output_dir <- infer_bam_output_dir(current_section$assay)
      bam_output_file <- file.path(
        bam_output_dir,
        sprintf(
          "%s_rep%d_%s.bam",
          current_section$cell_type,
          replicate_index,
          bam_id
        )
      )

      records[[length(records) + 1]] <- data.frame(
        cell_type = current_section$cell_type,
        assay = current_section$assay,
        replicate = replicate_index,
        file_role = "raw_bam",
        file_id = bam_id,
        extension = "bam",
        destination = bam_output_file,
        stringsAsFactors = FALSE
      )
    }

    if (current_section$assay %in% c("ScriptSeq", "ATAC-seq")) {
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
        file_role = "processed",
        file_id = processed_id,
        extension = extension,
        destination = output_file,
        stringsAsFactors = FALSE
      )
    }
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

get_remote_size <- function(url) {
  curl_bin <- Sys.which("curl")
  if (identical(curl_bin, "")) {
    return(NA_real_)
  }

  headers <- tryCatch(
    system2(
      curl_bin,
      args = c("-I", "-L", "--silent", "--show-error", url),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) character(0)
  )
  if (length(headers) == 0) {
    return(NA_real_)
  }

  cl_lines <- headers[grepl("^content-length:", tolower(headers))]
  if (length(cl_lines) == 0) {
    return(NA_real_)
  }

  size <- suppressWarnings(as.numeric(trimws(sub("^[^:]+:", "", tail(cl_lines, 1)))))
  if (is.na(size)) {
    return(NA_real_)
  }
  size
}

download_with_curl_resume <- function(url, dest, max_attempts = 5) {
  curl_bin <- Sys.which("curl")
  if (identical(curl_bin, "")) {
    return(FALSE)
  }

  for (attempt in seq_len(max_attempts)) {
    message(sprintf("Attempt %d/%d: %s", attempt, max_attempts, url))
    status <- system2(
      curl_bin,
      args = c(
        "-L",
        "--fail",
        "--retry", "5",
        "--retry-all-errors",
        "--retry-delay", "5",
        "--continue-at", "-",
        "--output", dest,
        url
      )
    )

    if (status == 0) {
      return(TRUE)
    }

    Sys.sleep(min(30, attempt * 5))
  }

  FALSE
}

download_with_base_r <- function(url, dest, max_attempts = 3) {
  for (attempt in seq_len(max_attempts)) {
    message(sprintf("Fallback attempt %d/%d: %s", attempt, max_attempts, url))
    ok <- tryCatch(
      {
        download.file(url, destfile = dest, mode = "wb")
        TRUE
      },
      error = function(e) FALSE
    )
    if (ok) {
      return(TRUE)
    }
    if (file.exists(dest)) {
      file.remove(dest)
    }
    Sys.sleep(min(20, attempt * 3))
  }
  FALSE
}

for (i in seq_len(nrow(manifest))) {
  url <- manifest$url[i]
  dest <- manifest$destination[i]

  remote_size <- get_remote_size(url)
  if (file.exists(dest)) {
    local_size <- file.info(dest)$size
    if (!is.na(remote_size) && !is.na(local_size) && local_size == remote_size) {
      message("Skipping existing complete file: ", dest)
      next
    }
    message(
      "Existing file is incomplete or size unknown; resuming download: ",
      dest
    )
  }

  message("Downloading ", url)
  ok <- download_with_curl_resume(url, dest)
  if (!ok) {
    message("curl resume download failed, trying base R downloader: ", url)
    ok <- download_with_base_r(url, dest)
  }
  if (!ok) {
    stop("Download failed after retries: ", url)
  }
}
