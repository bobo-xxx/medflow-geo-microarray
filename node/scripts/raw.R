# raw.R — Raw data file processing for geo-microarray-processing
#
# Detects raw file type from filename patterns, routes to platform-specific
# Bioconductor pipelines, and produces log2-scale probe-level expression
# matrices. All processors converge to the same output format.
#
# Supported: Affymetrix CEL, Illumina IDAT, Agilent GPR, Agilent FE TXT,
# NimbleGen PAIR. Methylation (IDAT+BPM) is skipped.

#' Detect raw file type from filename patterns
#'
#' @param files Character vector of filenames
#' @return Character: "affymetrix", "illumina", "agilent_2c", "agilent_1c",
#'         "nimblegen", "methylation", or "unknown"
detect_raw_type <- function(files) {
  # Priority order: check specific patterns first, then generic
  # Methylation before Illumina (both have IDAT)
  if (any(grepl("[.]BPM$", files, ignore.case = TRUE))) return("methylation")
  if (any(grepl("[.]CEL([.]gz)?$", files, ignore.case = TRUE))) return("affymetrix")
  if (any(grepl("[.]idat([.]gz)?$", files, ignore.case = TRUE))) return("illumina")
  if (any(grepl("[.]GPR([.]gz)?$", files, ignore.case = TRUE))) return("agilent_2c")
  if (any(grepl("[.]PAIR([.]gz)?$", files, ignore.case = TRUE))) return("nimblegen")

  # Agilent FE: check TXT file headers for known column names
  txt_files <- grep("[.]txt$", files, value = TRUE, ignore.case = TRUE)
  if (length(txt_files) > 0) {
    header <- tryCatch(
      readLines(txt_files[1], n = 5),
      error = function(e) ""
    )
    if (any(grepl("ProbeName|GeneName|gTotalGeneSignal|gProcessedSignal",
                  header, ignore.case = TRUE))) {
      return("agilent_1c")
    }
  }

  "unknown"
}

#' Detect GPR file source type from header
#'
#' Reads first line of GPR file: ^ATF → GenePix, ^TYPE → Agilent.
#'
#' @param gpr_file Path to a GPR file
#' @return Character: "genepix" or "agilent"
detect_gpr_source <- function(gpr_file) {
  if (length(gpr_file) == 0 || !file.exists(gpr_file[1])) return("genepix")
  first_line <- tryCatch(
    readLines(gpr_file[1], n = 1),
    error = function(e) ""
  )
  if (grepl("^ATF", first_line)) return("genepix")
  if (grepl("^TYPE", first_line)) return("agilent")
  "genepix"  # default
}

#' Process raw microarray files
#'
#' Detects platform from filenames and routes to the appropriate processor.
#' All processors produce a log2-scale probe-level expression matrix.
#'
#' @param files Character vector of raw file paths
#' @param out_dir Output directory for intermediate files
#' @param gse_id GEO series ID (for filename generation)
#' @return List with status, expr_matrix, and metadata
process_raw_files <- function(files, out_dir, gse_id) {
  if (length(files) == 0) {
    return(list(status = "error", msg = "No raw files provided"))
  }

  raw_type <- detect_raw_type(files)
  message("Detected raw type: ", raw_type)

  switch(raw_type,
    "affymetrix" = process_affy(files, out_dir, gse_id),
    "illumina"   = process_illumina(files, out_dir, gse_id),
    "agilent_2c" = process_agilent_2c(files, out_dir, gse_id),
    "agilent_1c" = process_agilent_1c(files, out_dir, gse_id),
    "nimblegen"  = process_nimblegen(files, out_dir, gse_id),
    "methylation" = list(
      status = "skipped_methylation",
      msg = "Methylation array detected (BPM+IDAT), not expression data"
    ),
    list(status = "error",
      msg = paste("Unknown raw file type. Files:", paste(head(basename(files), 5), collapse = ", ")))
  )
}

#' Process Affymetrix CEL files with RMA
#'
#' @param files Paths to CEL files
#' @param out_dir Output directory
#' @param gse_id GEO series ID
#' @return List with expr_matrix and status
process_affy <- function(files, out_dir, gse_id) {
  if (!requireNamespace("oligo", quietly = TRUE)) {
    return(list(status = "error", msg = "Package 'oligo' required for CEL processing"))
  }

  message("Reading ", length(files), " CEL files...")
  raw <- oligo::read.celfiles(files)
  eset <- oligo::rma(raw)
  expr <- Biobase::exprs(eset)

  list(
    status      = "success",
    expr_matrix = expr,
    platform    = "Affymetrix",
    pipeline    = "oligo::rma()"
  )
}

#' Process Illumina IDAT files with neqc
#'
#' @param files Paths to IDAT files (expression, no BPM)
#' @param out_dir Output directory
#' @param gse_id GEO series ID
#' @return List with expr_matrix and status
process_illumina <- function(files, out_dir, gse_id) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(list(status = "error", msg = "Package 'limma' required for IDAT processing"))
  }

  # Filter to actual IDAT files only
  idat_files <- grep("[.]idat([.]gz)?$", files, value = TRUE, ignore.case = TRUE)
  if (length(idat_files) == 0) {
    return(list(status = "error", msg = "No IDAT files found"))
  }

  message("Reading ", length(idat_files), " IDAT files...")
  # neqc handles: normexp bg → offset +16 → QN → log2
  raw <- limma::read.idat(idat_files)
  expr <- limma::neqc(raw)

  list(
    status      = "success",
    expr_matrix = expr,
    platform    = "Illumina",
    pipeline    = "limma::neqc()"
  )
}

#' Process Agilent two-color GPR files
#'
#' @param files Paths to GPR files
#' @param out_dir Output directory
#' @param gse_id GEO series ID
#' @return List with expr_matrix and status
process_agilent_2c <- function(files, out_dir, gse_id) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(list(status = "error", msg = "Package 'limma' required for GPR processing"))
  }

  gpr_files <- grep("[.]GPR([.]gz)?$", files, value = TRUE, ignore.case = TRUE)
  if (length(gpr_files) == 0) {
    return(list(status = "error", msg = "No GPR files found"))
  }

  source_type <- detect_gpr_source(gpr_files)
  message("Reading ", length(gpr_files), " GPR files (source: ", source_type, ")...")

  rg <- limma::read.maimages(gpr_files, source = source_type)
  rg <- limma::backgroundCorrect(rg, method = "normexp")
  rg <- limma::normalizeWithinArrays(rg, method = "loess")
  expr <- limma::normalizeBetweenArrays(rg, method = "quantile")

  # Extract expression matrix (single-channel: use A-values or M-values)
  if (!is.null(rg$A)) {
    expr_matrix <- rg$A
  } else if (!is.null(rg$M)) {
    expr_matrix <- rg$M
  } else {
    expr_matrix <- as.matrix(rg$genes)
  }

  list(
    status      = "success",
    expr_matrix = expr_matrix,
    platform    = "Agilent_2C",
    pipeline    = "limma::normexp+loess+quantile"
  )
}

#' Process Agilent single-color Feature Extraction TXT files
#'
#' @param files Paths to Agilent FE TXT files
#' @param out_dir Output directory
#' @param gse_id GEO series ID
#' @return List with expr_matrix and status
process_agilent_1c <- function(files, out_dir, gse_id) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(list(status = "error", msg = "Package 'limma' required for Agilent FE processing"))
  }

  txt_files <- grep("[.]txt$", files, value = TRUE, ignore.case = TRUE)
  message("Reading ", length(txt_files), " Agilent FE TXT files...")

  rg <- limma::read.maimages(txt_files, source = "agilent", green.only = TRUE)
  expr <- limma::normalizeBetweenArrays(rg, method = "quantile")

  # Extract expression from green channel
  if (!is.null(rg$E)) {
    expr_matrix <- rg$E
  } else {
    expr_matrix <- as.matrix(rg$genes)
  }

  list(
    status      = "success",
    expr_matrix = expr_matrix,
    platform    = "Agilent_1C",
    pipeline    = "limma::read.maimages(source=agilent)+quantile"
  )
}

#' Process NimbleGen PAIR files
#'
#' @param files Paths to PAIR files
#' @param out_dir Output directory
#' @param gse_id GEO series ID
#' @return List with expr_matrix and status
process_nimblegen <- function(files, out_dir, gse_id) {
  if (!requireNamespace("oligo", quietly = TRUE)) {
    return(list(status = "error", msg = "Package 'oligo' required for NimbleGen processing"))
  }

  pair_files <- grep("[.]PAIR([.]gz)?$", files, value = TRUE, ignore.case = TRUE)
  if (length(pair_files) == 0) {
    return(list(status = "error", msg = "No PAIR files found"))
  }

  message("Reading ", length(pair_files), " PAIR files...")
  raw <- oligo::read.xysfiles(pair_files)
  eset <- oligo::rma(raw)
  expr <- Biobase::exprs(eset)

  list(
    status      = "success",
    expr_matrix = expr,
    platform    = "NimbleGen",
    pipeline    = "oligo::read.xysfiles()+rma()"
  )
}
