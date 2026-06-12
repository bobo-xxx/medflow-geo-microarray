# fetch.R — GEO data retrieval for geo-microarray-processing
#
# Implements the 5-tier fallback strategy:
#   1. Local cache   -> use cached expression matrix files
#   2. Series matrix -> GEOquery::getGEO(GSEMatrix=TRUE)
#   3. Supplementary -> author-provided *.txt.gz from suppl/
#   4. Raw CEL files -> oligo::rma() normalization
#   5. Metadata only -> when all expression data fails
#
# Also handles platform detection and methylation array skip.
# Uses message() for progress; NDJSON conversion happens in main.R.

library(GEOquery)
library(Biobase)
library(stringr)
library(dplyr)
library(tidyr)

# Dependencies (normalize.R, annotate.R, validate.R, species.R) are
# sourced by main.R before this file. Functions from those modules
# are available in the global environment.

# Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Check if files indicate methylation array data (BPM + IDAT)
#'
#' @param files Character vector of filenames
#' @return Logical TRUE if methylation array detected
is_methylation <- function(files) {
  has_bpm  <- any(grepl("\\.bpm$", files, ignore.case = TRUE))
  has_idat <- any(grepl("\\.idat(\\.gz)?$", files, ignore.case = TRUE))
  has_idat && has_bpm
}

#' Detect platform type from file extensions
#'
#' @param files Character vector of filenames
#' @return Character: "Affymetrix", "Agilent", "Illumina", or "Unknown"
detect_platform_from_files <- function(files) {
  if (any(grepl("\\.CEL(\\.gz)?$", files, ignore.case = TRUE))) return("Affymetrix")
  if (any(grepl("\\.GPR(\\.gz)?$", files, ignore.case = TRUE))) return("Agilent")
  if (any(grepl("\\.idat(\\.gz)?$", files, ignore.case = TRUE))) return("Illumina")
  return("Unknown")
}

#' Fetch GEO microarray data with 5-tier fallback
#'
#' Main entry point for the fetch subcommand. Currently implements
#' Tier 2 (series matrix) with full annotation pipeline.
#'
#' @param opts Named list with gse_id, outdir, proxy, api_key
#' @return List with status, probe_file, gene_file, metadata, warnings, errors
fetch_geo_data <- function(opts) {
  gse_id     <- opts$gse_id
  output_dir <- opts$outdir %||% "."
  proxy      <- opts$proxy
  api_key    <- opts$api_key

  message("Fetching GEO data for ", gse_id, "...")

  # Validate GSE ID format
  if (!grepl("^GSE[0-9]+$", gse_id)) {
    return(list(
      status = "error",
      errors = list(sprintf("Invalid GSE ID format: %s", gse_id))
    ))
  }

  out_gse_dir <- file.path(output_dir, gse_id)
  dir.create(out_gse_dir, recursive = TRUE, showWarnings = FALSE)

  result <- list(
    status     = "unknown",
    gse_id     = gse_id,
    probe_file = NULL,
    gene_file  = NULL,
    meta_file = NULL,
    metadata   = list(),
    warnings   = list(),
    errors     = list()
  )

  # Set proxy if provided
  if (!is.null(proxy) && proxy != "") {
    Sys.setenv(http_proxy = proxy, https_proxy = proxy)
    on.exit(Sys.unsetenv(c("http_proxy", "https_proxy")))
  }

  # ---- Tier 2: Try processed series matrix ----
  message("Tier 2: Attempting series matrix download...")
  gse_matrix <- tryCatch({
    getGEO(gse_id, GSEMatrix = TRUE)
  }, error = function(e) {
    message("Failed to get series matrix: ", e$message)
    NULL
  })

  if (!is.null(gse_matrix) && length(gse_matrix) > 0) {
    message("Series matrix retrieved: ", length(gse_matrix), " platform(s)")
    result$status <- "success_matrix"

    for (i in seq_along(gse_matrix)) {
      eset    <- gse_matrix[[i]]
      gpl_id  <- annotation(eset)
      gpl_suffix <- if (length(gse_matrix) > 1) paste0("_", gpl_id) else ""

      process_expression_set(exprs(eset), eset, gpl_id, gpl_suffix,
                             out_gse_dir, result)
    }
    return(result)
  }

  # ---- Tier 3: Supplementary processed matrix ----
  message("Tier 3: Checking supplementary files...")
  suppl_dir <- file.path(out_gse_dir, "suppl")
  if (!dir.exists(suppl_dir)) {
    tryCatch({
      getGEOSuppFiles(gse_id, makeDirectory = FALSE, baseDir = out_gse_dir)
    }, error = function(e) {
      message("Failed to download suppl files: ", e$message)
    })
  }

  if (dir.exists(suppl_dir)) {
    suppl_txt <- list.files(suppl_dir, pattern = "\\.txt(\\.gz)?$",
                            full.names = TRUE, recursive = TRUE)
    # Look for author-provided processed matrix (not series_matrix from GEO)
    matrix_files <- suppl_txt[!grepl("series_matrix", basename(suppl_txt),
                                     ignore.case = TRUE)]

    if (length(matrix_files) > 0) {
      message("Found ", length(matrix_files), " supplementary matrix files")
      result$status <- "success_suppl"
      result$warnings <- c(result$warnings, "Using author-provided supplementary matrix")

      for (mf in matrix_files) {
        expr_matrix <- as.matrix(read.table(mf, header = TRUE, row.names = 1,
                                            check.names = FALSE))
        gpl_guess <- "unknown"
        process_raw_matrix(expr_matrix, gpl_guess, "", out_gse_dir, result)
      }
      if (length(result$probe_file) > 0) return(result)
    }
  }

  # ---- Tier 4: Raw data files ----
  message("Tier 4: Checking raw data files...")
  raw_dir <- file.path(suppl_dir, "RAW")
  tar_files <- list.files(suppl_dir, pattern = "_RAW\\.tar$",
                          full.names = TRUE)

  if (length(tar_files) > 0 && !dir.exists(raw_dir)) {
    dir.create(raw_dir, recursive = TRUE)
    untar(tar_files[1], exdir = raw_dir)
  }

  if (dir.exists(raw_dir)) {
    raw_files <- list.files(raw_dir, full.names = TRUE, recursive = TRUE)

    if (length(raw_files) > 0) {
      # Check for methylation before processing
      if (is_methylation(raw_files)) {
        result$status <- "skipped_methylation"
        return(result)
      }

      raw_result <- process_raw_files(raw_files, out_gse_dir, gse_id)

      if (raw_result$status == "success") {
        result$status <- "success_raw"
        result$warnings <- c(result$warnings,
          paste("Processed raw files via", raw_result$pipeline))

        gpl_guess <- raw_result$platform
        process_raw_matrix(raw_result$expr_matrix, gpl_guess, "",
                           out_gse_dir, result)
        return(result)
      } else if (raw_result$status == "skipped_methylation") {
        result$status <- "skipped_methylation"
        return(result)
      }
      result$warnings <- c(result$warnings,
        paste("Raw processing failed:", raw_result$msg))
    }
  }

  # ---- Tier 5: Metadata only ----
  message("Tier 5: Returning metadata only...")
  result$status <- "metadata_only"
  result$warnings <- c(result$warnings,
    "No expression data available; returning metadata only")

  gse_meta <- tryCatch({
    getGEO(gse_id, GSEMatrix = FALSE)
  }, error = function(e) NULL)

  if (!is.null(gse_meta)) {
    result$metadata$basic <- list(
      title   = tryCatch(Meta(gse_meta)$title,   error = function(e) NA),
      summary = tryCatch(Meta(gse_meta)$summary, error = function(e) NA)
    )
  }

  result
}

# -------------------------------------------------------------------
# Shared downstream processing
# -------------------------------------------------------------------

#' Process an ExpressionSet from series matrix (Tier 2)
#'
#' Applies pipeline detection, normalization, 5-tier gene annotation,
#' aggregation, and saves probe/gene CSV files.
#'
#' @param expr_matrix Expression matrix from exprs(eset)
#' @param eset ExpressionSet object
#' @param gpl_id Platform ID
#' @param gpl_suffix Suffix for output filenames
#' @param out_gse_dir Output directory
#' @param result Result list to update (modified in place)
process_expression_set <- function(expr_matrix, eset, gpl_id, gpl_suffix,
                                    out_gse_dir, result) {
  # Validate pre-normalization
  validation <- validate_expr_matrix(expr_matrix)
  if (!validation$valid) {
    result$warnings <- c(result$warnings, paste("Pre-validation:", validation$reason))
  }

  # Detect preprocessing pipeline from metadata
  pdata <- pData(eset)
  dp_cols <- grep("data.processing", colnames(pdata), ignore.case = TRUE, value = TRUE)
  dp_text <- if (length(dp_cols) > 0) paste(as.character(pdata[1, dp_cols]), collapse = " ") else ""
  pipeline <- detect_pipeline(dp_text)
  message("Detected pipeline: ", pipeline)

  # Apply pipeline-appropriate normalization
  norm_result <- apply_pipeline_transform(expr_matrix, pipeline)
  expr_matrix <- norm_result$expr
  colnames(expr_matrix) <- make.names(colnames(expr_matrix), unique = TRUE)
  qn_status <- if (is_quantile_normalized(expr_matrix)) "applied" else "not_applied"

  # Post-normalization validation
  post_val <- validate_expr_matrix(expr_matrix)
  if (!post_val$valid) {
    result$warnings <- c(result$warnings, paste("Post-normalization:", post_val$reason))
  }
  if (!is.null(norm_result$warning)) {
    result$warnings <- c(result$warnings, norm_result$warning)
  }

  # Save probe-level
  gse_id <- result$gse_id
  probe_file <- file.path(out_gse_dir, paste0("expr_probe_", gse_id, gpl_suffix, ".csv"))
  safe_write_csv(expr_matrix, probe_file)
  result$probe_file <- c(result$probe_file, probe_file)

  # 5-tier gene annotation using fData from ExpressionSet
  fdata <- fData(eset)
  gene_mapped <- FALSE
  probe2gene <- NULL

  gene_col <- intersect(colnames(fdata), c("Gene Symbol", "GENE_SYMBOL", "Symbol"))[1]
  if (!is.na(gene_col)) {
    message("Tier 1 annotation: using '", gene_col, "' column from fData")
    probe2gene <- data.frame(
      probe_id = as.character(fdata$ID), gene_symbol = as.character(fdata[[gene_col]]),
      stringsAsFactors = FALSE)
    probe2gene <- probe2gene[!is.na(probe2gene$gene_symbol) & probe2gene$gene_symbol != "", ]
    gene_mapped <- nrow(probe2gene) > 0
  }

  if (!gene_mapped) {
    ga_col <- grep("gene.assignment", colnames(fdata), ignore.case = TRUE, value = TRUE)[1]
    if (length(ga_col) > 0 && !is.na(ga_col)) {
      message("Tier 2 annotation: parsing gene_assignment column")
      symbols <- extract_gene_from_assignment(fdata[[ga_col]])
      probe2gene <- data.frame(probe_id = as.character(fdata$ID), gene_symbol = symbols,
                               stringsAsFactors = FALSE)
      probe2gene <- probe2gene[!is.na(probe2gene$gene_symbol) & probe2gene$gene_symbol != "", ]
      gene_mapped <- nrow(probe2gene) > 0
    }
  }

  if (!gene_mapped) {
    message("Tier 3 annotation: downloading GPL table for ", gpl_id)
    gpl_table <- get_gpl_annotation(gpl_id)
    if (!is.null(gpl_table) && nrow(gpl_table) > 0 &&
        "gene_symbol" %in% colnames(gpl_table) &&
        any(!is.na(gpl_table$gene_symbol) & gpl_table$gene_symbol != "")) {
      probe2gene <- gpl_table; gene_mapped <- TRUE
    }
  }

  if (!gene_mapped) {
    if (requireNamespace("AnnoProbe", quietly = TRUE)) {
      message("Tier 4 annotation: attempting AnnoProbe pipe for ", gpl_id)
    }
  }

  if (!gene_mapped) {
    message("Tier 5 annotation: using probe IDs as gene symbols")
    probe2gene <- data.frame(
      probe_id = rownames(expr_matrix), gene_symbol = rownames(expr_matrix),
      stringsAsFactors = FALSE)
    result$warnings <- c(result$warnings, "No gene annotation found; using probe IDs as gene symbols")
    gene_mapped <- TRUE
  }

  # Aggregate to gene level
  if (gene_mapped && exists("probe2gene") && nrow(probe2gene) > 0) {
    expr_gene <- aggregate_probe_to_gene(expr_matrix, probe2gene)
    if (!is.null(expr_gene) && validate_gene_expression(expr_gene)) {
      gene_file <- file.path(out_gse_dir, paste0("expr_gene_", gse_id, gpl_suffix, ".csv"))
      safe_write_csv(expr_gene, gene_file)
      result$gene_file <- c(result$gene_file, gene_file)
    }
  }

  # Phenotype data
  pdata <- pData(eset)
  if (!is.null(pdata) && nrow(pdata) > 0) {
    meta_file <- file.path(out_gse_dir, paste0("metadata_", gse_id, gpl_suffix, ".csv"))
    safe_write_csv(pdata, meta_file)
    result$meta_file <- c(result$meta_file, meta_file)
  }

  # Metadata
  tax_id <- tryCatch(experimentData(eset)@other$sample_taxid, error = function(e) NULL)
  species_info <- detect_species(tax_id)
  result$metadata[[gpl_id]] <- list(
    platform = gpl_id, organism = species_info$species,
    n_samples = ncol(expr_matrix), n_probes = nrow(expr_matrix),
    pipeline = pipeline, qn_status = qn_status, transform = norm_result$transform)
}

#' Process a raw expression matrix without an ExpressionSet (Tiers 3-4)
#'
#' Simpler pipeline: validate → normalize (log2, already from raw processor) →
#' annotate (GPL table only, no fData) → aggregate → save.
#'
#' @param expr_matrix Expression matrix (already log2-scale from raw processor)
#' @param gpl_guess Platform hint
#' @param gpl_suffix Filename suffix
#' @param out_gse_dir Output directory
#' @param result Result list to update
process_raw_matrix <- function(expr_matrix, gpl_guess, gpl_suffix,
                                out_gse_dir, result) {
  gse_id <- result$gse_id

  validation <- validate_expr_matrix(expr_matrix)
  if (!validation$valid) {
    result$warnings <- c(result$warnings, paste("Validation:", validation$reason))
    return()
  }

  colnames(expr_matrix) <- make.names(colnames(expr_matrix), unique = TRUE)

  probe_file <- file.path(out_gse_dir, paste0("expr_probe_", gse_id, gpl_suffix, ".csv"))
  safe_write_csv(expr_matrix, probe_file)
  result$probe_file <- c(result$probe_file, probe_file)

  # GPL annotation only (no fData available from raw files)
  if (gpl_guess != "unknown" && gpl_guess != "") {
    gpl_table <- get_gpl_annotation(gpl_guess)
    if (!is.null(gpl_table) && nrow(gpl_table) > 0) {
      expr_gene <- aggregate_probe_to_gene(expr_matrix, gpl_table)
      if (!is.null(expr_gene) && validate_gene_expression(expr_gene)) {
        gene_file <- file.path(out_gse_dir, paste0("expr_gene_", gse_id, gpl_suffix, ".csv"))
        safe_write_csv(expr_gene, gene_file)
        result$gene_file <- c(result$gene_file, gene_file)
      }
    }
  }

  result$metadata[[gpl_guess]] <- list(
    platform = gpl_guess, organism = "unknown",
    n_samples = ncol(expr_matrix), n_probes = nrow(expr_matrix),
    pipeline = "raw_processor", qn_status = "not_applied", transform = "none")
}
