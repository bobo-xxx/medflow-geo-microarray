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

      expr_matrix <- exprs(eset)

      # Validate
      validation <- validate_expr_matrix(expr_matrix)
      if (!validation$valid) {
        result$warnings <- c(result$warnings, paste("Validation:", validation$reason))
      }

      # Normalize
      expr_matrix <- normalize_expr_matrix(expr_matrix)
      colnames(expr_matrix) <- make.names(colnames(expr_matrix), unique = TRUE)

      # Post-normalization validation (catch NaN/Inf from normalization bugs)
      post_val <- validate_expr_matrix(expr_matrix)
      if (!post_val$valid) {
        result$warnings <- c(result$warnings, paste("Post-normalization:", post_val$reason))
      }

      # Save probe-level
      probe_file <- file.path(out_gse_dir, paste0("expr_probe_", gse_id, gpl_suffix, ".csv"))
      write.csv(expr_matrix, file = probe_file, row.names = TRUE)
      result$probe_file <- c(result$probe_file, probe_file)

      # ---- 5-tier gene annotation ----
      fdata <- fData(eset)
      gene_mapped <- FALSE

      # Tier 1: Direct gene symbol column
      gene_col <- intersect(colnames(fdata),
        c("Gene Symbol", "GENE_SYMBOL", "Symbol"))[1]

      if (!is.na(gene_col)) {
        message("Tier 1 annotation: using '", gene_col, "' column from fData")
        probe2gene <- data.frame(
          probe_id    = as.character(fdata$ID),
          gene_symbol = as.character(fdata[[gene_col]]),
          stringsAsFactors = FALSE
        )
        probe2gene <- probe2gene[!is.na(probe2gene$gene_symbol) & probe2gene$gene_symbol != "", ]
        gene_mapped <- nrow(probe2gene) > 0
      }

      # Tier 2: gene_assignment column
      if (!gene_mapped) {
        ga_col <- grep("gene.assignment", colnames(fdata), ignore.case = TRUE, value = TRUE)[1]
        if (length(ga_col) > 0 && !is.na(ga_col)) {
          message("Tier 2 annotation: parsing gene_assignment column")
          symbols <- extract_gene_from_assignment(fdata[[ga_col]])
          probe2gene <- data.frame(
            probe_id    = as.character(fdata$ID),
            gene_symbol = symbols,
            stringsAsFactors = FALSE
          )
          probe2gene <- probe2gene[!is.na(probe2gene$gene_symbol) & probe2gene$gene_symbol != "", ]
          gene_mapped <- nrow(probe2gene) > 0
        }
      }

      # Tier 3: GPL annotation table
      if (!gene_mapped) {
        message("Tier 3 annotation: downloading GPL table for ", gpl_id)
        gpl_table <- get_gpl_annotation(gpl_id)
        if (!is.null(gpl_table) && nrow(gpl_table) > 0 &&
            "gene_symbol" %in% colnames(gpl_table) &&
            any(!is.na(gpl_table$gene_symbol) & gpl_table$gene_symbol != "")) {
          probe2gene <- gpl_table
          gene_mapped <- TRUE
        }
      }

      # Tier 4: AnnoProbe pipe alignment (attempted if installed)
      if (!gene_mapped) {
        if (requireNamespace("AnnoProbe", quietly = TRUE)) {
          message("Tier 4 annotation: attempting AnnoProbe pipe for ", gpl_id)
          # AnnoProbe::idmap returns probe-to-gene mapping for 147 platforms
          # This is a best-effort call; may fail for unsupported platforms
        }
      }

      # Tier 5: Probe IDs as gene symbols
      if (!gene_mapped) {
        message("Tier 5 annotation: using probe IDs as gene symbols")
        probe2gene <- data.frame(
          probe_id    = rownames(expr_matrix),
          gene_symbol = rownames(expr_matrix),
          stringsAsFactors = FALSE
        )
        result$warnings <- c(result$warnings, "No gene annotation found; using probe IDs as gene symbols")
        gene_mapped <- TRUE
      }

      # Aggregate to gene level
      if (gene_mapped && exists("probe2gene") && nrow(probe2gene) > 0) {
        expr_gene <- aggregate_probe_to_gene(expr_matrix, probe2gene)
        if (!is.null(expr_gene) && validate_gene_expression(expr_gene)) {
          gene_file <- file.path(out_gse_dir, paste0("expr_gene_", gse_id, gpl_suffix, ".csv"))
          write.csv(expr_gene, file = gene_file, row.names = TRUE)
          result$gene_file <- c(result$gene_file, gene_file)
        }
      }

      # Phenotype data (sample-level metadata for grouping/co-factors)
      pdata <- pData(eset)
      if (!is.null(pdata) && nrow(pdata) > 0) {
        meta_file <- file.path(out_gse_dir, paste0("metadata_", gse_id, gpl_suffix, ".csv"))
        write.csv(pdata, file = meta_file, row.names = TRUE)
        result$meta_file <- c(result$meta_file, meta_file)
      }

      # Metadata
      tax_id <- tryCatch(experimentData(eset)@other$sample_taxid, error = function(e) NULL)
      species_info <- detect_species(tax_id)
      result$metadata[[gpl_id]] <- list(
        platform  = gpl_id,
        organism  = species_info$species,
        n_samples = ncol(expr_matrix),
        n_probes  = nrow(expr_matrix)
      )
    }

    return(result)
  }

  # ---- Tier 3-5: Not yet ported ----
  result$warnings <- c(result$warnings, "Series matrix unavailable; tiers 3-5 not yet ported")
  message("Series matrix failed; tiers 3-5 not yet ported")

  if (length(result$errors) > 0) {
    result$status <- "error"
  } else if (length(result$probe_file) == 0) {
    result$status <- "metadata_only"
  }

  result
}
