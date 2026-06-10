#!/usr/bin/env Rscript
#
# main.R — Single entry point for geo-microarray-processing node
#
# Usage:
#   Rscript scripts/main.R fetch --gse-id GSE100155 --outdir ./output
#   Rscript scripts/main.R qc    --input ./output/expr_gene.csv
#   Rscript scripts/main.R clean --input ./output/expr_gene.csv
#
# The first positional argument is the subcommand.
# All parameters declared in SKILL.md frontmatter are accepted.
# Output is NDJSON to stdout.

source("scripts/report.R")

# -------------------------------------------------------------------
# Argument parsing
# -------------------------------------------------------------------

#' Parse command-line arguments
#'
#' @param args Character vector of CLI args (default: commandArgs)
#' @return Named list of parsed values
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {

  if (length(args) == 0) {
    cat("Usage: Rscript scripts/main.R <subcommand> [options]\n")
    cat("\nSubcommands:\n")
    cat("  fetch   Download and process GEO microarray data\n")
    cat("  qc      Quality check on expression matrix\n")
    cat("  clean   Clean and normalize expression matrix\n")
    cat("\nOptions:\n")
    cat("  --gse-id ID     GEO Series identifier (fetch)\n")
    cat("  --outdir DIR    Output directory (default: .)\n")
    cat("  --input FILE    Input expression matrix (qc, clean)\n")
    cat("  --proxy URL     HTTP/HTTPS proxy for GEO access\n")
    cat("  --api-key KEY   NCBI API key for higher rate limits\n")
    quit(status = 1)
  }

  subcommand <- args[1]

  valid_subcommands <- c("fetch", "qc", "clean")
  if (!subcommand %in% valid_subcommands) {
    report_error(sprintf("Unknown subcommand '%s'. Valid: %s",
      subcommand, paste(valid_subcommands, collapse = ", ")))
  }

  # Parse --key=value and --key value pairs
  opts <- list(
    subcommand      = subcommand,
    gse_id          = NULL,
    outdir          = ".",
    input           = NULL,
    proxy           = NULL,
    api_key         = {v <- Sys.getenv("NCBI_API_KEY"); if (v == "") NULL else v}
  )

  remaining <- args[-1]
  i <- 1
  while (i <= length(remaining)) {
    key <- remaining[i]

    if (key == "--gse-id") {
      i <- i + 1
      if (i > length(remaining)) {
        report_error("--gse-id requires a value")
      }
      opts$gse_id <- remaining[i]
    } else if (key == "--outdir") {
      i <- i + 1
      if (i <= length(remaining)) opts$outdir <- remaining[i]
    } else if (key == "--input") {
      i <- i + 1
      if (i <= length(remaining)) opts$input <- remaining[i]
    } else if (key == "--proxy") {
      i <- i + 1
      if (i <= length(remaining)) opts$proxy <- remaining[i]
    } else if (key == "--api-key") {
      i <- i + 1
      if (i <= length(remaining)) opts$api_key <- remaining[i]
    } else if (startsWith(key, "--gse-id=")) {
      opts$gse_id <- sub("^--gse-id=", "", key)
    } else if (startsWith(key, "--outdir=")) {
      opts$outdir <- sub("^--outdir=", "", key)
    } else if (startsWith(key, "--input=")) {
      opts$input <- sub("^--input=", "", key)
    } else if (startsWith(key, "--proxy=")) {
      opts$proxy <- sub("^--proxy=", "", key)
    } else if (startsWith(key, "--api-key=")) {
      opts$api_key <- sub("^--api-key=", "", key)
    } else {
      report_error(sprintf("Unknown option: %s", key))
    }
    i <- i + 1
  }

  # Validate required args per subcommand
  if (opts$subcommand == "fetch" && is.null(opts$gse_id)) {
    report_error("fetch subcommand requires --gse-id")
  }
  if (opts$subcommand %in% c("qc", "clean") && is.null(opts$input)) {
    report_error(sprintf("%s subcommand requires --input", opts$subcommand))
  }

  return(opts)
}

# -------------------------------------------------------------------
# Subcommand stubs (filled in by port-fetch-geo, add-qc, add-clean)
# -------------------------------------------------------------------

#' Fetch subcommand — stub
do_fetch <- function(opts) {
  report_info(sprintf("Fetching GEO data for %s...", opts$gse_id))
  report_info("Subcommand 'fetch' not yet implemented")
  report_result("success_matrix", files = list(), metadata = list())
}

#' QC subcommand — stub
do_qc <- function(opts) {
  report_info(sprintf("Running QC on %s...", opts$input))
  report_info("Subcommand 'qc' not yet implemented")
  report_result("pass", files = list(), metadata = list())
}

#' Clean subcommand — stub
do_clean <- function(opts) {
  report_info(sprintf("Cleaning %s...", opts$input))
  report_info("Subcommand 'clean' not yet implemented")
  report_result("success", files = list(), metadata = list())
}

# -------------------------------------------------------------------
# Main dispatch
# -------------------------------------------------------------------

main <- function() {
  opts <- parse_args()

  switch(opts$subcommand,
    fetch = do_fetch(opts),
    qc    = do_qc(opts),
    clean = do_clean(opts)
  )
}

if (sys.nframe() == 0) {
  main()
}
