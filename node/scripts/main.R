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

# Resolve script directory for relative sourcing (works from any CWD)
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))

# Source modules in dependency order (fetch.R depends on all others)
source(file.path(script_dir, "report.R"))
source(file.path(script_dir, "normalize.R"))
source(file.path(script_dir, "validate.R"))
source(file.path(script_dir, "species.R"))
source(file.path(script_dir, "annotate.R"))
source(file.path(script_dir, "pipeline.R"))
source(file.path(script_dir, "raw.R"))
source(file.path(script_dir, "exceptions.R"))
source(file.path(script_dir, "fetch.R"))
source(file.path(script_dir, "qc.R"))
source(file.path(script_dir, "clean.R"))

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
    output          = NULL,
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
    } else if (key == "--output") {
      i <- i + 1
      if (i <= length(remaining)) opts$output <- remaining[i]
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
    } else if (startsWith(key, "--output=")) {
      opts$output <- sub("^--output=", "", key)
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

#' Fetch subcommand
do_fetch <- function(opts) {
  report_info(sprintf("Fetching GEO data for %s...", opts$gse_id))
  result <- fetch_geo_data(opts)

  if (result$status == "error") {
    err_msg <- paste(unlist(result$errors), collapse = "; ")
    report_error(err_msg)
    return(invisible(NULL))
  }

  # Build file list for NDJSON report
  files <- list()
  if (!is.null(result$probe_file)) {
    for (f in result$probe_file) {
      files <- c(files, list(list(path = f)))
    }
  }
  if (!is.null(result$gene_file)) {
    for (f in result$gene_file) {
      files <- c(files, list(list(path = f)))
    }
  }
  if (!is.null(result$meta_file)) {
    for (f in result$meta_file) {
      files <- c(files, list(list(path = f)))
    }
  }

  if (length(result$warnings) > 0) {
    for (w in result$warnings) {
      report_info(sprintf("Warning: %s", w))
    }
  }

  report_result(result$status, files = files, metadata = result$metadata)
}

#' QC subcommand
do_qc <- function(opts) {
  report_info(sprintf("Running QC on %s...", opts$input))
  result <- run_qc(opts$input)

  if (result$status == "error") {
    report_error(result$msg)
    return(invisible(NULL))
  }

  for (flag in result$flags) {
    note <- if (!is.null(flag$note)) paste0(" (", flag$note, ")") else ""
    msg  <- sprintf("QC flag [%s]: %s = %s (threshold: %s)%s",
                    toupper(flag$level), flag$metric, flag$value, flag$threshold, note)
    report_info(msg)
  }

  report_result(result$status, decision = result$decision,
                metrics = result$metrics, flags = result$flags)
}

#' Clean subcommand
do_clean <- function(opts) {
  report_info(sprintf("Cleaning %s...", opts$input))
  output <- if (!is.null(opts$output)) opts$output else sub("\\.csv$", "_clean.csv", opts$input)
  result <- run_clean(opts$input, output)

  if (result$status == "error") {
    report_error(result$msg)
    return(invisible(NULL))
  }

  report_info(sprintf("Input scale: %s, applied: %s, output scale: %s",
    result$input_scale, result$applied_transform, result$output_scale))

  report_result(result$status,
    files = list(list(path = result$output_path, rows = result$n_rows, cols = result$n_cols)),
    metadata = list(
      input_scale       = result$input_scale,
      output_scale      = result$output_scale,
      applied_transform = result$applied_transform
    ))
}

# -------------------------------------------------------------------
# Main dispatch
# -------------------------------------------------------------------

main <- function() {
  opts <- parse_args()

  # Environment check (after arg parsing — help text shouldn't require packages)
  env <- check_environment()
  if (env$status == "error") {
    report_exception_ndjson("E801_ENV_PKG", "env_bug", "halt", env$msg)
    quit(status = 3)
  }
  for (w in env$warnings) {
    report_info(sprintf("Env warning: %s", w))
  }

  switch(opts$subcommand,
    fetch = do_fetch(opts),
    qc    = do_qc(opts),
    clean = do_clean(opts)
  )
}

if (sys.nframe() == 0) {
  main()
}
