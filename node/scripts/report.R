# report.R — NDJSON reporting helpers for geo-microarray-processing
#
# All output to stdout is valid NDJSON. Each function writes one line.
# Use report_info() for progress, report_result() for final output,
# and report_error() for terminal failures.

# Exception accumulator (filled by report_exception_ndjson during run)
.exceptions <- list()

#' Write an info-level NDJSON message to stdout
#'
#' @param msg Character string with progress message
#' @param ... Additional named fields to include in the JSON object
report_info <- function(msg, ...) {
  extra <- list(...)
  obj <- c(list(level = "info", msg = msg), extra)
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
}

#' Write a result-level NDJSON message to stdout
#'
#' @param status Character status code (success_matrix, error, etc.)
#' @param files List of file info (each with path, rows, cols)
#' @param metadata List of metadata (platform, organism, n_samples)
#' @param ... Additional named fields
report_result <- function(status, files = list(), metadata = list(), ...) {
  extra <- list(...)
  obj <- c(list(level = "result", status = status,
                files = files, metadata = metadata), extra)
  # Remove empty lists for cleaner output
  if (length(files) == 0) obj$files <- NULL
  if (length(metadata) == 0) obj$metadata <- NULL
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
}

#' Write a prompt-level NDJSON message for interactive input
#'
#' @param code Exception code (e.g., "A3_FAILED")
#' @param msg Human-readable prompt message
#' @param options List of option lists, each with choice, label, optional prompt
#' @param timeout_sec Timeout in seconds (default 300)
report_prompt <- function(code, msg, options, timeout_sec = 300) {
  obj <- list(
    level   = "prompt",
    code    = code,
    msg     = msg,
    options = options,
    timeout = timeout_sec
  )
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
  # In non-interactive mode, return default (first option)
  if (!interactive() || !isatty(stdin())) {
    return(if (length(options) > 0) options[[1]] else NULL)
  }
  response <- tryCatch(readLines(stdin(), n = 1), error = function(e) NULL)
  if (!is.null(response) && nchar(response) > 0) {
    return(jsonlite::fromJSON(response))
  }
  if (length(options) > 0) options[[1]] else NULL
}

#' Write run provenance file for agent inspection
#'
#' Creates .run_result.json in the output directory with full provenance:
#' parameters, output metadata, exceptions, warnings, file list, timestamps.
#'
#' @param out_dir Output directory
#' @param result Result list from fetch_geo_data/run_qc/run_clean
#' @param params Parameter list (subcommand, gse_id, outdir, etc.)
#' @param exit_code Integer exit code
#' @param times Character vector of c(started_at, finished_at) ISO8601 timestamps
write_run_result <- function(out_dir, result, params, exit_code, times) {
  # Build output metadata from per-platform results
  output <- list()
  if (!is.null(result$metadata) && length(result$metadata) > 0) {
    for (name in names(result$metadata)) {
      output[[name]] <- result$metadata[[name]]
    }
  }

  # Build file list
  files <- list()
  for (key in c("probe_file", "gene_file", "meta_file")) {
    if (!is.null(result[[key]])) {
      for (f in result[[key]]) {
        files <- c(files, list(list(path = f)))
      }
    }
  }

  # Clean params (remove NULL, internal values)
  clean_params <- list()
  for (k in names(params)) {
    if (!is.null(params[[k]]) && k != "outdir") {
      clean_params[[k]] <- params[[k]]
    }
  }

  obj <- list(
    node         = "geo-microarray-processing",
    subcommand   = if (is.null(params$subcommand)) "unknown" else params$subcommand,
    gse_id       = params$gse_id,
    status       = if (is.null(result$status)) "unknown" else result$status,
    exit_code    = exit_code,
    started_at   = if (is.na(times[1])) NA_character_ else times[1],
    finished_at  = if (is.na(times[2])) NA_character_ else times[2],
    parameters   = clean_params,
    output       = output,
    exceptions   = if (length(.exceptions) > 0) .exceptions else list(),
    warnings     = if (length(result$warnings) > 0) as.list(result$warnings) else list(),
    files        = files
  )

  json_path <- file.path(out_dir, ".run_result.json")
  writeLines(jsonlite::toJSON(obj, auto_unbox = TRUE, pretty = TRUE), json_path)
}

#' Write an error-level NDJSON message and exit
#'
#' @param msg Error message
#' @param exit_code Integer exit code (default: 1)
report_error <- function(msg, exit_code = 1) {
  obj <- list(level = "error", msg = msg)
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
  quit(status = exit_code)
}
