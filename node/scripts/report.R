# report.R — NDJSON reporting helpers for geo-microarray-processing
#
# All output to stdout is valid NDJSON. Each function writes one line.
# Use report_info() for progress, report_result() for final output,
# and report_error() for terminal failures.

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

#' Write an error-level NDJSON message and exit
#'
#' @param msg Error message
#' @param exit_code Integer exit code (default: 1)
report_error <- function(msg, exit_code = 1) {
  obj <- list(level = "error", msg = msg)
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
  quit(status = exit_code)
}
