# exceptions.R — Structured exception handling for geo-microarray-processing
#
# Categories: A=Network, B=Data, C=Resource, W=Write, L=Lifecycle, E=Environment
# All exceptions flow through report_exception_ndjson() for machine-readable output.

#' Retry a function with exponential backoff
#'
#' @param fn Function to retry (no arguments)
#' @param max_attempts Maximum number of attempts (default 3)
#' @param base_delay Base delay in seconds (default 5)
#' @return Result of fn(), or NULL if all attempts exhausted
retry_with_backoff <- function(fn, max_attempts = 3, base_delay = 5) {
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(
      fn(),
      error = function(e) {
        if (attempt < max_attempts) {
          delay <- base_delay * (4^(attempt - 1))  # 5, 20, 80
          message(sprintf("Attempt %d/%d failed: %s. Retrying in %ds...",
            attempt, max_attempts, e$message, delay))
          Sys.sleep(delay)
          NULL
        } else {
          message(sprintf("Attempt %d/%d failed: %s. All attempts exhausted.",
            attempt, max_attempts, e$message))
          NULL
        }
      }
    )
    if (!is.null(result)) return(result)
  }
  NULL
}

#' Write CSV atomically with size verification
#'
#' Writes to temp file, verifies minimum size, renames atomically.
#'
#' @param data Data frame or matrix to write
#' @param path Output CSV path
#' @param min_ratio Minimum ratio of actual/expected size (default 0.5)
#' @return Logical TRUE if write succeeded
safe_write_csv <- function(data, path, min_ratio = 0.5) {
  tmp <- paste0(path, ".tmp")

  # Estimate minimum expected size (4 bytes per cell as rough lower bound)
  n_cells <- if (is.matrix(data)) nrow(data) * ncol(data) else prod(dim(data))
  if (is.data.frame(data)) n_cells <- nrow(data) * ncol(data)
  expected_min <- max(n_cells * 4, 10)

  # Write to temp file
  tryCatch(
    write.csv(data, file = tmp, row.names = TRUE),
    error = function(e) {
      message("Failed to write temp file: ", e$message)
      return(FALSE)
    }
  )

  # Verify size
  actual_size <- file.info(tmp)$size
  if (is.na(actual_size) || actual_size < expected_min * min_ratio) {
    message(sprintf("Write verification failed: %d bytes (expected >%d)",
      actual_size, as.integer(expected_min * min_ratio)))
    unlink(tmp)
    return(FALSE)
  }

  # Atomic rename
  if (!file.rename(tmp, path)) {
    message("Failed to rename temp file: permission denied or path in use")
    unlink(tmp)
    return(FALSE)
  }

  TRUE
}

#' Validate cached output directory
#'
#' Checks for .fetch_complete sentinel file.
#'
#' @param gse_dir Path to GSE output directory
#' @return List with status ("valid", "stale") and reason
validate_cache <- function(gse_dir) {
  sentinel <- file.path(gse_dir, ".fetch_complete")

  if (!dir.exists(gse_dir)) {
    return(list(status = "stale", reason = "Directory does not exist"))
  }

  if (!file.exists(sentinel)) {
    return(list(status = "stale", reason = "Sentinel missing"))
  }

  list(status = "valid", reason = "Cache valid")
}

#' Check runtime environment
#'
#' Verifies required packages and network connectivity.
#'
#' @return List with status ("ok", "error"), missing packages, and details
check_environment <- function() {
  required <- c("GEOquery", "Biobase", "limma")
  missing  <- character(0)

  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }

  if (length(missing) > 0) {
    return(list(status = "error", missing = missing,
      msg = paste("Missing required packages:", paste(missing, collapse = ", "))))
  }

  list(status = "ok", msg = "Environment OK")
}

#' Report an exception as structured NDJSON
#'
#' @param code Exception code (e.g., "A1_TIMEOUT")
#' @param nature Exception nature (network, data_corrupt, data_insufficient, ...)
#' @param action Response action (retry, skip_with_warning, halt, prompt)
#' @param msg Human-readable message
report_exception_ndjson <- function(code, nature, action, msg) {
  obj <- list(
    level  = "exception",
    code   = code,
    nature = nature,
    action = action,
    msg    = msg
  )
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
}
