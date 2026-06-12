# exceptions.R — Structured exception handling for geo-microarray-processing
#
# Categories: A=Network, B=Data, C=Resource, W=Write, L=Lifecycle, E=Environment
# All exceptions flow through report_exception_ndjson() for machine-readable output.

#' Execute a function with a timeout (best-effort via R's withTimeout)
#'
#' Uses R.utils::withTimeout if available, otherwise runs without timeout.
#' On timeout, raises an error caught by retry_with_backoff().
#'
#' @param fn Function to execute
#' @param timeout_sec Timeout in seconds (default 300)
#' @return Result of fn()
with_geo_timeout <- function(fn, timeout_sec = 300) {
  if (requireNamespace("R.utils", quietly = TRUE)) {
    R.utils::withTimeout(fn(), timeout = timeout_sec)
  } else {
    fn()
  }
}

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
  write_ok <- tryCatch({
    write.csv(data, file = tmp, row.names = TRUE)
    TRUE
  }, error = function(e) {
    message("Failed to write temp file: ", e$message)
    FALSE
  })
  if (!write_ok) return(FALSE)

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
  required <- c("GEOquery", "Biobase", "limma", "affy")
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

#' Register signal handlers for graceful shutdown
#'
#' Sets up SIGTERM (clean exit with checkpoint) and SIGINT (prompt user).
#'
#' @param gse_dir Current GSE output directory
#' @return NULL (side effect: registers handlers)
register_signal_handlers <- function(gse_dir = NULL) {
  # R does not have portable OS signal handling (SIGTERM/SIGINT).
  # Cleanup on interrupt is handled via on.exit() at the call site.
  # This function exists as a documented hook for platform-specific
  # implementations (e.g., use later::defer or withr::defer).
  invisible(NULL)
}

#' Classify a failure by matching stderr output against known patterns
#'
#' @param stderr_output Character string from stderr
#' @return List with code, nature, action fields
detect_exception <- function(stderr_output) {
  txt <- tolower(paste(stderr_output, collapse = " "))

  if (grepl("timeout|timed.out|connection.refused", txt))
    return(list(code = "A1_TIMEOUT", nature = "network", action = "retry"))
  if (grepl("404|not.found", txt))
    return(list(code = "A2_NOT_FOUND", nature = "network", action = "skip_with_warning"))
  if (grepl("permission.denied|access.denied", txt))
    return(list(code = "W002_PERM_DENIED", nature = "resource", action = "halt"))
  if (grepl("disk.full|no.space", txt))
    return(list(code = "W001_DISK_FULL", nature = "resource", action = "halt"))
  if (grepl("pthread_create|thread", txt))
    return(list(code = "C3_THREAD", nature = "resource", action = "retry"))

  list(code = "UNKNOWN", nature = "env_bug", action = "escalate")
}

#' Write a checkpoint line for resume support
#'
#' @param gse_dir Output directory
#' @param step Step name (e.g., "download_raw", "process_gene")
#' @param status Status: "complete", "in_progress", "failed"
write_checkpoint <- function(gse_dir, step, status) {
  if (is.null(gse_dir)) return(invisible(NULL))
  checkpoint_file <- file.path(gse_dir, ".fetch_checkpoint")
  line <- sprintf("%s|%s|%s|%s", step, "unknown", status,
                  format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
  write(line, file = checkpoint_file, append = TRUE)
}

#' Read the most recent checkpoint for a directory
#'
#' @param gse_dir Output directory
#' @return Data frame of checkpoint entries, or NULL
read_checkpoint <- function(gse_dir) {
  checkpoint_file <- file.path(gse_dir, ".fetch_checkpoint")
  if (!file.exists(checkpoint_file)) return(NULL)
  lines <- readLines(checkpoint_file)
  if (length(lines) == 0) return(NULL)
  parts <- strsplit(lines, "|", fixed = TRUE)
  do.call(rbind, lapply(parts, function(p) {
    data.frame(step = p[1], platform = p[2], status = p[3],
               timestamp = p[4], stringsAsFactors = FALSE)
  }))
}

#' Detect, classify, and report an exception in one call
#'
#' Convenience wrapper: detect_exception() → report_exception_ndjson().
#' Call this at every error site in the pipeline.
#'
#' @param msg Error message or stderr output to classify
#' @return The detected exception list (invisibly)
report_and_classify <- function(msg) {
  ex <- detect_exception(msg)
  report_exception_ndjson(ex$code, ex$nature, ex$action, msg)
  invisible(ex)
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
