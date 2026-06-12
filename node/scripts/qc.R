# qc.R — Quality control for expression matrices in geo-microarray-processing
#
# Observation-only: computes metrics, flags issues, returns a decision.
# Never transforms data. The 'clean' subcommand handles transformations.
#
# Decision order (first match wins):
#   veto > rerun > caution > pass

#' Compute QC metrics from an expression matrix CSV file
#'
#' Reads a CSV, computes 7 metrics, and detects expression scale.
#'
#' @param input_path Path to CSV file (genes as rows, samples as columns)
#' @return List with n_samples, n_genes, missing_pct, zero_var_pct,
#'         has_extreme_values, min_correlation, expr_scale, expr_range
compute_qc_metrics <- function(input_path) {
  if (!file.exists(input_path)) {
    stop("File not found: ", input_path)
  }

  # Read CSV with first column as row names
  expr <- as.matrix(read.csv(input_path, row.names = 1, check.names = FALSE))

  n_genes   <- nrow(expr)
  n_samples <- ncol(expr)

  # Missing values
  total_cells <- n_genes * n_samples
  n_missing   <- sum(is.na(expr))
  missing_pct <- if (total_cells > 0) (n_missing / total_cells) * 100 else 0

  # Zero-variance genes
  row_vars     <- apply(expr, 1, var, na.rm = TRUE)
  zero_var_pct <- (sum(row_vars == 0 | is.na(row_vars)) / n_genes) * 100

  # Extreme values
  has_extreme_values <- any(abs(expr) > 1e100, na.rm = TRUE)

  # Pairwise sample correlation (min)
  min_correlation <- NA_real_
  if (n_samples >= 2 && n_genes > 1) {
    cor_matrix      <- suppressWarnings(cor(expr, use = "pairwise.complete.obs"))
    min_correlation <- min(cor_matrix[lower.tri(cor_matrix)], na.rm = TRUE)
  }

  # Expression scale (reuse normalize.R's detect_expr_type)
  expr_scale <- detect_expr_type(as.vector(expr))

  # Expression range
  expr_range <- if (n_genes > 0) {
    max(expr, na.rm = TRUE) - min(expr, na.rm = TRUE)
  } else NA_real_

  list(
    n_samples          = n_samples,
    n_genes            = n_genes,
    missing_pct        = round(missing_pct, 2),
    zero_var_pct       = round(zero_var_pct, 2),
    has_extreme_values = has_extreme_values,
    min_correlation    = if (is.na(min_correlation)) NA_real_ else round(min_correlation, 4),
    expr_scale         = expr_scale,
    expr_range         = round(expr_range, 2)
  )
}

#' Flag QC metrics against thresholds
#'
#' Each metric has a warning threshold and optionally a hard threshold.
#' Returns a list of flags, each with metric, value, threshold, and level.
#'
#' @param metrics List from compute_qc_metrics()
#' @return List of flag lists
flag_qc_metrics <- function(metrics) {
  flags <- list()

  # missing_pct: >10% warning, >50% veto
  if (metrics$missing_pct > 10) {
    level <- if (metrics$missing_pct > 50) "veto" else "warning"
    flags <- c(flags, list(list(
      metric = "missing_pct", value = metrics$missing_pct,
      threshold = if (level == "veto") 50 else 10, level = level
    )))
  }

  # zero_var_pct: >20% caution
  if (metrics$zero_var_pct > 20) {
    flags <- c(flags, list(list(
      metric = "zero_var_pct", value = metrics$zero_var_pct,
      threshold = 20, level = "caution"
    )))
  }

  # extreme_values: any → caution
  if (isTRUE(metrics$has_extreme_values)) {
    flags <- c(flags, list(list(
      metric = "extreme_values", value = TRUE,
      threshold = 1e100, level = "caution"
    )))
  }

  # n_samples: <6 → caution
  if (metrics$n_samples < 6) {
    flags <- c(flags, list(list(
      metric = "n_samples", value = metrics$n_samples,
      threshold = 6, level = "caution",
      note = "Low sample count may reduce statistical power"
    )))
  }

  # n_genes: <5000 → rerun, <100 → veto
  if (metrics$n_genes < 5000) {
    level <- if (metrics$n_genes < 100) "veto" else "rerun"
    flags <- c(flags, list(list(
      metric = "n_genes", value = metrics$n_genes,
      threshold = if (level == "veto") 100 else 5000, level = level
    )))
  }

  # min_correlation: <0.8 → warning, <0.5 → rerun
  if (!is.na(metrics$min_correlation) && metrics$min_correlation < 0.8) {
    level <- if (metrics$min_correlation < 0.5) "rerun" else "warning"
    flags <- c(flags, list(list(
      metric = "min_correlation", value = metrics$min_correlation,
      threshold = if (level == "rerun") 0.5 else 0.8, level = level
    )))
  }

  # expr_scale: not log2 → caution (recommend clean)
  if (metrics$expr_scale != "log") {
    flags <- c(flags, list(list(
      metric = "expr_scale", value = metrics$expr_scale,
      threshold = "log", level = "caution",
      note = sprintf("Expression is %s-scale; recommend 'clean' subcommand", metrics$expr_scale)
    )))
  }

  flags
}

#' Determine QC decision from flags
#'
#' Priority: veto > rerun > caution > pass
#'
#' @param flags List of flag lists from flag_qc_metrics()
#' @return Character: "pass", "caution", "rerun", or "veto"
qc_decision <- function(flags) {
  levels <- vapply(flags, `[[`, "", "level")

  if ("veto" %in% levels) return("veto")
  if ("rerun" %in% levels) return("rerun")
  if (length(flags) > 0) return("caution")
  "pass"
}

#' Run QC on an expression matrix CSV file
#'
#' Main entry point for the QC subcommand. Loads, computes, flags, decides.
#'
#' @param input_path Path to expression matrix CSV file
#' @return List with status, decision, metrics, flags
run_qc <- function(input_path) {
  if (!file.exists(input_path)) {
    report_and_classify(paste("File not found:", input_path))
    return(list(status = "error", msg = paste("File not found:", input_path)))
  }

  result <- tryCatch(
    compute_qc_metrics(input_path),
    error = function(e) {
      return(list(status = "error", msg = e$message))
    }
  )

  if (!is.null(result$status) && result$status == "error") {
    return(result)
  }

  flags    <- flag_qc_metrics(result)
  decision <- qc_decision(flags)

  list(
    status   = "qc_complete",
    decision = decision,
    metrics  = result,
    flags    = flags
  )
}
