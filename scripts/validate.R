# validate.R — Validation functions for geo-microarray-processing
#
# Provides data integrity checks for expression matrices, gene-level
# expression data, and CEL files.

EXTREME_VALUE_THRESHOLD <- 1e50

#' Validate probe-level expression matrix
#'
#' Checks dimensions, data type, and extreme value thresholds.
#'
#' @param expr_matrix Numeric matrix (probes x samples)
#' @return List with 'valid' (logical), 'reason' (character), 'n_rows', 'n_cols'
validate_expr_matrix <- function(expr_matrix) {
  if (is.null(expr_matrix)) {
    return(list(valid = FALSE, reason = "Expression matrix is NULL"))
  }

  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    return(list(valid = FALSE, reason = "Expression data must be matrix or data frame"))
  }

  n_rows <- nrow(expr_matrix)
  n_cols <- ncol(expr_matrix)

  if (n_rows == 0) {
    return(list(valid = FALSE, reason = "Expression matrix has 0 rows"))
  }
  if (n_cols == 0) {
    return(list(valid = FALSE, reason = "Expression matrix has 0 columns"))
  }

  # Check for extreme values
  if (is.matrix(expr_matrix)) {
    max_val <- max(abs(expr_matrix), na.rm = TRUE)
  } else {
    numeric_cols <- vapply(expr_matrix, is.numeric, logical(1))
    max_val <- max(abs(as.matrix(expr_matrix[, numeric_cols, drop = FALSE])), na.rm = TRUE)
  }

  if (max_val > EXTREME_VALUE_THRESHOLD) {
    return(list(valid = FALSE,
      reason = sprintf("Extreme values detected (max abs = %.1e > threshold %.1e)",
        max_val, EXTREME_VALUE_THRESHOLD)))
  }

  list(valid = TRUE, reason = "OK", n_rows = n_rows, n_cols = n_cols)
}

#' Validate gene-level expression data frame
#'
#' @param expr_gene Data frame with gene symbols and expression columns
#' @return Logical TRUE if valid
validate_gene_expression <- function(expr_gene) {
  if (is.null(expr_gene)) return(FALSE)
  if (!is.data.frame(expr_gene)) return(FALSE)
  if (nrow(expr_gene) == 0 || ncol(expr_gene) < 2) return(FALSE)
  TRUE
}

#' Validate CEL file integrity
#'
#' Checks if a CEL file exists and is non-empty.
#'
#' @param cel_file Path to the CEL file
#' @return List with 'valid' (logical) and 'reason' (character)
validate_cel_integrity <- function(cel_file) {
  if (!file.exists(cel_file)) {
    return(list(valid = FALSE, reason = "File not found"))
  }

  file_size <- file.info(cel_file)$size
  if (is.na(file_size) || file_size == 0) {
    return(list(valid = FALSE, reason = "File is empty or corrupted"))
  }

  list(valid = TRUE, reason = "OK")
}
