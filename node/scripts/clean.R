# clean.R — Expression matrix normalization for geo-microarray-processing
#
# Applies normalization transforms to expression matrices:
#   raw      -> log2(x + 1e-6)
#   centered -> pass-through (already log2-scale, negative values meaningful)
#   log      -> pass-through
#
# Reuses normalize_expr_matrix() from normalize.R.

#' Clean/normalize an expression matrix
#'
#' Reads input CSV, detects expression scale, applies normalization,
#' writes cleaned CSV. Reports transform details.
#'
#' @param input  Path to input expression matrix CSV
#' @param output Path for cleaned output CSV (default: input_clean.csv)
#' @return List with status, input_scale, output_scale, applied_transform,
#'         n_rows, n_cols, input_path, output_path
run_clean <- function(input, output = NULL) {
  if (!file.exists(input)) {
    return(list(status = "error", msg = paste("File not found:", input)))
  }

  if (is.null(output)) {
    output <- sub("\\.csv$", "_clean.csv", input)
  }

  # Read
  expr <- as.matrix(read.csv(input, row.names = 1, check.names = FALSE))

  # Detect scale before
  input_scale <- detect_expr_type(as.vector(expr))

  # Apply normalization
  expr_clean <- normalize_expr_matrix(expr)

  # Detect scale after
  output_scale <- detect_expr_type(as.vector(expr_clean))

  # Determine what was applied
  applied_transform <- if (input_scale == "raw") {
    "log2(x+1e-6)"
  } else {
    "none"
  }

  # Write
  write.csv(expr_clean, file = output, row.names = TRUE)

  list(
    status            = "success",
    input_scale       = input_scale,
    output_scale      = output_scale,
    applied_transform = applied_transform,
    n_rows            = nrow(expr_clean),
    n_cols            = ncol(expr_clean),
    input_path        = input,
    output_path       = output
  )
}
