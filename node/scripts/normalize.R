# normalize.R — Expression matrix normalization for geo-microarray-processing
#
# Detects expression type (raw/centered/log) and applies appropriate
# transformation. Only raw-scale data is transformed; centered and
# log-scale data pass through unchanged.

#' Detect expression type from quantile and mean analysis
#'
#' @param x Numeric vector of expression values
#' @return Character: "raw", "centered", or "log"
detect_expr_type <- function(x) {
  q <- quantile(x, probs = c(0.01, 0.5, 0.99), na.rm = TRUE)
  mean_x <- mean(x, na.rm = TRUE)

  if (q[3] > 100) return("raw")
  if (abs(mean_x) < 0.5 && q[3] < 10 && q[1] > -10) return("centered")
  return("log")
}

#' Normalize raw expression matrix to log2 scale
#'
#' Applies log2(x + 1e-6) to raw data. Centered and log-scale data
#' pass through unchanged — centered data is already log2-transformed
#' and negative values carry meaningful biological information.
#'
#' @param x Numeric matrix (probes x samples)
#' @return Normalized numeric matrix of same dimensions
normalize_expr_matrix <- function(x) {
  type <- detect_expr_type(x)

  if (type == "raw") x <- log2(x + 1e-6)

  x
}
