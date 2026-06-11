# pipeline.R — Metadata-driven preprocessing pipeline detection
#
# Parses GEO data_processing text to identify the submitter's preprocessing
# pipeline, then applies the platform-appropriate normalization transform.
#
# Pipeline detection keywords are matched case-insensitively in priority
# order. The first match wins.

#' Detect preprocessing pipeline from GEO data_processing text
#'
#' @param txt Character string from pData data_processing column(s)
#' @return Character pipeline identifier: "rma", "mas5_gcos", "gcrma",
#'         "sst_rma", "genomestudio_avg", "neqc", "lumi", "agilent_fe",
#'         or "unknown"
detect_pipeline <- function(txt) {
  if (is.null(txt) || nchar(txt) == 0) return("unknown")

  txt <- tolower(paste(as.character(txt), collapse = " "))

  # Priority order: specific methods before generic ones
  # neqc must be checked before GenomeStudio (both may appear)
  if (grepl("neqc|normexp.*offset|offset.*16.*normexp", txt)) return("neqc")
  if (grepl("lumi|vst|variance.stabiliz", txt)) return("lumi")
  if (grepl("gcrma|gc.content.*rma|gc-rma", txt)) return("gcrma")
  if (grepl("sst-rma|tac.*software|transcriptome.analysis.console|signal.space", txt)) return("sst_rma")
  if (grepl("\\brma\\b|robust.multi.array|robust.multichip", txt)) return("rma")
  if (grepl("gcos.*target|mas5|target.intensity.*scal|scaled.*target.intensity", txt)) return("mas5_gcos")
  if (grepl("feature.extraction.*agilent|agilent.*feature.extraction|agilent.*loess", txt)) return("agilent_fe")
  if (grepl("genome.?studio|average.normalization.*genome|beadstudio", txt)) return("genomestudio_avg")

  "unknown"
}

#' Detect whether an expression matrix has been quantile-normalized
#'
#' Checks if sample distributions are identical using a 5-percentile
#' coefficient of variation test. After QN, all percentiles are equal
#' across samples, making CV approach zero.
#'
#' @param expr Numeric matrix (probes x samples)
#' @param tol CV threshold (default 0.002)
#' @return Logical TRUE if QN detected
is_quantile_normalized <- function(expr, tol = 0.002) {
  if (ncol(expr) < 2) return(FALSE)

  pcts <- c(0.25, 0.50, 0.75, 0.90, 0.95)
  vals <- apply(expr, 2, quantile, probs = pcts, na.rm = TRUE)
  cvs  <- apply(vals, 1, function(r) sd(r, na.rm = TRUE) / abs(mean(r, na.rm = TRUE)))

  all(cvs < tol)
}

#' Apply pipeline-appropriate normalization transform
#'
#' @param expr Numeric matrix (probes x samples)
#' @param pipeline Character pipeline identifier from detect_pipeline()
#' @return List with expr (transformed matrix) and transform (description)
apply_pipeline_transform <- function(expr, pipeline) {
  if (pipeline %in% c("rma", "gcrma", "neqc", "lumi", "sst_rma", "agilent_fe")) {
    list(expr = expr, transform = "none")
  } else if (pipeline == "mas5_gcos") {
    list(expr = log2(expr + 1e-6), transform = "log2(x+1e-6)")
  } else if (pipeline == "genomestudio_avg") {
    if (min(expr, na.rm = TRUE) < 0) {
      shift <- min(expr, na.rm = TRUE)
      expr  <- expr - shift
      list(expr = log2(expr + 1e-6), transform = "shift+log2(x+1e-6)")
    } else {
      list(expr = log2(expr + 1e-6), transform = "log2(x+1e-6)")
    }
  } else {
    # unknown — no transform, warn downstream
    list(expr = expr, transform = "none", warning = "Unknown pipeline; no transform applied")
  }
}
