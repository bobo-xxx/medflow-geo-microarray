# add-qc-subcommand Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `qc` subcommand that loads an expression matrix CSV, computes quality metrics, flags against warning/hard thresholds, and returns a pass/caution/rerun/veto decision via NDJSON. Observation-only — never transforms data.

**Architecture:** `node/scripts/qc.R` is a standalone module sourced by `main.R`. It reads a CSV, runs `detect_expr_type()` from `normalize.R`, computes 7 metrics, flags any that exceed thresholds, and returns a structured result. `main.R`'s `do_qc()` calls it and converts to NDJSON. The `clean` subcommand (next change) will be the one that applies transformations based on QC findings.

**Tech Stack:** R >= 4.3, stats (cor), testthat

---

## File Structure

| File | Purpose |
|---|---|
| `node/scripts/qc.R` | `run_qc(input_path)` — load, compute metrics, flag, decide |
| `tests/testthat/test-qc.R` | Unit tests with synthetic CSV fixtures |
| `node/scripts/main.R` | Replace `do_qc` stub with real implementation |

---

## QC Metrics and Thresholds

| Metric | Warning | Hard | Decision on Hard |
|---|---|---|---|
| missing_pct | > 10% | > 50% | `veto` |
| zero_var_pct | > 20% | — | `caution` |
| extreme_values | any > 1e100 | — | `caution` |
| n_samples | < 6 | — | `caution` |
| n_genes | < 5000 | < 100 | `rerun` |
| min_correlation | < 0.8 | < 0.5 | `rerun` |
| expr_scale | not log2 | — | `caution` (recommend clean) |

Decision rules (first match wins in priority order):
1. Any `veto`-level flag → `veto`
2. Any `rerun`-level flag → `rerun`
3. Any `caution`-level flag → `caution`
4. No flags → `pass`

---

### Task 1: Create tests/testthat/test-qc.R

**Files:**
- Create: `tests/testthat/test-qc.R`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-qc.R`:

```r
library(testthat)

source("../../node/scripts/normalize.R")  # for detect_expr_type
source("../../node/scripts/qc.R")

describe("compute_qc_metrics", {

  it("computes all 7 metrics for valid matrix", {
    set.seed(42)
    m <- matrix(runif(500 * 8, 2, 14), nrow = 500, ncol = 8)
    colnames(m) <- paste0("sample_", 1:8)
    write.csv(m, "test_input.csv", row.names = TRUE)

    result <- compute_qc_metrics("test_input.csv")
    expect_equal(result$n_samples, 8)
    expect_equal(result$n_genes, 500)
    expect_equal(result$missing_pct, 0)
    expect_true(result$expr_range > 0)
    expect_true(!is.na(result$min_correlation))
    unlink("test_input.csv")
  })

  it("detects missing values correctly", {
    m <- matrix(runif(100 * 6, 2, 14), nrow = 100, ncol = 6)
    m[1:30, 1] <- NA  # 30% missing in first column
    colnames(m) <- paste0("sample_", 1:6)
    write.csv(m, "test_input.csv", row.names = TRUE)

    result <- compute_qc_metrics("test_input.csv")
    # 30 NAs out of 600 total = 5%
    expect_gt(result$missing_pct, 4)
    expect_lt(result$missing_pct, 6)
    unlink("test_input.csv")
  })

  it("detects zero-variance genes", {
    m <- matrix(runif(100 * 6, 2, 14), nrow = 100, ncol = 6)
    m[1:25, ] <- 5  # 25% zero-variance
    colnames(m) <- paste0("sample_", 1:6)
    write.csv(m, "test_input.csv", row.names = TRUE)

    result <- compute_qc_metrics("test_input.csv")
    expect_equal(result$zero_var_pct, 25)
    unlink("test_input.csv")
  })

  it("detects extreme values", {
    m <- matrix(runif(100 * 6, 2, 14), nrow = 100, ncol = 6)
    m[1, 1] <- 1e150
    colnames(m) <- paste0("sample_", 1:6)
    write.csv(m, "test_input.csv", row.names = TRUE)

    result <- compute_qc_metrics("test_input.csv")
    expect_true(result$has_extreme_values)
    unlink("test_input.csv")
  })

  it("detects expression scale", {
    # Raw scale
    m <- matrix(runif(100 * 6, 0, 500), nrow = 100, ncol = 6)
    colnames(m) <- paste0("sample_", 1:6)
    write.csv(m, "test_input.csv", row.names = TRUE)
    result <- compute_qc_metrics("test_input.csv")
    expect_equal(result$expr_scale, "raw")
    unlink("test_input.csv")

    # Log scale
    m2 <- matrix(runif(100 * 6, 2, 14), nrow = 100, ncol = 6)
    colnames(m2) <- paste0("sample_", 1:6)
    write.csv(m2, "test_input.csv", row.names = TRUE)
    result2 <- compute_qc_metrics("test_input.csv")
    expect_equal(result2$expr_scale, "log")
    unlink("test_input.csv")
  })
})

describe("flag_qc_metrics", {

  it("returns no flags for clean data", {
    metrics <- list(
      missing_pct = 0, zero_var_pct = 0, has_extreme_values = FALSE,
      n_samples = 12, n_genes = 20000, min_correlation = 0.95,
      expr_scale = "log", expr_range = 12
    )
    flags <- flag_qc_metrics(metrics)
    expect_equal(length(flags), 0)
  })

  it("flags missing_pct > 10% as warning", {
    metrics <- list(
      missing_pct = 15, zero_var_pct = 0, has_extreme_values = FALSE,
      n_samples = 12, n_genes = 20000, min_correlation = 0.95,
      expr_scale = "log", expr_range = 12
    )
    flags <- flag_qc_metrics(metrics)
    miss_flag <- Filter(function(f) f$metric == "missing_pct", flags)
    expect_equal(length(miss_flag), 1)
    expect_equal(miss_flag[[1]]$level, "warning")
  })

  it("flags missing_pct > 50% veto", {
    metrics <- list(
      missing_pct = 55, zero_var_pct = 0, has_extreme_values = FALSE,
      n_samples = 12, n_genes = 20000, min_correlation = 0.95,
      expr_scale = "log", expr_range = 12
    )
    flags <- flag_qc_metrics(metrics)
    miss_flag <- Filter(function(f) f$metric == "missing_pct", flags)
    expect_equal(miss_flag[[1]]$level, "veto")
  })

  it("flags n_genes < 100 veto", {
    metrics <- list(
      missing_pct = 0, zero_var_pct = 0, has_extreme_values = FALSE,
      n_samples = 12, n_genes = 50, min_correlation = 0.95,
      expr_scale = "log", expr_range = 12
    )
    flags <- flag_qc_metrics(metrics)
    gene_flag <- Filter(function(f) f$metric == "n_genes", flags)
    expect_equal(gene_flag[[1]]$level, "veto")
  })

  it("flags min_correlation < 0.5 rerun", {
    metrics <- list(
      missing_pct = 0, zero_var_pct = 0, has_extreme_values = FALSE,
      n_samples = 12, n_genes = 20000, min_correlation = 0.3,
      expr_scale = "log", expr_range = 12
    )
    flags <- flag_qc_metrics(metrics)
    corr_flag <- Filter(function(f) f$metric == "min_correlation", flags)
    expect_equal(corr_flag[[1]]$level, "rerun")
  })

  it("flags n_samples < 6 caution", {
    metrics <- list(
      missing_pct = 0, zero_var_pct = 0, has_extreme_values = FALSE,
      n_samples = 4, n_genes = 20000, min_correlation = 0.95,
      expr_scale = "log", expr_range = 12
    )
    flags <- flag_qc_metrics(metrics)
    sample_flag <- Filter(function(f) f$metric == "n_samples", flags)
    expect_equal(sample_flag[[1]]$level, "caution")
  })
})

describe("qc_decision", {

  it("returns pass when no flags", {
    expect_equal(qc_decision(list()), "pass")
  })

  it("returns veto when any veto flag present", {
    flags <- list(
      list(metric = "missing_pct", level = "warning"),
      list(metric = "n_genes", level = "veto")
    )
    expect_equal(qc_decision(flags), "veto")
  })

  it("returns rerun when rerun flag present (no veto)", {
    flags <- list(
      list(metric = "missing_pct", level = "warning"),
      list(metric = "min_correlation", level = "rerun")
    )
    expect_equal(qc_decision(flags), "rerun")
  })

  it("returns caution for warning-level flags only", {
    flags <- list(
      list(metric = "missing_pct", level = "warning"),
      list(metric = "n_samples", level = "caution")
    )
    expect_equal(qc_decision(flags), "caution")
  })
})

describe("run_qc", {

  it("returns structured QC result with decision", {
    set.seed(42)
    m <- matrix(runif(500 * 8, 2, 14), nrow = 500, ncol = 8)
    colnames(m) <- paste0("sample_", 1:8)
    write.csv(m, "test_qc_input.csv", row.names = TRUE)

    result <- run_qc("test_qc_input.csv")
    expect_equal(result$status, "qc_complete")
    expect_true(result$decision %in% c("pass", "caution", "rerun", "veto"))
    expect_true(is.list(result$metrics))
    expect_equal(result$metrics$n_samples, 8)
    expect_equal(result$metrics$n_genes, 500)
    unlink("test_qc_input.csv")
  })

  it("returns error status for non-existent file", {
    result <- run_qc("/nonexistent/file.csv")
    expect_equal(result$status, "error")
  })
})
```

- [ ] **Step 2: Run test, verify FAILS**

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate geo-microarray-processing
Rscript -e 'testthat::test_file("tests/testthat/test-qc.R")'
```

Expected: FAIL — "cannot open file '../../node/scripts/qc.R'"

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-qc.R
git commit -m "test: add QC subcommand tests (compute_qc_metrics, flag_qc_metrics, qc_decision)

17 tests: 5 compute metrics (missing, zero-var, extreme, scale),
7 flag thresholds (warning/hard across metrics), 2 decision logic,
3 integration (run_qc).
Covers pass/caution/rerun/veto decision chain.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Create node/scripts/qc.R

**Files:**
- Create: `node/scripts/qc.R`

- [ ] **Step 1: Write node/scripts/qc.R**

Create `node/scripts/qc.R`:

```r
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

  n_genes  <- nrow(expr)
  n_samples <- ncol(expr)

  # Missing values
  total_cells <- n_genes * n_samples
  n_missing <- sum(is.na(expr))
  missing_pct <- if (total_cells > 0) (n_missing / total_cells) * 100 else 0

  # Zero-variance genes
  row_vars <- apply(expr, 1, var, na.rm = TRUE)
  zero_var_pct <- (sum(row_vars == 0 | is.na(row_vars)) / n_genes) * 100

  # Extreme values
  has_extreme_values <- any(abs(expr) > 1e100, na.rm = TRUE)

  # Pairwise sample correlation (min)
  min_correlation <- NA_real_
  if (n_samples >= 2 && n_genes > 1) {
    cor_matrix <- suppressWarnings(cor(expr, use = "pairwise.complete.obs"))
    # Lower triangle, excluding diagonal
    min_correlation <- min(cor_matrix[lower.tri(cor_matrix)], na.rm = TRUE)
  }

  # Expression scale (reuse normalize.R's detect_expr_type)
  expr_scale <- detect_expr_type(as.vector(expr))

  # Expression range (log2 values)
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

  # extreme_values: any caution
  if (metrics$has_extreme_values) {
    flags <- c(flags, list(list(
      metric = "extreme_values", value = TRUE,
      threshold = 1e100, level = "caution"
    )))
  }

  # n_samples: <6 caution
  if (metrics$n_samples < 6) {
    flags <- c(flags, list(list(
      metric = "n_samples", value = metrics$n_samples,
      threshold = 6, level = "caution",
      note = "Low sample count may reduce statistical power"
    )))
  }

  # n_genes: <5000 warning, <100 veto
  if (metrics$n_genes < 5000) {
    level <- if (metrics$n_genes < 100) "veto" else "rerun"
    flags <- c(flags, list(list(
      metric = "n_genes", value = metrics$n_genes,
      threshold = if (level == "veto") 100 else 5000, level = level
    )))
  }

  # min_correlation: <0.8 warning, <0.5 rerun
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
      note = sprintf("Expression is %s-scale; recommend 'clean' subcommand for normalization", metrics$expr_scale)
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
    return(list(status = "error", msg = paste("File not found:", input_path)))
  }

  metrics <- tryCatch(
    compute_qc_metrics(input_path),
    error = function(e) {
      return(list(status = "error", msg = e$message))
    }
  )

  if (!is.null(metrics$status) && metrics$status == "error") {
    return(metrics)
  }

  flags    <- flag_qc_metrics(metrics)
  decision <- qc_decision(flags)

  list(
    status   = "qc_complete",
    decision = decision,
    metrics  = metrics,
    flags    = flags
  )
}
```

- [ ] **Step 2: Run test, verify PASS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-qc.R")'
```

Expected: 17/17 PASS

- [ ] **Step 3: Commit**

```bash
git add node/scripts/qc.R
git commit -m "feat: add qc.R with compute_qc_metrics, flag_qc_metrics, qc_decision

Observation-only QC: 7 metrics (missing_pct, zero_var_pct, extreme_values,
n_samples, n_genes, min_correlation, expr_scale) with warning/hard thresholds.
Decision chain: veto > rerun > caution > pass.
17 tests covering all metrics, thresholds, and decision logic.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: Wire qc.R into main.R

**Files:**
- Modify: `node/scripts/main.R`

- [ ] **Step 1: Replace do_qc stub in main.R**

Read current main.R and replace the `do_qc` stub:

```r
#' QC subcommand
do_qc <- function(opts) {
  report_info(sprintf("Running QC on %s...", opts$input))
  result <- run_qc(opts$input)

  if (result$status == "error") {
    report_error(result$msg)
    return(invisible(NULL))
  }

  # Report each flag
  for (flag in result$flags) {
    level <- flag$level
    note  <- if (!is.null(flag$note)) paste0(" (", flag$note, ")") else ""
    msg   <- sprintf("QC flag [%s]: %s = %s (threshold: %s)%s",
                     toupper(level), flag$metric, flag$value, flag$threshold, note)
    report_info(msg)
  }

  report_result(
    status   = result$status,
    decision = result$decision,
    metrics  = result$metrics,
    flags    = result$flags
  )
}
```

Also add `source(file.path(script_dir, "qc.R"))` after the other source lines.

- [ ] **Step 2: Smoke test**

```bash
# Generate a test CSV
Rscript -e '
  set.seed(42)
  m <- matrix(runif(500 * 8, 2, 14), nrow = 500, ncol = 8)
  colnames(m) <- paste0("sample_", 1:8)
  write.csv(m, "test_qc.csv", row.names = TRUE)
'

# Run QC
Rscript node/scripts/main.R qc --input test_qc.csv 2>&1

# Clean up
rm test_qc.csv
```

Expected: NDJSON output with `{"level":"result","status":"qc_complete","decision":"pass",...}`

- [ ] **Step 3: Commit**

```bash
git add node/scripts/main.R
git commit -m "feat: wire qc subcommand to qc.R with NDJSON reporting

Replace do_qc stub with run_qc() call. QC output includes per-flag
NDJSON info lines and a result line with decision + metrics + flags.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: Full test suite verification

- [ ] **Step 1: Run all tests**

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate geo-microarray-processing
Rscript -e 'testthat::test_dir("tests/testthat/")'
```

Expected: 0 FAIL, ~103 tests (86 existing + 17 new)

- [ ] **Step 2: Verify full QC pipeline**

```bash
# Generate test data with known issues
Rscript -e '
  set.seed(42)
  m <- matrix(runif(500 * 8, 2, 14), nrow = 500, ncol = 8)
  m[1:100, 1] <- NA    # 100/4000 = 2.5% missing
  m[101:200, ] <- 5     # 20% zero-variance
  m[1, 2] <- 1e150      # extreme value
  colnames(m) <- paste0("sample_", 1:8)
  write.csv(m, "test_qc_bad.csv", row.names = TRUE)
'
Rscript node/scripts/main.R qc --input test_qc_bad.csv 2>&1
rm test_qc_bad.csv
```

Expected: NDJSON with `decision: "caution"` (missing < 10%, zero_var = 20% at threshold, extreme value detected, scale is log)

- [ ] **Step 3: Commit final state**

```bash
git add -A
git commit -m "test: full test suite — ~103 tests pass with new QC module

QC subcommand: compute_qc_metrics (7 metrics), flag_qc_metrics (thresholds),
qc_decision (veto > rerun > caution > pass). Observation-only.
End-to-end verified with clean and problematic data.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```
