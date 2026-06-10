# add-clean-subcommand Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `clean` subcommand that loads an expression matrix CSV, applies normalization transforms (log2 for raw, shift for centered), and writes a cleaned CSV with NDJSON reporting.

**Architecture:** `node/scripts/clean.R` reuses `normalize_expr_matrix()` from normalize.R. It reads input CSV → detects scale → normalizes → writes output CSV → reports via NDJSON. The transform is deterministic: raw → log2(x+1), centered → shift to non-negative, log → pass-through. `main.R`'s `do_clean()` calls it.

**Tech Stack:** R >= 4.3, stats

---

## File Structure

| File | Purpose |
|---|---|
| `node/scripts/clean.R` | `run_clean(input, output)` — load, normalize, write |
| `tests/testthat/test-clean.R` | Unit tests with synthetic CSV fixtures |
| `node/scripts/main.R` | Replace `do_clean` stub |

---

### Task 1: Create node/scripts/clean.R + tests/testthat/test-clean.R

**Files:**
- Create: `node/scripts/clean.R`
- Create: `tests/testthat/test-clean.R`

- [ ] **Step 1: Write test-clean.R**

```r
library(testthat)

source("../../node/scripts/normalize.R")
source("../../node/scripts/clean.R")

describe("run_clean", {

  it("normalizes raw-scale data to log2", {
    set.seed(42)
    m <- matrix(runif(500 * 6, 0, 500), nrow = 500, ncol = 6)
    colnames(m) <- paste0("sample_", 1:6)
    write.csv(m, "test_raw.csv", row.names = TRUE)

    result <- run_clean("test_raw.csv", "test_clean.csv")
    expect_equal(result$status, "success")
    expect_equal(result$input_scale, "raw")
    expect_equal(result$output_scale, "log")

    cleaned <- as.matrix(read.csv("test_clean.csv", row.names = 1))
    expect_equal(dim(cleaned), dim(m))
    # log2(x+1) should reduce max
    expect_lt(max(cleaned, na.rm = TRUE), max(m, na.rm = TRUE))

    unlink("test_raw.csv")
    unlink("test_clean.csv")
  })

  it("preserves log-scale data unchanged", {
    set.seed(42)
    m <- matrix(runif(500 * 6, 2, 14), nrow = 500, ncol = 6)
    colnames(m) <- paste0("sample_", 1:6)
    write.csv(m, "test_log.csv", row.names = TRUE)

    result <- run_clean("test_log.csv", "test_log_clean.csv")
    expect_equal(result$status, "success")
    expect_equal(result$input_scale, "log")
    expect_equal(result$output_scale, "log")
    expect_equal(result$applied_transform, "none")

    cleaned <- as.matrix(read.csv("test_log_clean.csv", row.names = 1))
    expect_equal(cleaned, m)

    unlink("test_log.csv")
    unlink("test_log_clean.csv")
  })

  it("shifts centered data to non-negative", {
    set.seed(42)
    m <- matrix(rnorm(500 * 6, mean = 0.1, sd = 1.5), nrow = 500, ncol = 6)
    colnames(m) <- paste0("sample_", 1:6)
    write.csv(m, "test_centered.csv", row.names = TRUE)

    result <- run_clean("test_centered.csv", "test_shifted.csv")
    expect_equal(result$status, "success")
    expect_equal(result$input_scale, "centered")
    expect_equal(result$applied_transform, "shift")

    cleaned <- as.matrix(read.csv("test_shifted.csv", row.names = 1))
    expect_gte(min(cleaned, na.rm = TRUE), 0)

    unlink("test_centered.csv")
    unlink("test_shifted.csv")
  })

  it("returns error for non-existent input", {
    result <- run_clean("/nonexistent/input.csv", "/tmp/out.csv")
    expect_equal(result$status, "error")
  })

  it("preserves row and column names", {
    set.seed(42)
    m <- matrix(runif(100 * 4, 0, 200), nrow = 100, ncol = 4)
    rownames(m) <- paste0("gene_", 1:100)
    colnames(m) <- paste0("sample_", 1:4)
    write.csv(m, "test_named.csv", row.names = TRUE)

    result <- run_clean("test_named.csv", "test_named_clean.csv")
    expect_equal(result$status, "success")

    cleaned <- as.matrix(read.csv("test_named_clean.csv", row.names = 1))
    expect_equal(rownames(cleaned), rownames(m))
    expect_equal(colnames(cleaned), colnames(m))

    unlink("test_named.csv")
    unlink("test_named_clean.csv")
  })

  it("handles data with NA values", {
    m <- matrix(runif(100 * 6, 0, 300), nrow = 100, ncol = 6)
    m[1:10, 1] <- NA
    colnames(m) <- paste0("s", 1:6)
    write.csv(m, "test_na.csv", row.names = TRUE)

    result <- run_clean("test_na.csv", "test_na_clean.csv")
    expect_equal(result$status, "success")

    unlink("test_na.csv")
    unlink("test_na_clean.csv")
  })
})
```

- [ ] **Step 2: Run test, verify FAILS**

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate geo-microarray-processing
Rscript -e 'testthat::test_file("tests/testthat/test-clean.R")'
```

Expected: FAIL

- [ ] **Step 3: Write node/scripts/clean.R**

```r
# clean.R — Expression matrix normalization for geo-microarray-processing
#
# Applies normalization transforms to expression matrices:
#   raw      -> log2(x + 1)
#   centered -> shift to non-negative
#   log      -> pass-through
#
# Reuses normalize_expr_matrix() from normalize.R.

#' Clean/normalize an expression matrix
#'
#' Reads input CSV, detects expression scale, applies normalization,
#' writes cleaned CSV. Reports transform details.
#'
#' @param input  Path to input expression matrix CSV
#' @param output Path for cleaned output CSV
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
    "log2(x+1)"
  } else if (input_scale == "centered" && min(expr, na.rm = TRUE) < 0) {
    "shift"
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
```

- [ ] **Step 4: Run test, verify PASS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-clean.R")'
```

Expected: ~7/7 PASS

- [ ] **Step 5: Commit**

```bash
git add node/scripts/clean.R tests/testthat/test-clean.R
git commit -m "feat: add clean.R with run_clean — normalize expression matrices

Applies normalize_expr_matrix() transforms: raw->log2(x+1),
centered->shift, log->pass-through. Preserves row/col names,
handles NA values. Reports input/output scale and applied transform.
~7 tests covering raw, log, centered, NA, missing file, and name preservation.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Wire clean.R into main.R

**Files:**
- Modify: `node/scripts/main.R`

- [ ] **Step 1: Add source and replace do_clean stub**

Add `source(file.path(script_dir, "clean.R"))` after the qc.R source line.

Replace the do_clean stub:

```r
#' Clean subcommand
do_clean <- function(opts) {
  report_info(sprintf("Cleaning %s...", opts$input))
  output <- if (!is.null(opts$output)) opts$output else sub("\\.csv$", "_clean.csv", opts$input)
  result <- run_clean(opts$input, output)

  if (result$status == "error") {
    report_error(result$msg)
    return(invisible(NULL))
  }

  report_info(sprintf("Input scale: %s, applied: %s, output scale: %s",
    result$input_scale, result$applied_transform, result$output_scale))

  report_result(result$status,
    files = list(list(path = result$output_path, rows = result$n_rows, cols = result$n_cols)),
    metadata = list(
      input_scale       = result$input_scale,
      output_scale      = result$output_scale,
      applied_transform = result$applied_transform
    ))
}
```

Also add `--output` to the arg parser's `opts` list and parsing loop.

- [ ] **Step 2: Run full test suite**

```bash
Rscript -e 'testthat::test_dir("tests/testthat/")'
```

Expected: 0 FAIL, ~121 tests

- [ ] **Step 3: Commit**

```bash
git add node/scripts/main.R
git commit -m "feat: wire clean subcommand to clean.R with NDJSON reporting

Replace do_clean stub with run_clean() call. Cleans output includes
transform details and output file path. Added --output arg support.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```
