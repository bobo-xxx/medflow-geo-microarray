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
