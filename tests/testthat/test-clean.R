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
    expect_equal(as.vector(cleaned), as.vector(m))
    expect_equal(dim(cleaned), dim(m))

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
