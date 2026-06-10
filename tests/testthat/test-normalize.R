library(testthat)

source("../../node/scripts/normalize.R")

describe("detect_expr_type", {

  it("returns 'raw' when 99th percentile > 100", {
    x <- c(rep(1, 98), rep(500, 2))
    expect_equal(detect_expr_type(x), "raw")
  })

  it("returns 'centered' when mean near 0 and range narrow", {
    x <- rnorm(100, mean = 0.1, sd = 1.5)
    expect_equal(detect_expr_type(x), "centered")
  })

  it("returns 'log' for typical microarray values (mid-range)", {
    x <- runif(100, 2, 14)
    expect_equal(detect_expr_type(x), "log")
  })

  it("handles NA values", {
    x <- c(runif(95, 2, 14), rep(NA, 5))
    result <- detect_expr_type(x)
    expect_true(result %in% c("raw", "centered", "log"))
  })
})

describe("normalize_expr_matrix", {

  it("applies log2(x + 1e-6) when type is raw", {
    m <- matrix(runif(20 * 5, 0, 500), nrow = 20, ncol = 5)
    result <- normalize_expr_matrix(m)
    expect_lt(max(result, na.rm = TRUE), max(m, na.rm = TRUE))
  })

  it("preserves centered data unchanged (already log2-scale)", {
    m <- matrix(rnorm(20 * 5, mean = 0.1, sd = 1.5), nrow = 20, ncol = 5)
    result <- normalize_expr_matrix(m)
    # Values should be unchanged (no shift applied)
    expect_equal(as.vector(result), as.vector(m))
  })

  it("preserves log-scale data unchanged", {
    m <- matrix(runif(20 * 5, 2, 14), nrow = 20, ncol = 5)
    result <- normalize_expr_matrix(m)
    expect_equal(result, m)
  })

  it("preserves matrix dimensions", {
    m <- matrix(runif(100 * 5, 0, 200), nrow = 100, ncol = 5)
    result <- normalize_expr_matrix(m)
    expect_equal(dim(result), dim(m))
  })
})
