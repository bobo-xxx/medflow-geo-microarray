library(testthat)

source("../../scripts/normalize.R")

describe("detect_expr_type", {

  it("returns 'raw' when 99th percentile > 100", {
    x <- c(rep(1, 98), rep(500, 2))  # 99pct ~ 500
    expect_equal(detect_expr_type(x), "raw")
  })

  it("returns 'centered' when mean near 0 and range narrow", {
    x <- rnorm(100, mean = 0.1, sd = 1.5)
    expect_equal(detect_expr_type(x), "centered")
  })

  it("returns 'log' for typical microarray values (mid-range)", {
    x <- runif(100, 2, 14)  # typical log2 expression range
    expect_equal(detect_expr_type(x), "log")
  })

  it("handles NA values", {
    x <- c(runif(95, 2, 14), rep(NA, 5))
    result <- detect_expr_type(x)
    expect_true(result %in% c("raw", "centered", "log"))
  })
})

describe("normalize_expr_matrix", {

  it("applies log2(x+1) when type is raw", {
    m <- matrix(runif(20 * 5, 0, 500), nrow = 20, ncol = 5)
    result <- normalize_expr_matrix(m)
    expect_true(all(result >= 0, na.rm = TRUE))
    expect_lt(max(result, na.rm = TRUE), max(m, na.rm = TRUE))
  })

  it("shifts negative values to non-negative", {
    m <- matrix(rnorm(20 * 5, -3, 2), nrow = 20, ncol = 5)
    result <- normalize_expr_matrix(m)
    expect_gte(min(result, na.rm = TRUE), 0)
  })

  it("preserves original values when already log-scale and non-negative", {
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
