library(testthat)

source("../../scripts/validate.R")

describe("validate_expr_matrix", {

  it("accepts a valid expression matrix", {
    m <- matrix(runif(100 * 5, 0, 20), nrow = 100, ncol = 5)
    colnames(m) <- paste0("sample_", 1:5)
    result <- validate_expr_matrix(m)
    expect_true(result$valid)
    expect_equal(result$n_rows, 100)
    expect_equal(result$n_cols, 5)
  })

  it("rejects NULL input", {
    result <- validate_expr_matrix(NULL)
    expect_false(result$valid)
  })

  it("rejects empty matrix (0 rows)", {
    m <- matrix(numeric(0), nrow = 0, ncol = 5)
    result <- validate_expr_matrix(m)
    expect_false(result$valid)
  })

  it("rejects empty matrix (0 columns)", {
    m <- matrix(runif(10), nrow = 10, ncol = 0)
    result <- validate_expr_matrix(m)
    expect_false(result$valid)
  })

  it("rejects matrix with extreme values (> 1e50)", {
    m <- matrix(runif(100 * 5, 0, 20), nrow = 100, ncol = 5)
    m[1, 1] <- 1e100
    result <- validate_expr_matrix(m)
    expect_false(result$valid)
    expect_match(result$reason, "[Ee]xtreme")
  })
})

describe("validate_gene_expression", {

  it("accepts valid gene expression data frame", {
    expr_gene <- data.frame(
      gene_symbol = c("GAPDH", "ACTB", "TP53"),
      sample_1 = c(10.5, 8.2, 12.1),
      sample_2 = c(11.0, 7.8, 11.9)
    )
    expect_true(validate_gene_expression(expr_gene))
  })

  it("rejects NULL", {
    expect_false(validate_gene_expression(NULL))
  })

  it("rejects empty data frame", {
    expect_false(validate_gene_expression(data.frame()))
  })
})

describe("validate_cel_integrity", {

  it("rejects non-existent file", {
    result <- validate_cel_integrity("/nonexistent/file.CEL")
    expect_false(result$valid)
  })
})
