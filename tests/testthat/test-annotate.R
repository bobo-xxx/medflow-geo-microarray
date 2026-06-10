library(testthat)

source("../../node/scripts/annotate.R")

describe("extract_gene_from_assignment", {

  it("extracts gene symbol from standard Affymetrix assignment format", {
    x <- "NM_001101 // ACTB // actin beta // ---"
    result <- extract_gene_from_assignment(x)
    expect_equal(result, "ACTB")
  })

  it("returns NA for malformed assignment (single field)", {
    x <- "NM_001101"
    result <- extract_gene_from_assignment(x)
    expect_true(is.na(result))
  })

  it("handles vector input", {
    x <- c(
      "NM_001101 // ACTB // actin beta // ---",
      "NM_002046 // GAPDH // description // ---",
      "NR_003286 // MALAT1 // long non-coding // ---"
    )
    result <- extract_gene_from_assignment(x)
    expect_equal(result, c("ACTB", "GAPDH", "MALAT1"))
  })

  it("handles empty string", {
    result <- extract_gene_from_assignment("")
    expect_true(is.na(result))
  })

  it("trims whitespace from gene symbols", {
    x <- "NM_001101 //  ACTB  // actin beta // ---"
    result <- extract_gene_from_assignment(x)
    expect_equal(result, "ACTB")
  })
})

describe("aggregate_probe_to_gene", {

  it("aggregates multiple probes to gene level by mean", {
    expr_matrix <- matrix(
      c(10, 20, 30, 5, 15, 25),
      nrow = 3, ncol = 2,
      dimnames = list(c("probe_1", "probe_2", "probe_3"), c("sample_A", "sample_B"))
    )

    gpl_table <- data.frame(
      probe_id = c("probe_1", "probe_2", "probe_3"),
      gene_symbol = c("GENE1", "GENE1", "GENE2"),
      stringsAsFactors = FALSE
    )

    result <- aggregate_probe_to_gene(expr_matrix, gpl_table)

    expect_equal(nrow(result), 2)  # GENE1, GENE2
    expect_equal(result["GENE1", "sample_A"], 15)
    expect_equal(result["GENE2", "sample_A"], 30)
  })

  it("handles ' /// ' separated gene symbols", {
    expr_matrix <- matrix(
      c(10, 20),
      nrow = 2, ncol = 1,
      dimnames = list(c("probe_1", "probe_2"), c("sample_X"))
    )

    gpl_table <- data.frame(
      probe_id = c("probe_1", "probe_2"),
      gene_symbol = c("GENE1 /// GENE2", "GENE3"),
      stringsAsFactors = FALSE
    )

    result <- aggregate_probe_to_gene(expr_matrix, gpl_table)
    expect_equal(result["GENE1", "sample_X"], 10)
    expect_equal(result["GENE2", "sample_X"], 10)
    expect_equal(result["GENE3", "sample_X"], 20)
  })

  it("returns NULL for NULL expr_matrix", {
    expect_null(aggregate_probe_to_gene(
      NULL, data.frame(probe_id = "p", gene_symbol = "G")))
  })

  it("returns NULL for NULL gpl_table", {
    m <- matrix(10, nrow = 1, ncol = 1, dimnames = list("p1", "s1"))
    expect_null(aggregate_probe_to_gene(m, NULL))
  })
})
