library(testthat)

source("tests/testthat/helpers.R")

describe("main.R CLI dispatch", {

  it("prints help when no arguments given", {
    result <- system2("Rscript", c("scripts/main.R"),
                      stdout = TRUE, stderr = TRUE)
    expect_true(any(grepl("Usage|usage|fetch|subcommand",
                          result, ignore.case = TRUE)),
      "Should print usage information when no args")
  })

  it("exits with error for unknown subcommand", {
    exit_code <- system2("Rscript", c("scripts/main.R", "unknown"),
                         stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 1)
  })

  it("exits with error when fetch missing --gse-id", {
    exit_code <- system2("Rscript", c("scripts/main.R", "fetch"),
                         stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 1)
  })

  it("accepts fetch subcommand with --gse-id", {
    # Stub implementation — just validates it doesn't crash
    exit_code <- system2("Rscript",
      c("scripts/main.R", "fetch", "--gse-id", "GSE100155"),
      stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 0)
  })

  it("accepts qc subcommand with --input", {
    exit_code <- system2("Rscript",
      c("scripts/main.R", "qc", "--input", "test.csv"),
      stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 0)
  })

  it("accepts clean subcommand with --input", {
    exit_code <- system2("Rscript",
      c("scripts/main.R", "clean", "--input", "test.csv"),
      stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 0)
  })
})

describe("NDJSON output", {

  it("produces valid JSON lines on stdout", {
    result <- system2("Rscript",
      c("scripts/main.R", "fetch", "--gse-id", "GSE100155"),
      stdout = TRUE, stderr = FALSE)

    # Each line should be parseable JSON
    for (line in result) {
      parsed <- tryCatch(jsonlite::fromJSON(line),
                         error = function(e) NULL)
      expect_false(is.null(parsed),
        sprintf("Line is not valid JSON: '%s'", substr(line, 1, 80)))
    }
  })

  it("has 'level' field in every NDJSON line", {
    result <- system2("Rscript",
      c("scripts/main.R", "fetch", "--gse-id", "GSE100155"),
      stdout = TRUE, stderr = FALSE)

    for (line in result) {
      parsed <- jsonlite::fromJSON(line)
      expect_true("level" %in% names(parsed),
        sprintf("NDJSON line missing 'level' field: '%s'", substr(line, 1, 80)))
    }
  })
})
