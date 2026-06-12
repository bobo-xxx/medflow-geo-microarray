library(testthat)
library(jsonlite)

source("../../node/scripts/exceptions.R")
source("../../node/scripts/report.R")

describe("write_run_result", {
  it("writes valid JSON with required fields", {
    d <- file.path(tempdir(), "run_result_test")
    dir.create(d)

    result <- list(
      status = "success_matrix",
      metadata = list(GPL570 = list(
        platform = "GPL570", organism = "Homo sapiens",
        n_samples = 12, n_probes = 54675,
        annotation_tier = 1, annotation_method = "fData:Gene Symbol",
        annotation_warning = NULL
      )),
      warnings = list("test warning"),
      probe_file = list("expr_probe.csv"),
      gene_file = list("expr_gene.csv"),
      meta_file = list("metadata.csv")
    )
    params <- list(subcommand = "fetch", gse_id = "GSE12345", outdir = d)

    write_run_result(d, result, params, 0, c("start","end"))

    path <- file.path(d, ".run_result.json")
    expect_true(file.exists(path))

    parsed <- fromJSON(path)
    expect_equal(parsed$node, "geo-microarray-processing")
    expect_equal(parsed$subcommand, "fetch")
    expect_equal(parsed$gse_id, "GSE12345")
    expect_equal(parsed$status, "success_matrix")
    expect_equal(parsed$exit_code, 0)
    expect_true("parameters" %in% names(parsed))
    expect_true("output" %in% names(parsed))
    expect_true("files" %in% names(parsed))
    expect_true("exceptions" %in% names(parsed))
    expect_true("warnings" %in% names(parsed))
    expect_true("started_at" %in% names(parsed))
    expect_true("finished_at" %in% names(parsed))

    unlink(d, recursive = TRUE)
  })

  it("includes annotation_tier in output section", {
    d <- file.path(tempdir(), "run_result_test2")
    dir.create(d)

    result <- list(
      status = "success_matrix",
      metadata = list(GPL570 = list(
        platform = "GPL570", organism = "Homo sapiens",
        n_samples = 12, n_probes = 54675,
        annotation_tier = 5, annotation_method = "probe_ids",
        annotation_warning = "Tier 1: no Symbol; Tier 2: no assignment; Tier 3: GPL empty; Tier 4: no BioC DB"
      )),
      warnings = list(), probe_file = list(), gene_file = NULL, meta_file = list()
    )
    params <- list(subcommand = "fetch", gse_id = "GSE12345", outdir = d)

    write_run_result(d, result, params, 0, c("start","end"))

    parsed <- fromJSON(file.path(d, ".run_result.json"))
    platform_data <- parsed$output[[1]]
    expect_equal(platform_data$annotation_tier, 5)
    expect_equal(platform_data$annotation_method, "probe_ids")
    expect_match(platform_data$annotation_warning, "Tier 1:")

    unlink(d, recursive = TRUE)
  })

  it("no gene file when annotation_tier is 5", {
    d <- file.path(tempdir(), "run_result_test3")
    dir.create(d)

    result <- list(
      status = "success_matrix",
      metadata = list(GPL570 = list(
        annotation_tier = 5, annotation_method = "probe_ids"
      )),
      warnings = list(), probe_file = list("probe.csv"),
      gene_file = NULL, meta_file = list()
    )
    params <- list(subcommand = "fetch", gse_id = "GSE12345", outdir = d)

    write_run_result(d, result, params, 0, c("start","end"))

    parsed <- fromJSON(file.path(d, ".run_result.json"))
    file_paths <- if (is.data.frame(parsed$files)) parsed$files$path else vapply(parsed$files, `[[`, "", "path")
    has_gene <- any(grepl("gene", file_paths))
    expect_false(has_gene)

    unlink(d, recursive = TRUE)
  })
})
