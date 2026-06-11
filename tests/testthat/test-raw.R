library(testthat)

source("../../node/scripts/raw.R")

describe("detect_raw_type", {

  it("detects Affymetrix from CEL files", {
    files <- c("GSM1.CEL.gz", "GSM2.CEL", "GSM3.CEL.gz")
    expect_equal(detect_raw_type(files), "affymetrix")
  })

  it("detects methylation from BPM + IDAT", {
    files <- c("sample1.idat", "sample2.idat", "manifest.bpm")
    expect_equal(detect_raw_type(files), "methylation")
  })

  it("detects Illumina from IDAT only (no BPM)", {
    files <- c("sample1.idat.gz", "sample2.idat")
    expect_equal(detect_raw_type(files), "illumina")
  })

  it("detects Agilent 2-color from GPR files", {
    files <- c("sample1.GPR", "sample2.GPR.gz")
    expect_equal(detect_raw_type(files), "agilent_2c")
  })

  it("detects Agilent FE from TXT with ProbeName/GeneName headers", {
    tmp <- tempfile(fileext = ".txt")
    writeLines(c("ProbeName\tGeneName\tgTotalGeneSignal\tsample1"), tmp)
    files <- c(tmp)  # use actual file path
    expect_equal(detect_raw_type(files), "agilent_1c")
    unlink(tmp)
  })

  it("returns unknown for unrecognized files", {
    files <- c("readme.txt", "notes.csv")
    expect_equal(detect_raw_type(files), "unknown")
  })

  it("checks TXT headers correctly for non-Agilent TXT", {
    tmp <- tempfile(fileext = ".txt")
    writeLines(c("gene_id\texpression\tpvalue"), tmp)
    files <- c(tmp)
    expect_equal(detect_raw_type(files), "unknown")
    unlink(tmp)
  })

  it("BPM check wins over IDAT (methylation before illumina)", {
    files <- c("sample1.idat", "manifest.BPM")
    expect_equal(detect_raw_type(files), "methylation")
  })
})

describe("detect_gpr_source", {

  it("returns genepix for ATF header", {
    tmp <- tempfile(fileext = ".GPR")
    writeLines('ATF\t1', tmp)
    expect_equal(detect_gpr_source(tmp), "genepix")
    unlink(tmp)
  })

  it("returns agilent for TYPE header", {
    tmp <- tempfile(fileext = ".GPR")
    writeLines('TYPE=GenePixResults', tmp)
    expect_equal(detect_gpr_source(tmp), "agilent")
    unlink(tmp)
  })

  it("returns genepix as default for unrecognized header", {
    tmp <- tempfile(fileext = ".GPR")
    writeLines('some\tother\theader', tmp)
    expect_equal(detect_gpr_source(tmp), "genepix")
    unlink(tmp)
  })
})

describe("process_raw_files", {

  it("returns error for empty file list", {
    result <- process_raw_files(character(0), tempdir(), "GSE12345")
    expect_equal(result$status, "error")
  })

  it("returns error for unknown type", {
    result <- process_raw_files(c("readme.txt"), tempdir(), "GSE12345")
    expect_equal(result$status, "error")
    expect_match(result$msg, "Unknown")
  })

  it("returns skip for methylation", {
    result <- process_raw_files(
      c("sample.idat", "manifest.BPM"), tempdir(), "GSE12345")
    expect_equal(result$status, "skipped_methylation")
  })
})

describe("Processor error handling", {

  it("process_illumina returns error when no IDAT in file list", {
    result <- process_illumina(c("notes.txt"), tempdir(), "GSE12345")
    expect_equal(result$status, "error")
    expect_match(result$msg, "No IDAT")
  })

  it("process_agilent_2c returns error when no GPR in file list", {
    result <- process_agilent_2c(c("notes.txt"), tempdir(), "GSE12345")
    expect_equal(result$status, "error")
    expect_match(result$msg, "No GPR")
  })

  it("all processors are defined as functions", {
    expect_true(is.function(process_affy))
    expect_true(is.function(process_illumina))
    expect_true(is.function(process_agilent_2c))
    expect_true(is.function(process_agilent_1c))
  })

  it("processors return list with required fields on error path", {
    r <- process_illumina(c("notes.txt"), tempdir(), "GSE12345")
    expect_equal(r$status, "error")
    expect_true("msg" %in% names(r))
  })
})

describe("process_affy — CEL via RMA", {

  it("applies ReadAffy → rma → exprs pipeline via affy (preferred)", {
    mock_matrix <- matrix(runif(100 * 4, 2, 14), nrow = 100, ncol = 4)

    local_mocked_bindings(
      requireNamespace = function(pkg, ...) pkg == "affy",
      .package = "base"
    )
    local_mocked_bindings(
      ReadAffy = function(filenames) list(),
      rma = function(raw, ...) list(),
      .package = "affy"
    )
    local_mocked_bindings(
      exprs = function(object) mock_matrix,
      .package = "Biobase"
    )

    result <- process_affy(c("test1.CEL", "test2.CEL"), tempdir(), "GSE12345")
    expect_equal(result$status, "success")
    expect_true(is.matrix(result$expr_matrix))
    expect_equal(result$platform, "Affymetrix")
    expect_match(result$pipeline, "affy::rma")
    expect_equal(dim(result$expr_matrix), c(100, 4))
  })

  it("falls back to oligo with limma processing when affy not available", {
    mock_matrix <- matrix(runif(100 * 4, 2, 14), nrow = 100, ncol = 4)

    local_mocked_bindings(
      requireNamespace = function(pkg, ...) pkg %in% c("oligo", "limma"),
      .package = "base"
    )
    local_mocked_bindings(
      read.celfiles = function(files) list(),
      .package = "oligo"
    )
    local_mocked_bindings(
      exprs = function(object) mock_matrix,
      .package = "Biobase"
    )
    local_mocked_bindings(
      normalizeBetweenArrays = function(x, method) x,
      .package = "limma"
    )

    result <- process_affy(c("test.CEL"), tempdir(), "GSE12345")
    expect_equal(result$status, "success")
    expect_match(result$pipeline, "oligo::read.celfiles")
  })

  it("returns error when neither affy nor oligo available", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) FALSE,
      .package = "base"
    )
    result <- process_affy(c("test.CEL"), tempdir(), "GSE12345")
    expect_equal(result$status, "error")
  })
})

describe("process_illumina — IDAT via neqc", {

  it("applies read.idat → neqc pipeline", {
    mock_matrix <- matrix(runif(100 * 4, 2, 14), nrow = 100, ncol = 4)
    # Create dummy idat files so gunzip/re-scan works
    tmpd <- file.path(tempdir(), "idat_test")
    dir.create(tmpd, showWarnings=FALSE)
    idat1 <- file.path(tmpd, "test1.idat"); file.create(idat1)
    idat2 <- file.path(tmpd, "test2.idat"); file.create(idat2)

    local_mocked_bindings(
      read.idat = function(files, bgxfile, ...) mock_matrix,
      neqc = function(x) list(E = log2(mock_matrix + 1e-6)),
      .package = "limma"
    )

    result <- process_illumina(c(idat1, idat2), tempdir(), "GSE12345")
    expect_equal(result$status, "success")
    expect_equal(result$platform, "Illumina")
    expect_match(result$pipeline, "neqc")
    expect_equal(dim(result$expr_matrix), c(100, 4))
    unlink(tmpd, recursive=TRUE)
  })

  it("returns error when requireNamespace fails for limma", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) pkg != "limma",
      .package = "base"
    )
    result <- process_illumina(c("test.idat"), tempdir(), "GSE12345")
    expect_equal(result$status, "error")
    expect_match(result$msg, "limma")
  })
})

describe("process_agilent_2c — GPR via normexp+loess+quantile", {

  it("applies full two-color pipeline", {
    mock_matrix <- matrix(runif(100 * 4, 2, 14), nrow = 100, ncol = 4)
    mock_rg <- list(A = mock_matrix)

    local_mocked_bindings(
      read.maimages = function(files, source, ...) mock_rg,
      backgroundCorrect = function(rg, method) rg,
      normalizeWithinArrays = function(rg, method) rg,
      normalizeBetweenArrays = function(rg, method) mock_matrix,
      .package = "limma"
    )

    result <- process_agilent_2c(c("test1.GPR", "test2.GPR"), tempdir(), "GSE12345")
    expect_equal(result$status, "success")
    expect_equal(result$platform, "Agilent_2C")
    expect_match(result$pipeline, "normexp")
    expect_match(result$pipeline, "loess")
    expect_match(result$pipeline, "quantile")
  })

  it("returns error when requireNamespace fails for limma", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) pkg != "limma",
      .package = "base"
    )
    result <- process_agilent_2c(c("test.GPR"), tempdir(), "GSE12345")
    expect_equal(result$status, "error")
  })
})

describe("process_agilent_1c — Agilent FE single-color", {

  it("applies source=agilent green.only pipeline", {
    mock_matrix <- matrix(runif(100 * 4, 2, 14), nrow = 100, ncol = 4)
    mock_rg <- list(E = mock_matrix)

    local_mocked_bindings(
      read.maimages = function(files, source, green.only, ...) mock_rg,
      normalizeBetweenArrays = function(rg, method) list(E = mock_matrix),
      .package = "limma"
    )

    result <- process_agilent_1c(c("test1.txt", "test2.txt"), tempdir(), "GSE12345")
    expect_equal(result$status, "success")
    expect_equal(result$platform, "Agilent_1C")
    expect_match(result$pipeline, "agilent")
    expect_match(result$pipeline, "quantile")
  })

  it("returns error when requireNamespace fails for limma", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) pkg != "limma",
      .package = "base"
    )
    result <- process_agilent_1c(c("test.txt"), tempdir(), "GSE12345")
    expect_equal(result$status, "error")
  })
})

