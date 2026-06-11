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

  it("detects NimbleGen from PAIR files", {
    files <- c("sample1.PAIR", "sample2.PAIR.gz")
    expect_equal(detect_raw_type(files), "nimblegen")
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

  it("process_nimblegen returns error when no PAIR in file list", {
    result <- process_nimblegen(c("notes.txt"), tempdir(), "GSE12345")
    expect_equal(result$status, "error")
    expect_match(result$msg, "No PAIR")
  })

  it("all processors are defined as functions", {
    expect_true(is.function(process_affy))
    expect_true(is.function(process_illumina))
    expect_true(is.function(process_agilent_2c))
    expect_true(is.function(process_agilent_1c))
    expect_true(is.function(process_nimblegen))
  })

  it("processors return list with required fields on error path", {
    r <- process_illumina(c("notes.txt"), tempdir(), "GSE12345")
    expect_equal(r$status, "error")
    expect_true("msg" %in% names(r))
  })
})
