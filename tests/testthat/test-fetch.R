library(testthat)

# Source dependencies in order (fetch.R depends on all of these)
source("../../scripts/normalize.R")
source("../../scripts/validate.R")
source("../../scripts/species.R")
source("../../scripts/annotate.R")
source("../../scripts/fetch.R")

describe("is_methylation", {

  it("returns TRUE when both BPM and IDAT files present", {
    files <- c("sample1.idat", "sample2.idat", "manifest.bpm")
    expect_true(is_methylation(files))
  })

  it("returns FALSE when only IDAT files (no BPM)", {
    files <- c("sample1.idat", "sample2.idat")
    expect_false(is_methylation(files))
  })

  it("returns FALSE for CEL files", {
    files <- c("sample1.CEL", "sample2.CEL.gz")
    expect_false(is_methylation(files))
  })

  it("case insensitive for .BPM", {
    files <- c("sample1.idat", "manifest.BPM")
    expect_true(is_methylation(files))
  })
})

describe("detect_platform_from_files", {

  it("detects Affymetrix from CEL files", {
    files <- c("GSM1.CEL.gz", "GSM2.CEL", "GSM3.CEL.gz")
    expect_equal(detect_platform_from_files(files), "Affymetrix")
  })

  it("detects Agilent from GPR files", {
    files <- c("sample1.GPR", "sample2.GPR.gz")
    expect_equal(detect_platform_from_files(files), "Agilent")
  })

  it("detects Illumina from idat-only files", {
    files <- c("sample1.idat", "sample2.idat")
    expect_equal(detect_platform_from_files(files), "Illumina")
  })

  it("returns Unknown for unrecognized files", {
    files <- c("sample1.txt", "sample2.csv")
    expect_equal(detect_platform_from_files(files), "Unknown")
  })
})

describe("do_fetch", {

  it("returns error status for invalid GSE ID", {
    result <- fetch_geo_data(list(
      gse_id = "GSE_INVALID_123456789",
      outdir = tempdir(),
      proxy = NULL,
      api_key = NULL
    ))
    expect_equal(result$status, "error")
  })

  it("returns error for non-GSE format ID", {
    result <- fetch_geo_data(list(
      gse_id = "INVALID",
      outdir = tempdir(),
      proxy = NULL,
      api_key = NULL
    ))
    expect_equal(result$status, "error")
  })
})
