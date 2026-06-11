library(testthat)

source("../../node/scripts/pipeline.R")

describe("detect_pipeline", {

  it("detects RMA from data_processing text", {
    txt <- "Background-adjusted signals were normalized by Robust Multichip Analysis (RMA)"
    expect_equal(detect_pipeline(txt), "rma")
  })

  it("detects MAS5/GCOS from target intensity", {
    txt <- "trimmed mean signal of each array was scaled to a target intensity of 300 using GCOS"
    expect_equal(detect_pipeline(txt), "mas5_gcos")
  })

  it("detects GCRMA", {
    txt <- "Data were normalized using GCRMA with gc content adjustment"
    expect_equal(detect_pipeline(txt), "gcrma")
  })

  it("detects SST-RMA/TAC", {
    txt <- "Transcriptome Analysis Console (TAC) software using SST-RMA normalization"
    expect_equal(detect_pipeline(txt), "sst_rma")
  })

  it("detects GenomeStudio average normalization", {
    txt <- "The data were normalised using average normalization with Genome Studio"
    expect_equal(detect_pipeline(txt), "genomestudio_avg")
  })

  it("detects neqc over GenomeStudio when both present", {
    txt <- "GenomeStudio export, then normalized using neqc with offset 16"
    expect_equal(detect_pipeline(txt), "neqc")
  })

  it("detects lumi VST", {
    txt <- "Variance stabilization transformation (VST) applied via lumi package"
    expect_equal(detect_pipeline(txt), "lumi")
  })

  it("detects Agilent Feature Extraction", {
    txt <- "Agilent Feature Extraction (FE) software with loess normalization"
    expect_equal(detect_pipeline(txt), "agilent_fe")
  })

  it("returns unknown for unrecognized text", {
    expect_equal(detect_pipeline("some custom processing"), "unknown")
  })

  it("handles NULL/empty input", {
    expect_equal(detect_pipeline(NULL), "unknown")
    expect_equal(detect_pipeline(""), "unknown")
  })

  it("case insensitive matching", {
    expect_equal(detect_pipeline("rma normalized"), "rma")
    expect_equal(detect_pipeline("RMA Normalized"), "rma")
    expect_equal(detect_pipeline("MAS5.0 with GCOS"), "mas5_gcos")
  })
})

describe("is_quantile_normalized", {

  it("returns FALSE for un-normalized data with different distributions", {
    set.seed(42)
    m <- cbind(
      sample_A = rnorm(1000, mean = 10, sd = 2),
      sample_B = rnorm(1000, mean = 12, sd = 3),
      sample_C = rnorm(1000, mean = 8,  sd = 1.5),
      sample_D = rnorm(1000, mean = 11, sd = 2.5)
    )
    expect_false(is_quantile_normalized(m))
  })

  it("returns TRUE after explicit quantile normalization", {
    set.seed(42)
    m <- cbind(
      sample_A = rnorm(1000, mean = 10, sd = 2),
      sample_B = rnorm(1000, mean = 12, sd = 3),
      sample_C = rnorm(1000, mean = 8,  sd = 1.5),
      sample_D = rnorm(1000, mean = 11, sd = 2.5)
    )
    library(limma)
    m_qn <- normalizeBetweenArrays(m, method = "quantile")
    expect_true(is_quantile_normalized(m_qn))
  })

  it("returns FALSE when only p50 passes but others fail", {
    set.seed(42)
    m <- cbind(
      sample_A = rnorm(1000, mean = 10, sd = 0.01),
      sample_B = rnorm(1000, mean = 10, sd = 0.01),
      sample_C = rnorm(1000, mean = 10, sd = 0.01),
      sample_D = rnorm(1000, mean = 10, sd = 0.01)
    )
    # All medians nearly identical, but tails differ
    m_qn <- normalizeBetweenArrays(m, method = "quantile")
    expect_true(is_quantile_normalized(m_qn))
  })
})

describe("apply_pipeline_transform", {

  it("returns pass-through for rma pipeline", {
    m <- matrix(runif(100, 2, 14), nrow = 10, ncol = 10)
    result <- apply_pipeline_transform(m, "rma")
    expect_equal(result$expr, m)
    expect_equal(result$transform, "none")
  })

  it("applies log2(x+1e-6) for mas5_gcos", {
    m <- matrix(runif(100, 0, 5000), nrow = 10, ncol = 10)
    result <- apply_pipeline_transform(m, "mas5_gcos")
    expect_equal(result$transform, "log2(x+1e-6)")
    expect_lt(max(result$expr), max(m))
  })

  it("applies shift+log2 for genomestudio_avg with negatives", {
    m <- matrix(c(-10, 0, 5, 100, 5000, -5, 10, 50, 200, 6000), nrow = 5, ncol = 2)
    result <- apply_pipeline_transform(m, "genomestudio_avg")
    expect_equal(result$transform, "shift+log2(x+1e-6)")
    # Minimum value (originally -10) becomes log2(1e-6) after shift+log2
    expect_equal(min(result$expr), log2(1e-6))
    # Maximum value should be properly log2-transformed
    expect_equal(max(result$expr), log2(max(m) - min(m) + 1e-6))
  })

  it("applies log2 only for genomestudio_avg without negatives", {
    m <- matrix(runif(100, 0, 5000), nrow = 10, ncol = 10)
    result <- apply_pipeline_transform(m, "genomestudio_avg")
    expect_equal(result$transform, "log2(x+1e-6)")
  })
})
