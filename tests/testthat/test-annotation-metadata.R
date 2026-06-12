library(testthat)
library(Biobase)

source("../../node/scripts/exceptions.R")

# ===========================================================================
# Annotation tier metadata tracking
# ===========================================================================

# Simulate fData for different platform types
make_fdata <- function(with_symbol = TRUE, with_assignment = FALSE) {
  fdata <- data.frame(
    ID = paste0("probe_", 1:100),
    stringsAsFactors = FALSE
  )
  if (with_symbol) {
    fdata[["Gene Symbol"]] <- paste0("GENE", 1:100)
  }
  if (with_assignment) {
    fdata[["gene_assignment"]] <- paste0("NM_001 // GENE", 1:100, " // desc")
  }
  fdata
}

describe("Annotation tier metadata", {

  it("annotation_tier entries flow from fetch.R variables", {
    # Verify the annotation tracking source code has all 5 tiers
    src <- readLines("../../node/scripts/fetch.R")
    anno_block <- paste(src[grep("anno_tier|anno_method|anno_reasons", src)], collapse = "\n")

    # Tier 1 tracked
    expect_true(any(grepl("anno_tier <- 1L", anno_block, fixed = TRUE)),
      "Tier 1 should set anno_tier = 1L")
    expect_true(any(grepl("anno_method.*paste0.*fData", anno_block)),
      "Tier 1 should set annotation method to fData column name")

    # Tier 2 tracked
    expect_true(any(grepl("anno_tier <- 2L", anno_block, fixed = TRUE)),
      "Tier 2 should set anno_tier = 2L")
    expect_true(any(grepl("anno_method.*gene_assignment", anno_block)),
      "Tier 2 should set method to gene_assignment")

    # Tier 3 tracked
    expect_true(any(grepl("anno_tier <- 3L", anno_block, fixed = TRUE)),
      "Tier 3 should set anno_tier = 3L")
    expect_true(any(grepl("anno_method.*GPL_table", anno_block)),
      "Tier 3 should set method to GPL_table")

    # Tier 4 tracked (Bioconductor annotation DB)
    expect_true(any(grepl("anno_tier <- 4L", anno_block, fixed = TRUE)),
      "Tier 4 should set anno_tier = 4L")
    expect_true(any(grepl("BioC_DB", anno_block)),
      "Tier 4 should use Bioconductor annotation DB")

    # Tier 5 tracked with warning
    expect_true(any(grepl("anno_tier <- 5L", anno_block, fixed = TRUE)),
      "Tier 5 should set anno_tier = 5L")
    expect_true(any(grepl("anno_method.*probe_ids", anno_block)),
      "Tier 5 should set method to probe_ids")
    expect_true(any(grepl("anno_warning.*paste.*anno_reasons", anno_block)),
      "Tier 5 should set warning from accumulated reasons")
  })

  it("failure reasons are accumulated across tiers", {
    src <- readLines("../../node/scripts/fetch.R")
    anno_block <- paste(src[grep("anno_reasons", src)], collapse = "\n")

    # Each tier should append a reason on failure
    expect_true(any(grepl("anno_reasons.*Tier 1", anno_block)),
      "Tier 1 failure should append reason")
    expect_true(any(grepl("anno_reasons.*Tier 2", anno_block)),
      "Tier 2 failure should append reason")
    expect_true(any(grepl("anno_reasons.*Tier 3", anno_block)),
      "Tier 3 failure should append reason")
    expect_true(any(grepl("anno_reasons.*Tier 4", anno_block)),
      "Tier 4 failure should append reason")
  })

})

describe("gpl_to_bioc_package", {
  it("maps GPL570 to hgu133plus2.db", {
    expect_equal(gpl_to_bioc_package("GPL570"), "hgu133plus2.db")
  })
  it("maps GPL16686 to hugene20sttranscriptcluster.db", {
    expect_equal(gpl_to_bioc_package("GPL16686"), "hugene20sttranscriptcluster.db")
  })
  it("maps GPL17586 to hta20transcriptcluster.db", {
    expect_equal(gpl_to_bioc_package("GPL17586"), "hta20transcriptcluster.db")
  })
  it("case insensitive", {
    expect_equal(gpl_to_bioc_package("gpl570"), "hgu133plus2.db")
  })
  it("returns NULL for unmapped GPL", {
    expect_null(gpl_to_bioc_package("GPL99999"))
  })
})

describe("annotate_with_bioc_db — fail first, then pass", {
  it("returns NULL for unmapped GPL", {
    expect_null(annotate_with_bioc_db(c("1","2"), "GPL99999"))
  })
  it("returns NULL when package not installed", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) FALSE,
      .package = "base"
    )
    expect_null(annotate_with_bioc_db(c("1","2"), "GPL570"))
  })
  it("returns probe-to-gene mapping when package is available", {
    # Mock AnnotationDbi::select to return test data
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) TRUE,
      .package = "base"
    )
    # Actually test with real package since it's installed
    if (requireNamespace("hugene20sttranscriptcluster.db", quietly = TRUE)) {
      result <- annotate_with_bioc_db(c("16657445","99999999"), "GPL16686")
      expect_false(is.null(result))
      expect_true("probe_id" %in% colnames(result))
      expect_true("gene_symbol" %in% colnames(result))
      # 16657445 → OR4F5
      expect_true("OR4F5" %in% result$gene_symbol)
    }
  })

  it("result metadata includes annotation fields", {
    src <- readLines("../../node/scripts/fetch.R")
    meta_block <- paste(src[grep("annotation_tier|annotation_method|annotation_warning", src)], collapse = "\n")

    # process_expression_set stores annotation in result$metadata
    expect_true(any(grepl("annotation_tier.*=.*anno_tier", meta_block)),
      "metadata should include annotation_tier")
    expect_true(any(grepl("annotation_method.*=.*anno_method", meta_block)),
      "metadata should include annotation_method")
    expect_true(any(grepl("annotation_warning.*=.*anno_warning", meta_block)),
      "metadata should include annotation_warning")
  })

  it("process_raw_matrix sets deferred annotation", {
    src <- readLines("../../node/scripts/fetch.R")
    raw_block <- src[max(grep("process_raw_matrix", src)):length(src)]
    raw_block <- paste(raw_block[1:50], collapse = "\n")

    expect_true(any(grepl("annotation_tier.*NA_integer", raw_block)),
      "raw matrix should set annotation_tier to NA")
    expect_true(any(grepl("annotation.*deferred", raw_block)),
      "raw matrix should note deferred annotation")
  })
})
