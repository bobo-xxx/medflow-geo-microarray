library(testthat)

source("../../node/scripts/exceptions.R")

describe("retry_with_backoff", {
  it("returns result on first success", {
    call_count <- 0
    fn <- function() { call_count <<- call_count + 1; "ok" }
    result <- retry_with_backoff(fn, max_attempts = 3, base_delay = 0)
    expect_equal(result, "ok"); expect_equal(call_count, 1)
  })
  it("succeeds on 2nd attempt", {
    call_count <- 0
    fn <- function() { call_count <<- call_count + 1; if (call_count < 2) stop("fail") else "ok" }
    result <- retry_with_backoff(fn, max_attempts = 3, base_delay = 0)
    expect_equal(result, "ok"); expect_equal(call_count, 2)
  })
  it("returns NULL after all exhausted", {
    call_count <- 0
    fn <- function() { call_count <<- call_count + 1; stop("always fail") }
    result <- retry_with_backoff(fn, max_attempts = 3, base_delay = 0)
    expect_null(result); expect_equal(call_count, 3)
  })
})

describe("safe_write_csv", {
  it("writes CSV and returns TRUE on success", {
    data <- data.frame(gene = c("A","B"), s1 = c(1.1, 2.2), s2 = c(3.3, 4.4))
    path <- file.path(tempdir(), "test_output.csv")
    result <- safe_write_csv(data, path)
    expect_true(result); expect_true(file.exists(path)); unlink(path)
  })
})

describe("check_environment", {
  it("returns OK when all requirements met", {
    local_mocked_bindings(requireNamespace = function(pkg, ...) TRUE, .package = "base")
    expect_equal(check_environment()$status, "ok")
  })
  it("returns error when GEOquery missing", {
    local_mocked_bindings(requireNamespace = function(pkg, ...) pkg != "GEOquery", .package = "base")
    r <- check_environment()
    expect_equal(r$status, "error"); expect_match(paste(r$missing, collapse=""), "GEOquery")
  })
  it("checks affy in required packages", {
    local_mocked_bindings(requireNamespace = function(pkg, ...) pkg != "affy", .package = "base")
    r <- check_environment()
    expect_equal(r$status, "error"); expect_match(paste(r$missing, collapse=""), "affy")
  })
})

describe("detect_exception", {
  it("detects timeout", { r <- detect_exception("connection timed out"); expect_equal(r$code, "A1_TIMEOUT") })
  it("detects 404", { r <- detect_exception("HTTP 404"); expect_equal(r$code, "A2_NOT_FOUND") })
  it("detects all attempts exhausted", { r <- detect_exception("All attempts exhausted"); expect_equal(r$code, "A3_FAILED") })
  it("detects methylation BPM", { r <- detect_exception("Methylation BPM present"); expect_equal(r$code, "B2_METHYLATION") })
  it("detects file not found", { r <- detect_exception("File not found: /path"); expect_equal(r$code, "B3_EMPTY") })
  it("detects all tiers failed", { r <- detect_exception("All data retrieval methods failed"); expect_equal(r$code, "T5_ALL_FAILED") })
  it("detects metadata only fallback", { r <- detect_exception("returning metadata only"); expect_equal(r$code, "T5_ALL_FAILED") })
  it("detects disk full", { r <- detect_exception("No space left on device"); expect_equal(r$code, "W001_DISK_FULL") })
  it("detects perm denied", { r <- detect_exception("Permission denied"); expect_equal(r$code, "W002_PERM_DENIED") })
  it("detects thread error", { r <- detect_exception("pthread_create failed"); expect_equal(r$code, "C3_THREAD") })
  it("detects unreachable", { r <- detect_exception("ftp.ncbi.nlm.nih.gov unreachable"); expect_equal(r$code, "E803_ENV_NET") })
  it("returns UNKNOWN for unrecognized", { r <- detect_exception("random error"); expect_equal(r$code, "UNKNOWN") })
})

describe("write_checkpoint / read_checkpoint", {
  it("writes and reads checkpoints", {
    d <- file.path(tempdir(), "chkpt_test"); dir.create(d)
    write_checkpoint(d, "download_raw", "complete")
    write_checkpoint(d, "process_gene", "in_progress")
    chk <- read_checkpoint(d)
    expect_equal(nrow(chk), 2); expect_equal(chk$step[1], "download_raw")
    unlink(d, recursive = TRUE)
  })
  it("returns NULL for missing checkpoint file", {
    expect_null(read_checkpoint(tempdir()))
  })
})

describe("report_exception_ndjson", {
  it("emits valid NDJSON with required fields", {
    output <- capture.output(
      report_exception_ndjson("A1_TIMEOUT", "network", "retry", "Connection timed out"), type = "output")
    parsed <- jsonlite::fromJSON(output)
    expect_equal(parsed$level, "exception"); expect_equal(parsed$code, "A1_TIMEOUT")
  })
})

describe("validate_cache", {
  it("returns valid when sentinel exists", {
    d <- file.path(tempdir(), "cache_test"); dir.create(d)
    writeLines("ok", file.path(d, ".fetch_complete"))
    expect_equal(validate_cache(d)$status, "valid"); unlink(d, recursive = TRUE)
  })
  it("returns stale when sentinel missing", {
    d <- file.path(tempdir(), "cache_test"); dir.create(d)
    expect_equal(validate_cache(d)$status, "stale"); unlink(d, recursive = TRUE)
  })
})

# Flow tests
describe("Flow: report_and_classify wiring", {
  it("classifies and emits NDJSON for timeout", {
    output <- capture.output(report_and_classify("Connection timed out after 60s"), type = "output")
    parsed <- jsonlite::fromJSON(output)
    expect_equal(parsed$code, "A1_TIMEOUT"); expect_equal(parsed$action, "retry")
  })
  it("handles unknown errors", {
    output <- capture.output(report_and_classify("some random failure"), type = "output")
    parsed <- jsonlite::fromJSON(output)
    expect_equal(parsed$code, "UNKNOWN"); expect_equal(parsed$action, "escalate")
  })
})

describe("Flow: safe_write_csv pipeline", {
  it("writes expression matrix and verifies", {
    m <- matrix(runif(100 * 6, 2, 14), nrow = 100, ncol = 6)
    path <- file.path(tempdir(), "flow_test.csv")
    result <- safe_write_csv(m, path)
    expect_true(result); expect_true(file.exists(path))
    back <- as.matrix(read.csv(path, row.names = 1))
    expect_equal(dim(back), dim(m)); unlink(path)
  })
})

describe("detect_exception: remaining patterns", {
  it("detects unknown format", {
    r <- detect_exception("unknown format: cannot read file")
    expect_equal(r$code, "B1_FORMAT")
    expect_equal(r$nature, "data_corrupt")
  })
  it("detects corrupt data", {
    r <- detect_exception("corrupted CEL file detected")
    expect_equal(r$code, "B1_FORMAT")
  })
  it("detects env net unreachable", {
    r <- detect_exception("ftp.ncbi.nlm.nih.gov unreachable: DNS failure")
    expect_equal(r$code, "E803_ENV_NET")
    expect_equal(r$action, "halt")
  })
})

describe("Flow: E801_ENV_PKG fires via report_exception_ndjson", {
  it("emits valid NDJSON for env pkg error", {
    output <- capture.output(
      report_exception_ndjson("E801_ENV_PKG", "env_bug", "halt",
        "Missing required packages: GEOquery"),
      type = "output"
    )
    parsed <- jsonlite::fromJSON(output)
    expect_equal(parsed$code, "E801_ENV_PKG")
    expect_equal(parsed$action, "halt")
  })
})

describe("Flow: retry wraps fallible ops", {
  it("returns result on eventual success", {
    attempts <- 0
    result <- retry_with_backoff(function() {
      attempts <<- attempts + 1; if (attempts < 3) stop("transient") else 42
    }, max_attempts = 3, base_delay = 0)
    expect_equal(result, 42); expect_equal(attempts, 3)
  })
})
