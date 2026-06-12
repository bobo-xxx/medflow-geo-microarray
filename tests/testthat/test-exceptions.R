library(testthat)

source("../../node/scripts/exceptions.R")

# ===========================================================================
# Unit tests
# ===========================================================================

describe("retry_with_backoff", {

  it("returns result on first success", {
    call_count <- 0
    fn <- function() { call_count <<- call_count + 1; "ok" }
    result <- retry_with_backoff(fn, max_attempts = 3, base_delay = 0)
    expect_equal(result, "ok")
    expect_equal(call_count, 1)
  })

  it("retries on failure and succeeds on 2nd attempt", {
    call_count <- 0
    fn <- function() {
      call_count <<- call_count + 1
      if (call_count < 2) stop("fail") else "ok"
    }
    result <- retry_with_backoff(fn, max_attempts = 3, base_delay = 0)
    expect_equal(result, "ok")
    expect_equal(call_count, 2)
  })

  it("returns NULL after all attempts exhausted", {
    call_count <- 0
    fn <- function() { call_count <<- call_count + 1; stop("always fail") }
    result <- retry_with_backoff(fn, max_attempts = 3, base_delay = 0)
    expect_null(result)
    expect_equal(call_count, 3)
  })
})

describe("safe_write_csv", {

  it("writes CSV and returns TRUE on success", {
    data <- data.frame(gene = c("A","B"), s1 = c(1.1, 2.2), s2 = c(3.3, 4.4))
    path <- file.path(tempdir(), "test_output.csv")
    result <- safe_write_csv(data, path)
    expect_true(result)
    expect_true(file.exists(path))
    unlink(path)
  })

  it("returns FALSE when write to unwritable path fails", {
    data <- data.frame(a = 1)
    path <- "/root/test_perm.csv"
    result <- safe_write_csv(data, path)
    expect_false(result)
  })

  it("does not leave temp file on success", {
    data <- data.frame(x = 1:3)
    path <- file.path(tempdir(), "test_notmp.csv")
    result <- safe_write_csv(data, path)
    expect_false(file.exists(paste0(path, ".tmp")))
    unlink(path)
  })
})

describe("check_environment", {

  it("returns OK when all requirements met", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) TRUE,
      .package = "base"
    )
    result <- check_environment()
    expect_equal(result$status, "ok")
  })

  it("returns error when GEOquery missing", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) pkg != "GEOquery",
      .package = "base"
    )
    result <- check_environment()
    expect_equal(result$status, "error")
    expect_match(paste(result$missing, collapse = ""), "GEOquery")
  })

  it("checks affy in required packages", {
    # affy should be in the required list
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) pkg != "affy",
      .package = "base"
    )
    result <- check_environment()
    expect_equal(result$status, "error")
    expect_match(paste(result$missing, collapse = ""), "affy")
  })
})

describe("detect_exception", {

  it("detects timeout pattern", {
    r <- detect_exception("connection timed out after 60 seconds")
    expect_equal(r$code, "A1_TIMEOUT")
    expect_equal(r$nature, "network")
  })

  it("detects 404 not found", {
    r <- detect_exception("HTTP 404 Not Found")
    expect_equal(r$code, "A2_NOT_FOUND")
  })

  it("detects permission denied", {
    r <- detect_exception("Permission denied")
    expect_equal(r$code, "W002_PERM_DENIED")
  })

  it("returns UNKNOWN for unrecognized errors", {
    r <- detect_exception("some random error message")
    expect_equal(r$code, "UNKNOWN")
  })
})

describe("write_checkpoint / read_checkpoint", {

  it("writes and reads checkpoints", {
    d <- file.path(tempdir(), "chkpt_test")
    dir.create(d)
    write_checkpoint(d, "download_raw", "complete")
    write_checkpoint(d, "process_gene", "in_progress")
    chk <- read_checkpoint(d)
    expect_equal(nrow(chk), 2)
    expect_equal(chk$step[1], "download_raw")
    expect_equal(chk$status[2], "in_progress")
    unlink(d, recursive = TRUE)
  })

  it("returns NULL for missing checkpoint file", {
    expect_null(read_checkpoint(tempdir()))
  })
})

describe("report_exception_ndjson", {

  it("emits valid NDJSON with required fields", {
    output <- capture.output(
      report_exception_ndjson("A1_TIMEOUT", "network", "retry",
                              "Connection timed out"),
      type = "output"
    )
    parsed <- jsonlite::fromJSON(output)
    expect_equal(parsed$level, "exception")
    expect_equal(parsed$code, "A1_TIMEOUT")
    expect_equal(parsed$nature, "network")
    expect_equal(parsed$action, "retry")
  })
})

describe("validate_cache", {

  it("returns valid when sentinel exists", {
    dir.create(d <- file.path(tempdir(), "cache_test"))
    writeLines("ok", file.path(d, ".fetch_complete"))
    result <- validate_cache(d)
    expect_equal(result$status, "valid")
    unlink(d, recursive = TRUE)
  })

  it("returns stale when sentinel missing", {
    dir.create(d <- file.path(tempdir(), "cache_test"))
    result <- validate_cache(d)
    expect_equal(result$status, "stale")
    unlink(d, recursive = TRUE)
  })
})

# ===========================================================================
# Flow / integration tests
# ===========================================================================

describe("Flow: safe_write_csv in pipeline context", {

  it("writes real expression matrix and verifies file is complete", {
    m <- matrix(runif(100 * 6, 2, 14), nrow = 100, ncol = 6)
    path <- file.path(tempdir(), "flow_test.csv")
    result <- safe_write_csv(m, path)
    expect_true(result)
    expect_true(file.exists(path))
    # Read back and verify dimensions
    back <- as.matrix(read.csv(path, row.names = 1))
    expect_equal(dim(back), dim(m))
    unlink(path)
  })
})

describe("Flow: check_environment → main pipeline guard", {

  it("prevents pipeline start when packages are missing", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) pkg != "GEOquery",
      .package = "base"
    )
    env <- check_environment()
    expect_equal(env$status, "error")
    # In main.R, this would trigger quit(status = 3)
    expect_true("GEOquery" %in% env$missing)
  })

  it("allows pipeline start when all packages present", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) TRUE,
      .package = "base"
    )
    env <- check_environment()
    expect_equal(env$status, "ok")
  })
})

describe("Flow: retry_with_backoff wraps fallible operations", {

  it("returns result when operation succeeds within retries", {
    attempts <- 0
    result <- retry_with_backoff(function() {
      attempts <<- attempts + 1
      if (attempts < 3) stop("transient error") else 42
    }, max_attempts = 3, base_delay = 0)
    expect_equal(result, 42)
    expect_equal(attempts, 3)
  })

  it("returns NULL when all retries exhausted (caller handles fallback)", {
    result <- retry_with_backoff(function() stop("permanent error"),
                                 max_attempts = 2, base_delay = 0)
    expect_null(result)
  })
})

describe("Flow: detect_exception guides error routing", {

  it("routes timeout to retry action", {
    r <- detect_exception("Operation timed out after 120s")
    expect_equal(r$action, "retry")
  })

  it("routes disk full to halt action", {
    r <- detect_exception("No space left on device")
    expect_equal(r$action, "halt")
  })
})
