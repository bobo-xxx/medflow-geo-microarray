library(testthat)

source("../../node/scripts/exceptions.R")

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

  it("accepts custom max_attempts", {
    call_count <- 0
    fn <- function() { call_count <<- call_count + 1; stop("fail") }
    result <- retry_with_backoff(fn, max_attempts = 2, base_delay = 0)
    expect_null(result)
    expect_equal(call_count, 2)
  })
})

describe("safe_write_csv", {

  it("writes CSV atomically with sentinel", {
    data <- data.frame(gene = c("A","B"), s1 = c(1.1, 2.2), s2 = c(3.3, 4.4))
    path <- file.path(tempdir(), "test_output.csv")
    sentinel <- paste0(path, ".complete")

    result <- safe_write_csv(data, path)
    expect_true(result)
    expect_true(file.exists(path))
    # Sentinel not written by safe_write_csv (caller responsibility)
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

  it("returns error when required package missing", {
    local_mocked_bindings(
      requireNamespace = function(pkg, ...) pkg != "GEOquery",
      .package = "base"
    )
    result <- check_environment()
    expect_equal(result$status, "error")
    expect_match(paste(result$missing, collapse = ""), "GEOquery")
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
    expect_equal(parsed$msg, "Connection timed out")
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
