### Task 1: Create node/scripts/exceptions.R

- [x] `retry_with_backoff(fn, max_attempts, base_delay)` — exponential backoff wrapper
- [x] `with_geo_timeout(fn, timeout_sec)` — timeout wrapper (best-effort via R.utils)
- [x] `safe_write_csv(data, path)` — atomic write with temp file + verify + rename
- [x] `validate_cache(gse_dir)` — .fetch_complete sentinel check
- [x] `check_environment()` — package availability (GEOquery, Biobase, limma, affy)
- [x] `register_signal_handlers()` — documented no-op (R lacks portable signal handling)
- [x] `detect_exception(stderr_output)` — classify failure by pattern
- [x] `write_checkpoint(gse_dir, step, status)` / `read_checkpoint(gse_dir)`
- [x] `report_exception_ndjson(code, nature, action, msg, ...)` — structured NDJSON

### Task 2: Wire exceptions into main.R, fetch.R, report.R

- [x] `main.R`: call `check_environment()` at startup (register_signal_handlers is a documented no-op)
- [x] `fetCh.R`: wrap Tier 2 getGEO in `retry_with_backoff()`
- [x] `fetch.R`: replace `write.csv()` with `safe_write_csv()` + return value checks
- [x] `fetch.R`: write `.fetch_complete` sentinel and `write_checkpoint()` during processing
- [x] `report.R`: add `report_prompt()` (report_exception_ndjson lives in exceptions.R)
- [ ] `fetch.R`: call `validate_cache()` before Tier 1 (deferred — requires filelist.txt parsing)
- [ ] `fetch.R`: wrap Tier 3-4 downloads in `retry_with_backoff()` (deferred — tiers rarely hit)

### Task 3: Update SKILL.md exceptions

- [x] Replace pattern-based exceptions with 12 structured codes (A1-A3, B1-B3, C3, W1-W2, E801, E803, T5)
- [ ] Add `node/scripts/exceptions.R` to SKILL.md entry description (deferred — minor)

### Task 4: Create tests/testthat/test-exceptions.R

- [x] `retry_with_backoff`: success on 1st, 2nd attempt, all exhausted (3 tests)
- [x] `safe_write_csv`: success write, unwritable path, temp file cleanup (3 tests)
- [x] `validate_cache`: valid sentinel, missing sentinel (2 tests)
- [x] `check_environment`: all OK, missing GEOquery, missing affy (3 tests)
- [x] `detect_exception`: timeout, 404, perm denied, unknown (4 tests)
- [x] `write_checkpoint/read_checkpoint`: write+read, missing file (2 tests)
- [x] `report_exception_ndjson`: valid NDJSON fields (1 test)
- [x] Flow tests: safe_write_csv pipeline, check_environment guard, retry wrapping, detect_exception routing (6 tests)

### Task 5: Full test suite verification

- [x] 259 tests pass (0 FAIL)
- [x] Exceptions module: 24 unit + flow tests
- [ ] End-to-end smoke test on real data (network-dependent, blocked in this env)
