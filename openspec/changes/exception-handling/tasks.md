### Task 1: Create node/scripts/exceptions.R

- [ ] `retry_with_backoff(fn, max_attempts, base_delay)` — exponential backoff wrapper
- [ ] `with_geo_timeout(fn, timeout_sec)` — timeout wrapper for GEO calls
- [ ] `safe_write_csv(data, path)` — atomic write with temp file + verify + rename
- [ ] `validate_cache(gse_dir)` — filelist.txt size check + .fetch_complete sentinel
- [ ] `check_environment()` — package availability, network reachability, R version
- [ ] `register_signal_handlers()` — SIGTERM/SIGINT handlers
- [ ] `detect_exception(stderr_output)` — classify failure by pattern
- [ ] `write_checkpoint(gse_dir, step, status)` / `read_checkpoint(gse_dir)`
- [ ] `report_exception_ndjson(code, nature, action, msg, ...)` — structured NDJSON

### Task 2: Wire exceptions into main.R, fetch.R, report.R

- [ ] `main.R`: call `check_environment()` and `register_signal_handlers()` at startup
- [ ] `fetch.R`: wrap all GEO calls in `with_geo_timeout()` + `retry_with_backoff()`
- [ ] `fetch.R`: replace `write.csv()` with `safe_write_csv()`
- [ ] `fetch.R`: call `validate_cache()` before Tier 1
- [ ] `report.R`: add `report_exception_ndjson()` and `report_prompt()` functions

### Task 3: Update SKILL.md exceptions

- [ ] Replace pattern-based exceptions with structured codes (A1-A4, B1-B4, C3, W1-W4, L1-L3, E1-E4)
- [ ] Add `node/scripts/exceptions.R` to entry description

### Task 4: Create tests/testthat/test-exceptions.R

- [ ] `retry_with_backoff`: success on 2nd attempt, all exhausted
- [ ] `safe_write_csv`: normal write, disk full simulation, permission denied
- [ ] `validate_cache`: valid sentinel, missing sentinel, size mismatch
- [ ] `check_environment`: all OK, missing package, no network
- [ ] `report_exception_ndjson`: valid NDJSON with all required fields
- [ ] Prompt contract: valid JSON format, option structure, timeout behavior

### Task 5: Full test suite verification

- [ ] All existing tests pass
- [ ] New exception tests pass
- [ ] End-to-end smoke test on real data
