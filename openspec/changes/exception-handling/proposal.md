## Why

The node currently has no structured exception handling. Network failures cause silent hangs. Partial downloads produce corrupted output. SIGTERM during a 30-minute download leaves stale files. Error messages are unstructured strings, not machine-readable codes. Agents consuming NDJSON output have no way to programmatically distinguish "skip this dataset" from "abort the pipeline."

## What Changes

- **New**: `node/scripts/exceptions.R` — centralized exception module with detection, retry, and NDJSON reporting
- **New**: Interactive prompt pattern — node writes prompt NDJSON, reads JSON response from stdin
- **New**: Cache validation — filelist.txt size comparison + `.fetch_complete` sentinel
- **New**: Atomic CSV writes via `safe_write_csv()` — temp file + verify + rename
- **New**: Environment startup check — package availability, network reachability, R version
- **Modified**: All processors — wrap GEO/network calls in `withTimeout()` + `retry_with_backoff()`
- **Modified**: `SKILL.md` exceptions — align with structured exception codes

## Capabilities

### New Capabilities
- `exception-handling`: Centralized exception detection, retry, backoff, atomic writes, cache validation, environment checks, lifecycle signal handling
- `interactive-prompts`: NDJSON prompt → stdin response pattern for user interaction on multi-failure

### Modified Capabilities
- None

## Impact

- `node/scripts/exceptions.R` — NEW
- `node/scripts/main.R` — add `check_environment()` call
- `node/scripts/fetch.R` — wrap network calls, use `safe_write_csv()`
- `node/scripts/raw.R` — wrap processing in timeout
- `node/scripts/report.R` — add `report_exception_ndjson()`
- `tests/testthat/test-exceptions.R` — NEW
