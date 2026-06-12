## Architecture

```
node/scripts/exceptions.R
├── retry_with_backoff(fn, max_attempts=3, base_delay=5)
├── with_geo_timeout(fn, timeout_sec=300)
├── safe_write_csv(data, path)
├── validate_cache(gse_dir)
├── check_environment()
├── register_signal_handlers()
├── detect_exception(stderr_output)
├── write_checkpoint(gse_dir, step, status)
├── read_checkpoint(gse_dir)
└── report_exception_ndjson(code, nature, action, ...)
```

## Exception Codes

All codes follow `<category><number>_<name>` pattern. Categories: A=Network, B=Data, C=Resource, W=Write, L=Lifecycle, E=Environment.

## Interactive Prompt Contract

```
Node → stdout: {"level":"prompt","code":"...","msg":"...","options":[...]}
Node ← stdin:  {"choice":"<option>","<param>":"<value>"}
Timeout: 300s → default action (first option)
Headless: no stdin → default action immediately
```

## Files

| File | Change |
|---|---|
| `node/scripts/exceptions.R` | NEW — all exception handling logic |
| `node/scripts/main.R` | MODIFIED — call check_environment(), register signal handlers |
| `node/scripts/fetch.R` | MODIFIED — wrap GEO calls, safe_write_csv, validate_cache |
| `node/scripts/report.R` | MODIFIED — add report_exception_ndjson() |
| `node/SKILL.md` | MODIFIED — update exceptions with structured codes |
| `tests/testthat/test-exceptions.R` | NEW |
