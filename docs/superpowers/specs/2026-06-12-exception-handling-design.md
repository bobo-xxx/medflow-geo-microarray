# Exception Handling — Design Doc

**Date:** 2026-06-12
**Status:** Design approved
**Pipeline:** superpowers brainstorming → OpenSpec change

## Overview

Comprehensive exception handling for the geo-microarray-processing node across all failure categories: network, data quality, resource, write-path, lifecycle, and environment.

## Architecture

```
node/scripts/exceptions.R
├── detect_exception()       — classify failure by pattern
├── retry_with_backoff()     — exponential backoff (network)
├── safe_write_csv()         — atomic write + verify
├── check_environment()      — startup package/network health
├── validate_cache()         — filelist.txt size check + .fetch_complete sentinel
└── report_exception_ndjson()— structured NDJSON with code, nature, action, prompt
```

## Categories

### A — Network Failures

| ID | Scenario | Detection | Action |
|---|---|---|---|
| A1 | Timeout | `withTimeout()` wrapper on every GEO call | Retry 3× with exponential backoff (5s, 20s, 80s) |
| A2 | Partial download | Compare actual file size to filelist.txt expected size | Re-download once; on 2nd failure, skip |
| A3 | 3 failures on same GSE | Counter per GSE | Prompt user: retry with proxy, try alternative GSE, skip, or abort |
| A4 | Multi-platform partial success | Per-platform status tracking | Retry failed platforms 3×; if still failing, prompt: continue_partial, retry_failed, abort |

### B — Data Quality Failures

| ID | Scenario | Detection | Action |
|---|---|---|---|
| B1 | File format / data type pre-check | `detect_raw_type()` + `detect_pipeline()` before processing | Route to correct processor or skip |
| B2 | Structured NDJSON exceptions | `report_exception_ndjson()` with code + nature + decision | Include in all result lines |
| B4 | Empty/broken matrix | `validate_expr_matrix()` (already implemented) | Reject, flag in warnings |

### C — Resource Failures

| ID | Scenario | Detection | Action |
|---|---|---|---|
| C3 | Thread creation failure | Catch `pthread_create` error from preprocessCore | Fall back to single-thread; warn via NDJSON |

### W — Write-Path Failures

| ID | Scenario | Detection | Action |
|---|---|---|---|
| W1 | Disk full | Temp file size < expected after write | halt (exit 1), clean up temp |
| W2 | Permission denied | `file.rename(tmp, final)` fails | halt (exit 1) |
| W3 | Path too long | `file.create(tmp)` fails | skip_with_warning |
| W4 | NFS stale handle | Temp file disappears after write | Retry 1×, then halt |

### L — Lifecycle Signals

| ID | Scenario | Detection | Action |
|---|---|---|---|
| L1 | SIGTERM | Signal handler | Stop download, delete temp, write checkpoint, exit |
| L2 | SIGINT | Signal handler | Stop download, keep partial for resume, prompt user |
| L3 | Normal completion | — | Delete checkpoint, write `.fetch_complete` sentinel |

### E — Environment Failures

| ID | Scenario | Detection | Action |
|---|---|---|---|
| E1 | Required package missing | `check_environment()` at startup | halt (exit 3), list missing packages |
| E2 | Package wrong version | `packageVersion()` | skip_with_warning |
| E3 | Network unreachable | `curl::nslookup("ftp.ncbi.nlm.nih.gov")` | halt |
| E4 | R version mismatch | `R.version$major` | skip_with_warning |

## Cache Validation (State Corruption)

```
GSE directory exists?
├── suppl/filelist.txt exists?
│   ├── Compare actual file sizes to expected
│   ├── Any mismatch → re-download that file
│   └── All match → intermediate files complete
├── .fetch_complete sentinel exists?
│   ├── Yes → use cached final outputs (Tier 1)
│   └── No → re-process from intermediate files
└── No directory → fresh fetch
```

Sentinel written atomically: `writeLines(..., tmp)` + `file.rename(tmp, ".fetch_complete")`.

## Checkpoint Format

`GSE12345/.fetch_checkpoint`:

```
download_raw|GPL96|complete|2026-06-12T10:30:00
process_gene|GPL96|complete|2026-06-12T10:30:15
download_raw|GPL97|in_progress|2026-06-12T10:31:00
```

On resume: skip `complete`, restart `in_progress` from beginning of that step.

## Graceful Degradation

```
All platforms succeed     → status: success
Some platforms succeed    → status: partial + per-platform breakdown
No platforms succeed      → status: error (fall to next tier)
```

NDJSON partial report includes platform-level status, summary counts, and interactive prompt options.

## Interactive Prompt Pattern

The node writes a prompt line to stdout, then reads a single JSON response line from stdin.
If no stdin is available (non-interactive/headless mode), the prompt times out after 300s
and the default action (first option) is taken.

```
Node writes:  {"level":"prompt",...} to stdout
Node reads:   {"choice":"retry_proxy","proxy":"http://localhost:2999"} from stdin
              (or times out → default action)
```

## Interactive Prompt Pattern

```json
{"level":"prompt","code":"A3_FAILED",
 "msg":"3 download attempts failed for GSE318047",
 "options":[
   {"choice":"retry_proxy","label":"Retry with proxy","prompt":"Enter proxy URL:"},
   {"choice":"alternative_gse","label":"Try a different GSE ID","prompt":"Enter GSE ID:"},
   {"choice":"skip","label":"Skip this dataset"},
   {"choice":"abort","label":"Abort processing"}
 ]}
```

## Atomic Write Pattern

All CSV writes go through `safe_write_csv(data, path)`:

```
1. Write to <path>.tmp
2. Verify tmp size > expected_min
3. file.rename(<path>.tmp, <path>)
4. If rename fails → W002_PERM_DENIED
5. If tmp too small → W001_DISK_FULL
```

## Implementation Order

1. `node/scripts/exceptions.R` — core module
2. Startup `check_environment()` — env validation
3. Network retry + interactive prompt (A1-A4)
4. Atomic write + sentinel (W1-W4, L1-L3)
5. Cache validation (state corruption)
6. Graceful degradation (partial success)
7. Structured NDJSON exceptions (B2)
