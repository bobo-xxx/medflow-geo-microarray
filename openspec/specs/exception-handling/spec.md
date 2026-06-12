# exception-handling Specification

## Purpose
TBD - created by archiving change exception-handling. Update Purpose after archive.
## Requirements
### Requirement: Network retry with exponential backoff

The system SHALL retry failed network operations with exponential backoff.

#### Scenario: Successful retry

- **WHEN** a GEO network call fails with timeout or connection error
- **THEN** the system SHALL retry up to 3 times with delays 5s, 20s, 80s
- **AND** report each attempt via NDJSON info line

#### Scenario: All retries exhausted

- **WHEN** all 3 retry attempts fail
- **THEN** the system SHALL emit a prompt-level NDJSON line with code A3_FAILED
- **AND** await user response or timeout to default action

### Requirement: Atomic CSV writes

The system SHALL write all CSV output files atomically.

#### Scenario: Successful atomic write

- **WHEN** writing a CSV output file
- **THEN** data SHALL be written to a `.tmp` file first
- **AND** verified for minimum expected size
- **AND** atomically renamed to the final path

#### Scenario: Disk full detected

- **WHEN** temp file size is less than 50% of expected minimum
- **THEN** the system SHALL report exception W001_DISK_FULL and halt

#### Scenario: Permission denied

- **WHEN** `file.rename(tmp, final)` fails
- **THEN** the system SHALL report exception W002_PERM_DENIED and halt

### Requirement: Cache validation

The system SHALL validate cached intermediate files before reusing them.

#### Scenario: Valid cache detected

- **WHEN** `suppl/filelist.txt` exists and all file sizes match
- **AND** `.fetch_complete` sentinel exists
- **THEN** the system SHALL use cached files without re-downloading

#### Scenario: Stale cache detected

- **WHEN** `.fetch_complete` sentinel is missing
- **OR** file sizes don't match filelist.txt
- **THEN** the system SHALL re-download mismatched files

### Requirement: Environment startup check

The system SHALL validate the runtime environment before processing.

#### Scenario: Required package missing

- **WHEN** any of GEOquery, Biobase, limma, or affy are not available
- **THEN** the system SHALL report exception E801_ENV_PKG and halt with exit code 3

#### Scenario: Network unreachable

- **WHEN** `ftp.ncbi.nlm.nih.gov` is not reachable via DNS
- **THEN** the system SHALL report exception E803_ENV_NET and halt

#### Scenario: Environment OK

- **WHEN** all required packages and network are available
- **THEN** processing SHALL proceed normally

### Requirement: Lifecycle signal handling

The system SHALL handle SIGTERM and SIGINT gracefully.

#### Scenario: SIGTERM received

- **WHEN** SIGTERM is received during download
- **THEN** the system SHALL stop the current download
- **AND** delete incomplete temp files
- **AND** write a checkpoint file with current progress
- **AND** exit

#### Scenario: SIGINT received

- **WHEN** SIGINT is received during download
- **THEN** the system SHALL stop the current download
- **AND** keep partial download for resume
- **AND** prompt the user whether to resume on next run
- **AND** exit

### Requirement: Graceful degradation for multi-platform

The system SHALL report partial success when some platforms fail.

#### Scenario: Partial multi-platform success

- **WHEN** some platforms succeed and some fail after retries
- **THEN** the system SHALL return status `partial`
- **AND** include per-platform status breakdown in NDJSON result
- **AND** prompt the user to continue, retry, or abort

### Requirement: Structured NDJSON exception reporting

The system SHALL report all exceptions as structured NDJSON with machine-readable codes. Exception codes SHALL be emitted at the error site, not deferred to the caller.

#### Scenario: Exception reported at error site

- **WHEN** any exception condition is detected
- **THEN** the NDJSON line SHALL be emitted immediately at the detection point
- **AND** include `code`, `nature`, `action`, and `msg` fields
- **AND** the `nature` field SHALL be one of: data_insufficient, data_corrupt, data_mismatch, env_bug, network, resource

### Requirement: Error sites SHALL emit structured NDJSON exceptions

Every error-producing operation in the node SHALL call `detect_exception()` to classify the error, then `report_exception_ndjson()` to emit a structured NDJSON line with machine-readable code, nature, and action fields.

#### Scenario: Network error at Tier 2

- **WHEN** `getGEO()` fails after all retries
- **THEN** the system SHALL call `detect_exception(stderr)` to classify
- **AND** emit `{"level":"exception","code":"A3_FAILED",...}` via `report_exception_ndjson()`

#### Scenario: Data corruption detected

- **WHEN** `validate_expr_matrix()` returns `valid=FALSE`
- **THEN** the system SHALL classify via `detect_exception(reason)` 
- **AND** emit structured NDJSON with the appropriate B-code

#### Scenario: Write failure detected

- **WHEN** `safe_write_csv()` returns FALSE
- **THEN** the system SHALL classify via `detect_exception("write failed")`
- **AND** emit structured NDJSON with code W001 or W002

#### Scenario: Processor fails in raw.R

- **WHEN** any processor returns `status="error"`
- **THEN** the system SHALL emit structured NDJSON with the processor's error message classified

#### Scenario: Methylation detected

- **WHEN** BPM + IDAT files are found
- **THEN** the system SHALL emit `{"level":"exception","code":"B2_METHYLATION","action":"skip_with_warning",...}`

#### Scenario: All tiers exhausted

- **WHEN** Tiers 2-4 all fail and Tier 5 is reached
- **THEN** the system SHALL emit `{"level":"exception","code":"T5_ALL_FAILED","action":"halt",...}`

#### Scenario: File not found in qc/clean

- **WHEN** input file does not exist for qc or clean subcommand
- **THEN** the system SHALL emit structured NDJSON with code B3_EMPTY

