## ADDED Requirements

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

The system SHALL report all exceptions as structured NDJSON with machine-readable codes.

#### Scenario: Exception reported

- **WHEN** any exception condition is detected
- **THEN** the NDJSON line SHALL include `code`, `nature`, `action`, and `msg` fields
- **AND** the `nature` field SHALL be one of: data_insufficient, data_corrupt, data_mismatch, env_bug, network, resource
