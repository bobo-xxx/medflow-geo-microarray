## ADDED Requirements

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

## MODIFIED Requirements

### Requirement: Report an exception as structured NDJSON

The system SHALL report all exceptions as structured NDJSON with machine-readable codes. Exception codes SHALL be emitted at the error site, not deferred to the caller.

#### Scenario: Exception reported at error site

- **WHEN** any exception condition is detected
- **THEN** the NDJSON line SHALL be emitted immediately at the detection point
- **AND** include `code`, `nature`, `action`, and `msg` fields
- **AND** the `nature` field SHALL be one of: data_insufficient, data_corrupt, data_mismatch, env_bug, network, resource
