## ADDED Requirements

### Requirement: Detect Illumina non-normalized TXT

The system SHALL detect Illumina GenomeStudio non-normalized TXT files from TargetID or AVG_Signal column headers.

#### Scenario: Illumina TXT detected

- **WHEN** a TXT file header contains TargetID or AVG_Signal columns
- **THEN** `detect_raw_type()` SHALL return `illumina_txt`

#### Scenario: Illumina TXT not confused with Agilent FE

- **WHEN** Agilent FE columns (ProbeName, gProcessedSignal) are detected
- **THEN** `detect_raw_type()` SHALL return `agilent_1c` (Agilent checked first)

### Requirement: Process Illumina non-normalized TXT

The system SHALL process Illumina GenomeStudio non-normalized TXT files via neqc.

#### Scenario: Interleaved columns parsed

- **WHEN** TXT has alternating value/Detection Pval columns
- **THEN** odd columns SHALL be treated as expression values
- **AND** even columns SHALL be treated as detection p-values

#### Scenario: neqc applied

- **WHEN** expression and detection matrices are constructed
- **THEN** `limma::neqc()` SHALL be applied (normexp bg → offset +16 → QN → log2)
- **AND** the result SHALL be a log2-scale probe-level matrix

#### Scenario: Error for missing TXT files

- **WHEN** no `.txt` files are provided
- **THEN** status SHALL be `error` with message "No TXT files found"
