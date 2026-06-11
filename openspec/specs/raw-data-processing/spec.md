# raw-data-processing Specification

## Purpose
TBD - created by archiving change implement-raw-data-processing. Update Purpose after archive.
## Requirements
### Requirement: Detect raw file type from filename patterns

The system SHALL detect the microarray platform from raw file extensions and headers.

#### Scenario: Affymetrix CEL detected

- **WHEN** files include `*.CEL` or `*.CEL.gz`
- **THEN** `detect_raw_type()` SHALL return `affymetrix`

#### Scenario: Illumina expression IDAT detected

- **WHEN** files include `*.idat` or `*.idat.gz` AND no `*.BPM` files present
- **THEN** `detect_raw_type()` SHALL return `illumina`

#### Scenario: Methylation detected (IDAT + BPM)

- **WHEN** files include both `*.BPM` and `*.idat`
- **THEN** `detect_raw_type()` SHALL return `methylation`

#### Scenario: Agilent two-color GPR detected

- **WHEN** files include `*.GPR` or `*.GPR.gz`
- **THEN** `detect_raw_type()` SHALL return `agilent_2c`

#### Scenario: Agilent FE single-color TXT detected

- **WHEN** files include `*.txt` with header containing ProbeName or GeneName or gTotalGeneSignal
- **THEN** `detect_raw_type()` SHALL return `agilent_1c`

#### Scenario: NimbleGen PAIR detected

- **WHEN** files include `*.PAIR` or `*.PAIR.gz`
- **THEN** `detect_raw_type()` SHALL return `nimblegen`

#### Scenario: Unknown files

- **WHEN** files match no known pattern
- **THEN** `detect_raw_type()` SHALL return `unknown`

### Requirement: Detect GPR file source from header

The system SHALL auto-detect whether a GPR file uses GenePix or Agilent format.

#### Scenario: GenePix format detected

- **WHEN** first line of GPR file starts with `ATF`
- **THEN** `detect_gpr_source()` SHALL return `genepix`

#### Scenario: Agilent format detected

- **WHEN** first line of GPR file starts with `TYPE`
- **THEN** `detect_gpr_source()` SHALL return `agilent`

#### Scenario: Default to GenePix

- **WHEN** GPR header is unrecognized
- **THEN** `detect_gpr_source()` SHALL return `genepix` as default

### Requirement: Route raw files to platform-specific processor

The system SHALL route detected platforms to the appropriate Bioconductor pipeline.

#### Scenario: Route to Affymetrix processor

- **WHEN** type is `affymetrix`
- **THEN** `process_raw_files()` SHALL call `process_affy()` which applies `oligo::rma()`

#### Scenario: Route to Illumina processor

- **WHEN** type is `illumina`
- **THEN** `process_raw_files()` SHALL call `process_illumina()` which applies `limma::neqc()`

#### Scenario: Skip methylation

- **WHEN** type is `methylation`
- **THEN** `process_raw_files()` SHALL return status `skipped_methylation`

#### Scenario: Error for unknown type

- **WHEN** type is `unknown`
- **THEN** `process_raw_files()` SHALL return status `error` with descriptive message

### Requirement: Shared downstream processing

The system SHALL apply the same post-processing pipeline to all raw data sources.

#### Scenario: Validate after processing

- **WHEN** raw files are processed into an expression matrix
- **THEN** the matrix SHALL be validated with `validate_expr_matrix()`

#### Scenario: Annotate with GPL table

- **WHEN** raw files produce a probe-level matrix
- **THEN** the system SHALL attempt probe-to-gene annotation via GPL annotation table

