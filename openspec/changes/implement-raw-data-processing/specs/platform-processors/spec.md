## ADDED Requirements

### Requirement: Process Affymetrix CEL files with RMA

The system SHALL process Affymetrix CEL files using the standard RMA pipeline.

#### Scenario: CEL files processed with RMA

- **WHEN** CEL files are provided
- **THEN** `process_affy()` SHALL call `oligo::read.celfiles()` then `oligo::rma()`
- **AND** return a log2-scale probe-level expression matrix

#### Scenario: oligo package not available

- **WHEN** the `oligo` package is not installed
- **THEN** `process_affy()` SHALL return status `error` with message

### Requirement: Process Illumina IDAT files with neqc

The system SHALL process Illumina expression IDAT files using the standard neqc pipeline.

#### Scenario: IDAT files processed with neqc

- **WHEN** expression IDAT files are provided (no BPM)
- **THEN** `process_illumina()` SHALL call `limma::read.idat()` then `limma::neqc()`
- **AND** return a log2-scale probe-level expression matrix

#### Scenario: limma package not available

- **WHEN** the `limma` package is not installed
- **THEN** `process_illumina()` SHALL return status `error` with message

### Requirement: Process Agilent two-color GPR files with limma

The system SHALL process Agilent two-color GPR files with normexp background correction, loess within-array normalization, and quantile between-array normalization.

#### Scenario: GPR files processed with normexp+loess+quantile

- **WHEN** GPR files are provided
- **THEN** `process_agilent_2c()` SHALL:
  1. `detect_gpr_source()` to determine GenePix or Agilent format
  2. `limma::read.maimages()` with detected source
  3. `limma::backgroundCorrect(method="normexp")`
  4. `limma::normalizeWithinArrays(method="loess")`
  5. `limma::normalizeBetweenArrays(method="quantile")`
- **AND** return a log2-scale expression matrix

### Requirement: Process Agilent single-color FE TXT files

The system SHALL process Agilent Feature Extraction single-color TXT files.

#### Scenario: Agilent FE TXT processed with green-only

- **WHEN** Agilent FE TXT files with ProbeName/GeneName columns are provided
- **THEN** `process_agilent_1c()` SHALL call `limma::read.maimages(source="agilent", green.only=TRUE)` then `limma::normalizeBetweenArrays(method="quantile")`
- **AND** return a log2-scale expression matrix

### Requirement: Process NimbleGen PAIR files with RMA

The system SHALL process NimbleGen PAIR files using the standard RMA pipeline.

#### Scenario: PAIR files processed with RMA

- **WHEN** PAIR files are provided
- **THEN** `process_nimblegen()` SHALL call `oligo::read.xys()` then `oligo::rma()`
- **AND** return a log2-scale probe-level expression matrix

### Requirement: All processors converge to same output format

The system SHALL ensure all platform processors return the same output structure.

#### Scenario: Consistent processor output

- **WHEN** any processor completes successfully
- **THEN** the return value SHALL be a list with:
  - `status`: `"success"`
  - `expr_matrix`: log2-scale probe-level numeric matrix
  - `platform`: character platform identifier
  - `pipeline`: character description of the pipeline applied
