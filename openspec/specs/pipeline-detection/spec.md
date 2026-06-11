# pipeline-detection Specification

## Purpose
TBD - created by archiving change metadata-driven-preprocessing. Update Purpose after archive.
## Requirements
### Requirement: Parse GEO data_processing metadata

The system SHALL parse the `data_processing` column(s) from the ExpressionSet's pData to identify the preprocessing pipeline applied by the submitter.

Keywords are matched case-insensitively in priority order. The first matching pipeline wins.

#### Scenario: RMA pipeline detected

- **WHEN** `data_processing` contains "RMA" or "robust multi-array"
- **THEN** the system SHALL identify the pipeline as `rma`

#### Scenario: MAS5/GCOS pipeline detected

- **WHEN** `data_processing` contains "GCOS" or "MAS5" or "target intensity"
- **THEN** the system SHALL identify the pipeline as `mas5_gcos`

#### Scenario: GCRMA pipeline detected

- **WHEN** `data_processing` contains "GCRMA" or "gc content" or "gc-rma"
- **THEN** the system SHALL identify the pipeline as `gcrma`

#### Scenario: SST-RMA/TAC pipeline detected

- **WHEN** `data_processing` contains "TAC" or "SST-RMA" or "signal space" or "Transcriptome Analysis Console"
- **THEN** the system SHALL identify the pipeline as `sst_rma`

#### Scenario: GenomeStudio average normalization detected

- **WHEN** `data_processing` contains "Genome Studio" or "GenomeStudio" or "average normalization" AND does NOT contain "neqc" or "normexp"
- **THEN** the system SHALL identify the pipeline as `genomestudio_avg`

#### Scenario: neqc pipeline detected

- **WHEN** `data_processing` contains "neqc" or "normexp" or "offset 16"
- **THEN** the system SHALL identify the pipeline as `neqc`

#### Scenario: lumi pipeline detected

- **WHEN** `data_processing` contains "lumi" or "VST" or "variance stabiliz"
- **THEN** the system SHALL identify the pipeline as `lumi`

#### Scenario: Agilent Feature Extraction detected

- **WHEN** `data_processing` contains "Feature Extraction" or "Agilent" with "loess"
- **THEN** the system SHALL identify the pipeline as `agilent_fe`

#### Scenario: No pipeline identified

- **WHEN** `data_processing` does not match any known pipeline keywords
- **THEN** the system SHALL return `unknown`

### Requirement: Detect quantile normalization

The system SHALL detect whether the expression matrix has already been quantile-normalized using a 5-percentile coefficient of variation check.

#### Scenario: QN detected

- **WHEN** CV of sample quantiles at p25, p50, p75, p90, p95 is less than 0.002
- **THEN** the system SHALL report QN status as `applied`

#### Scenario: QN not detected

- **WHEN** CV of sample quantiles at any of p25, p50, p75, p90, p95 exceeds 0.002
- **THEN** the system SHALL report QN status as `not_applied`

### Requirement: Verify metadata-to-data consistency

The system SHALL verify that the detected pipeline is consistent with the actual data distribution.

#### Scenario: RMA claim consistent with log-scale data

- **WHEN** pipeline is `rma`, `gcrma`, `neqc`, or `lumi`
- **AND** `detect_expr_type()` returns "log"
- **THEN** the system SHALL accept the pipeline claim

#### Scenario: MAS5/GCOS claim consistent with linear data

- **WHEN** pipeline is `mas5_gcos`
- **AND** `detect_expr_type()` returns "raw" and min >= 0
- **THEN** the system SHALL accept the pipeline claim

#### Scenario: GenomeStudio avg claim consistent with linear data with negatives

- **WHEN** pipeline is `genomestudio_avg`
- **AND** `detect_expr_type()` returns "raw" and min < 0
- **THEN** the system SHALL accept the pipeline claim

#### Scenario: Pipeline claim inconsistent with data

- **WHEN** pipeline claim does not match data distribution
- **THEN** the system SHALL report a WARNING in NDJSON
- **AND** fall back to `unknown` pipeline (raw file processing)

