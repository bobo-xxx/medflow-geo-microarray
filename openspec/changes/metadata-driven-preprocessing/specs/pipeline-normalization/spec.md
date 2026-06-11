## ADDED Requirements

### Requirement: Apply pipeline-appropriate normalization

The system SHALL apply the normalization transform appropriate to the detected preprocessing pipeline.

#### Scenario: RMA/GCRMA/neqc/lumi — pass-through

- **WHEN** pipeline is `rma`, `gcrma`, `neqc`, or `lumi`
- **THEN** the expression matrix SHALL be passed through unchanged
- **AND** the transform field SHALL be reported as `none`

#### Scenario: MAS5/GCOS — log2 only

- **WHEN** pipeline is `mas5_gcos`
- **THEN** the system SHALL apply `log2(x + 1e-6)` to the expression matrix
- **AND** the transform field SHALL be reported as `log2(x+1e-6)`

#### Scenario: GenomeStudio average normalization — shift then log2

- **WHEN** pipeline is `genomestudio_avg`
- **AND** the expression matrix contains negative values
- **THEN** the system SHALL shift all values by subtracting the global minimum: `x <- x - min(x)`
- **AND** then apply `log2(x + 1e-6)`
- **AND** the transform field SHALL be reported as `shift+log2(x+1e-6)`

#### Scenario: GenomeStudio without negatives — log2 only

- **WHEN** pipeline is `genomestudio_avg`
- **AND** the expression matrix has no negative values
- **THEN** the system SHALL apply `log2(x + 1e-6)`
- **AND** the transform field SHALL be reported as `log2(x+1e-6)`

#### Scenario: SST-RMA/TAC — pass-through

- **WHEN** pipeline is `sst_rma`
- **THEN** the expression matrix SHALL be passed through unchanged
- **AND** the `gene_assignment` column in fData SHALL be checked for gene symbol extraction

#### Scenario: Agilent FE — pass-through

- **WHEN** pipeline is `agilent_fe`
- **THEN** the expression matrix SHALL be passed through unchanged

#### Scenario: Unknown pipeline — no transform applied

- **WHEN** pipeline is `unknown`
- **THEN** the system SHALL NOT transform the expression matrix
- **AND** report a WARNING recommending raw file processing

### Requirement: Post-normalization validation

The system SHALL validate the expression matrix after normalization to detect NaN, Inf, or extreme values introduced by the transform.

#### Scenario: Clean normalization

- **WHEN** post-normalization validation passes
- **THEN** processing continues normally

#### Scenario: NaN detected after normalization

- **WHEN** NaN values are detected in the normalized matrix
- **THEN** the system SHALL report a WARNING with the NaN count
- **AND** the transform SHALL be reported as `failed`

### Requirement: Transform reporting

The system SHALL include the applied transform and detected pipeline in the NDJSON result metadata.

#### Scenario: Transform metadata reported

- **WHEN** normalization completes
- **THEN** the NDJSON result line SHALL include:
  - `pipeline_detected`: the identified pipeline name
  - `transform_applied`: the transform that was applied
  - `qn_status`: `applied` or `not_applied`
  - `transform_warnings`: list of any warnings
