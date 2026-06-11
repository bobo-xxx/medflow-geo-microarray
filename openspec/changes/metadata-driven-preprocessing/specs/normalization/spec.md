## MODIFIED Requirements

### Requirement: Normalize expression matrix

The system SHALL normalize expression matrices using pipeline-driven transforms instead of heuristic type detection. `detect_expr_type()` is retained for QC reporting and consistency verification only — it no longer drives the transform decision.

#### Scenario: Pipeline-driven normalization

- **WHEN** `normalize_expr_matrix()` is called
- **THEN** the system SHALL delegate to `apply_pipeline_transform()` for the transform decision
- **AND** use `detect_expr_type()` only for reporting and consistency verification

#### Scenario: Shift preserves expression ordering

- **WHEN** a shift is applied for negative values in GenomeStudio data
- **THEN** the shift SHALL be `x - min(x)` to preserve relative differences between all values
- **AND** no values SHALL be collapsed to a single floor

#### Scenario: log2 pseudo-count

- **WHEN** log2 transform is applied
- **THEN** the pseudo-count SHALL be `1e-6` (not 1)
- **AND** the formula SHALL be `log2(x + 1e-6)`
