## ADDED Requirements

### Requirement: Report annotation tier in result metadata

The system SHALL report which gene annotation tier produced the gene symbols in the result NDJSON metadata.

#### Scenario: Tier 1 success reported

- **WHEN** gene symbols are obtained from fData() direct column
- **THEN** `annotation_tier` SHALL be 1
- **AND** `annotation_method` SHALL be `"fData:<column_name>"`

#### Scenario: Tier 5 fallback reported

- **WHEN** all other tiers fail and probe IDs are used
- **THEN** `annotation_tier` SHALL be 5
- **AND** `annotation_method` SHALL be `"probe_ids"`
- **AND** `annotation_warning` SHALL contain the failure reasons for Tiers 1-4

#### Scenario: Raw data processing reports deferred annotation

- **WHEN** data comes from raw files (Tier 4 fetch) without an ExpressionSet
- **THEN** `annotation_tier` SHALL be NA
- **AND** `annotation_warning` SHALL indicate annotation is deferred
