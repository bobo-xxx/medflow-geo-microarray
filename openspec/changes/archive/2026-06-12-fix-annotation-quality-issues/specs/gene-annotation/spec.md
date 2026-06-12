## MODIFIED Requirements

### Requirement: Report annotation tier in result metadata

The system SHALL clean gene symbol prefixes (`entg|`, `ref|`, `gb|`) from all annotation sources.

#### Scenario: entg| prefix cleaned at Tier 3

- **WHEN** GPL table returns gene symbols with `entg|` prefix
- **THEN** `clean_gene_symbols()` SHALL strip the prefix before aggregation

#### Scenario: entg| prefix cleaned at Tier 4

- **WHEN** BioC annotation DB returns gene symbols with `entg|` prefix
- **THEN** `clean_gene_symbols()` SHALL strip the prefix before aggregation

## ADDED Requirements

### Requirement: GPL supplementary file fallback

The system SHALL attempt to download the GPL supplementary annotation file when the GPL table lacks gene symbol columns.

#### Scenario: Supplementary file provides gene symbols

- **WHEN** `Meta(gpl)$supplementary_file` references a valid annotation file
- **AND** the file is under 1GB
- **THEN** the system SHALL download and parse the file for gene symbol mappings

#### Scenario: Supplementary file unavailable or too large

- **WHEN** the supplementary file is over 1GB or download fails
- **THEN** the system SHALL fall through to the next annotation tier

### Requirement: Gene matrix suppressed at Tier 5

The system SHALL NOT produce a gene-level expression matrix when annotation falls to Tier 5.

#### Scenario: Tier 5 skips gene CSV

- **WHEN** `annotation_tier == 5` (probe IDs as gene names)
- **THEN** the gene CSV SHALL NOT be written
- **AND** only the probe matrix SHALL be produced as authoritative output
