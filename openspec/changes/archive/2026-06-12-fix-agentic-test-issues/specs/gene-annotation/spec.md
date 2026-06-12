## ADDED Requirements

### Requirement: AnnoProbe pipe alignment in Tier 4

The system SHALL call `AnnoProbe::idmap(gpl_id, type = "pipe")` when AnnoProbe is available, not just check `requireNamespace`.

#### Scenario: AnnoProbe provides gene symbols

- **WHEN** AnnoProbe is installed and the GPL is supported
- **THEN** `idmap()` SHALL be called with `type = "pipe"`
- **AND** the resulting probe-to-gene mapping SHALL be used for Tier 4 annotation

#### Scenario: AnnoProbe fails gracefully

- **WHEN** AnnoProbe `idmap()` throws an error
- **THEN** the system SHALL catch the error and fall through to Tier 5 (probe IDs)
