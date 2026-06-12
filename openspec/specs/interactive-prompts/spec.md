# interactive-prompts Specification

## Purpose
TBD - created by archiving change exception-handling. Update Purpose after archive.
## Requirements
### Requirement: Prompt via NDJSON

The system SHALL emit interactive prompts as structured NDJSON lines to stdout.

#### Scenario: Prompt emitted

- **WHEN** the node requires user input (e.g., after 3 failed retries)
- **THEN** a JSON line SHALL be written to stdout with `level: "prompt"`
- **AND** SHALL include `code`, `msg`, and `options` fields
- **AND** each option SHALL have `choice`, `label`, and optional `prompt` fields

### Requirement: Response via stdin

The system SHALL read a single JSON response line from stdin.

#### Scenario: Response received

- **WHEN** the user responds with a valid JSON choice
- **THEN** the system SHALL parse the choice and execute the corresponding action

#### Scenario: Timeout in interactive mode

- **WHEN** no response is received within 300 seconds
- **THEN** the system SHALL execute the default action (first option in the prompt)

#### Scenario: Headless mode

- **WHEN** stdin is not a terminal (non-interactive)
- **THEN** the system SHALL immediately execute the default action without waiting

### Requirement: Prompt options contract

All interactive prompts SHALL follow a consistent option structure.

#### Scenario: Option structure

- **WHEN** a prompt is emitted
- **THEN** each option SHALL have:
  - `choice`: machine-readable identifier (e.g., "retry_proxy")
  - `label`: human-readable description (e.g., "Retry with proxy")
  - `prompt`: (optional) input prompt for the parameter
- **AND** the first option in the list SHALL be the default

