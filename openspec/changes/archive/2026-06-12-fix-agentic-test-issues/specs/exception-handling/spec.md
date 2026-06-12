## MODIFIED Requirements

### Requirement: Structured NDJSON exception reporting

The system SHALL report exceptions at the correct NDJSON level based on action severity:

- `retry` → `level: "retry"` (diagnostic, not terminal)
- `skip_with_warning` → `level: "decision"` (skip, pipeline continues)
- `halt` → `level: "exception"` (terminal, calls `quit()`)
- `prompt` → `level: "prompt"` (awaiting user input)
- `escalate` → `level: "exception"` (escalate to human, no quit)

#### Scenario: Retry does not emit exception level

- **WHEN** a transient network error triggers a retry
- **THEN** the NDJSON line SHALL have `level: "retry"` not `level: "exception"`

#### Scenario: Halt terminates the process

- **WHEN** an exception with `action: "halt"` is reported
- **AND** `dry_run` is FALSE
- **THEN** the system SHALL call `quit(status = exit_code)`

