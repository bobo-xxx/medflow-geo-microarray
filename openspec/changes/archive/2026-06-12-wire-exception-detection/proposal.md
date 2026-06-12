## Why

`detect_exception()` and `report_exception_ndjson()` are defined and tested but never called at actual error sites. Only 1 of 12 SKILL.md exception patterns actually fires (E801_ENV_PKG). The other 11 are dead declarations. Every error-producing operation must go through detect → classify → report.

## What Changes

- **Modified**: `node/scripts/exceptions.R` — wire detection calls at every error site
- **Modified**: `node/scripts/fetch.R` — wrap Tier 2-5 error paths
- **Modified**: `node/scripts/raw.R` — wrap processor errors
- **Modified**: `node/scripts/qc.R` — wrap qc file errors
- **Modified**: `node/scripts/clean.R` — wrap clean file errors
- **Modified**: `openspec/specs/exception-handling/spec.md` — add requirement for error site wiring

## Capabilities

### Modified Capabilities
- `exception-handling`: Add requirement that every error site SHALL call detect_exception() → report_exception_ndjson()
