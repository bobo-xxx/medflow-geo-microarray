## Why

Real agentic framework testing across two datasets revealed 6 issues — 4 bugs, 2 design gaps.

## What Changes

**Bugs (fix)**
- **Issue 1**: `report_and_classify("disk full...")` fires unconditionally outside the `if` guard on `safe_write_csv()` return
- **Issue 2**: SKILL.md body promises Tier-1 org.db validation that was removed — stale documentation
- **Issue 3**: AnnoProbe Tier 4 annotation checks `requireNamespace` but never calls `idmap()` — dead code, GB_ACC platforms fall to probe IDs

**Design Gaps (spec clarification needed)**
- **Issue 5**: `action: halt` is metadata only — doesn't `quit()`. What should "halt" mean?
- **Issue 6**: `status: "success_matrix"` coexists with `exception` lines — no gating between result and exception streams

**Core Protocol (Issue 4)**
- Root nodes (no upstream) should use `bind: config` not `bind: upstream` — requires core protocol update

## Capabilities

### Modified Capabilities
- `exception-handling`: Clarify halt semantics, wire AnnoProbe, fix false-positive W001
- `gene-annotation`: Add AnnoProbe pipe as Tier 4 with actual `idmap()` call

## Open Questions for Spec Brainstorming
- What does `action: halt` actually do — quit? block? escalate?
- How should result status interact with exception lines — mutually exclusive or additive?
- Should `bind: upstream` vs `bind: config` be determined by the core or the node?
