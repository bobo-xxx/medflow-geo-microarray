### Task 1: Fix Issue 1 — Conditional report_and_classify

- [x] Move `report_and_classify("disk full...")` inside `if (!safe_write_csv)` brace block (2 sites)
- [x] Test: safe_write_csv returns TRUE on success

### Task 2: Fix Issue 2 — SKILL.md species text

- [x] Remove org.db validation promise from species section body text

### Task 3: Fix Issue 3 — Wire AnnoProbe idmap()

- [x] Call `AnnoProbe::idmap(gpl_id, type = "pipe")` in Tier 4 annotation
- [x] Build probe2gene data.frame from AnnoProbe result
- [x] tryCatch for graceful fallback on error

### Task 4: Fix Issue 4 — SKILL.md bind values

- [x] Change `bind: upstream` → `bind: config` for subcommand, --gse-id, --input

### Task 5: Fix Issues 5+6 — Action→level mapping + halt semantics

- [x] `report_exception_ndjson()` maps action to correct NDJSON level
- [x] `action: "halt"` calls `quit(status = exit_code)` (suppressed via `dry_run=TRUE`)
- [x] Tests: halt quit behavior, dry_run suppression, level mapping for all actions

### Task 6: Add regression tests

- [x] Issue 1: safe_write_csv success path
- [x] Issue 5: halt quit + dry_run
- [x] Issue 6: action→level mapping (retry, decision, exception, escalate)
- [x] 278 tests pass (0 FAIL)
