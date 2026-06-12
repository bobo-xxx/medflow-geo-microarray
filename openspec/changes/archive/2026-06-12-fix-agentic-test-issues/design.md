## Architecture

Six bug fixes spanning exception handling, annotation, documentation, and spec alignment.

### Fixes

| # | File | Change |
|---|---|---|
| 1 | `fetch.R:279-287,395-402` | `report_and_classify` inside `if (!safe_write_csv)` braces |
| 2 | `SKILL.md` | Remove org.db validation promise from species section |
| 3 | `fetch.R:325-339` | Call `AnnoProbe::idmap(gpl_id, type="pipe")` in Tier 4 |
| 4 | `SKILL.md` frontmatter | `bind: upstream` → `bind: config` for subcommand, gse-id, input |
| 5 | `exceptions.R:254-273` | `report_exception_ndjson` maps action→level, `quit()` on halt |
| 6 | `exceptions.R:254-273` | retry→level:retry, skip→level:decision, halt→level:exception |

### Action→Level Mapping

```
retry              → level: "retry"
skip_with_warning  → level: "decision"
halt               → level: "exception" + quit(status)
prompt             → level: "prompt"
escalate           → level: "exception" (no quit)
```
