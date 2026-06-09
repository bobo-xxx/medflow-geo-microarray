# Node Package

## Purpose

Defines the standard package format for an IRE node. Every node — standard (analysis) or bridge (data transformation) — uses the same structure. The SKILL.md frontmatter is the machine-readable contract between the node and AI agents.

## Requirements

### Requirement: Package Format

The node SHALL be a directory `node-name@version/` containing:

```
node-name@version/
├── SKILL.md                   # Frontmatter (agent contract) + body (human narrative)
├── env.yaml                   # Declarative conda/mamba environment
├── scripts/
│   ├── main.<ext>             # Single entry point (.py preferred)
│   ├── input_validation.<ext>  # Optional: executable input checks
│   ├── output_validation.<ext> # Optional: executable output checks
│   └── ...                    # Other helper scripts called by main
└── references/                # Optional: static assets, reference data
```

#### Scenario: Agent invocation

- **WHEN** an agent invokes this node
- **THEN** it SHALL call only `scripts/main.<ext>` with the entry from SKILL.md
- **AND** for multi-action nodes, the first positional argument is the subcommand

### Requirement: SKILL.md Frontmatter Contract

The YAML frontmatter of SKILL.md SHALL be the single machine-readable source of truth. AI agents structurally parse the frontmatter; they semantically read the body sections for context.

```yaml
---
name: <short-name>
description: <function + preconditions. Primary basis for agent selection.>
type: standard          # standard | bridge
spec-ref: null          # bridge nodes only: M2M spec reference

inputs:
  - name: <file-name>
    format: <pickle|csv|tsv|h5ad|...>
    semantic_type: <type-token>
    description: <plain-language>

outputs:
  - name: <file-name>
    format: <pickle|csv|tsv|h5ad|...>
    semantic_type: <type-token>
    columns: [<expected-columns>]    # optional

entry: scripts/main.<ext>

parameters:
  - name: --<flag>
    type: file|file_out|int|float|bool|string|choice
    required: true|false
    default: <value>               # static params only
    range: [<min>, <max>]          # float/int params, optional
    bind: upstream|config|static|framework
    description: <plain-language>

exceptions:
  - exit_code: <N>
    pattern: "<stderr-substring>"
    nature: data_insufficient|data_corrupt|data_mismatch|env_bug
    action: halt|skip_with_warning|escalate

hardware:
  memory_gb: <N>
  cpu: <N>
  gpu: true|false
  runtime: "<estimate>"
---
```

#### Scenario: Bind annotations

- **WHEN** a parameter has `bind: upstream` — it SHALL be wired from a prior node's output
- **WHEN** `bind: config` — it SHALL come from the protocol's `config:` block
- **WHEN** `bind: static` — it SHALL use the `default` value; omit from CLI unless overridden
- **WHEN** `bind: framework` — it SHALL be set by the orchestrator (e.g., `--outdir`)

#### Scenario: Exception actions

- **WHEN** a node fails and stderr matches a declared `pattern` with `action: halt` — stop the run
- **WHEN** `action: skip_with_warning` — mark skipped, record reason, continue pipeline
- **WHEN** `action: escalate` — pause run, present stderr to human

No `retry` action exists. Nodes SHALL handle transient errors internally.

### Requirement: Reporting

All progress and results SHALL be written to stdout as NDJSON:

```json
{"level": "info", "msg": "<progress>"}
{"level": "result", "files": [{"path": "...", "rows": 100}], "summary": "<...>"}
```

Gate nodes additionally include `"decision": "pass"|"caution"|"rerun"|"veto"` or `"metrics": {"<name>": <number>}` in the result line.

### Requirement: Script Resilience

Scripts SHALL handle network retries and algorithm fallbacks internally. The agent never writes replacement analysis code. Nodes SHALL contain no hardcoded secrets (API keys, tokens, credentials). Preferred language order: Python > R > shell.
