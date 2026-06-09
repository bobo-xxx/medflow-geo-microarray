# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with this node package.

## Project Identity

This is a node package in the IRE agentic bioinformatics workflow framework. It is a standalone git repository — one repo per node. The framework core (`ire-core`) links to this repo via its `registry.yaml`.

## Language

**English is the working language.** SKILL.md, env.yaml, code comments, commit messages, error messages, NDJSON reporting, and OpenSpec artifacts are all in English.

## Development Environment

All tools and languages are managed through **conda** as the top-level environment manager.

```bash
conda activate ire-product

# Python: use uv for package management inside conda env
conda install -c conda-forge uv
uv pip install <package>

# R: use R installed via conda (not system R)
conda install -c conda-forge r-base r-essentials
# R packages: prefer conda-forge, fall back to install.packages()
```

### Rules

1. **Always activate the conda environment**: `conda activate ire-product`
2. **Python packages**: use `uv pip install`. Lock with `uv pip freeze > requirements.txt`.
3. **R packages**: prefer `conda install -c conda-forge r-<package>`.
4. **All languages**: install via conda-forge first.
5. **No hardcoded secrets**: node packages are open-source. No API keys, tokens, or credentials.

### GitHub Token

The GitHub token lives in `.claude/settings.local.json` (never committed). Available as `$GITHUB_TOKEN`.

## Node Package Reference

This project follows the **node package v2 spec** at `openspec/specs/node-package/spec.md`. The core connection protocol is at `openspec/specs/core-connection/spec.md`.

**Quick reference:**

```
node-package@version/
├── SKILL.md                   # Frontmatter (agent contract) + body (human narrative)
├── env.yaml                   # Declarative conda/mamba environment
├── scripts/
│   ├── main.<ext>             # Single entry point (.py preferred)
│   ├── input_validation.<ext>  # Optional
│   ├── output_validation.<ext> # Optional
│   └── ...                    # Internal helpers called by main
└── references/                # Optional static assets
```

**SKILL.md key rules:**
- Frontmatter is the machine contract — agents structurally parse it
- Body sections are human/LLM narrative
- Parameters declare `bind` (upstream/config/static/framework), `type`, and optional `range`
- Exceptions declare `pattern`, `nature`, and `action` (halt/skip_with_warning/escalate)
- No `retry` action — nodes handle transient errors internally
- Reporting: NDJSON to stdout (`{"level":"info",...}` and `{"level":"result",...}`)

**Scripts key rules:**
- Single entry: `main.<ext>` (`.py` preferred, `.R` secondary)
- Multi-action: `python main.py qc --input ...` (subcommand dispatch)
- Resilient: network retries, try/fallback, clear stderr + non-zero exit
- No secrets: no API keys, tokens, or credentials

## OpenSpec Workflow

This project uses OpenSpec for specification-driven development.

| Command | Purpose |
|---------|---------|
| `/opsx:propose <name>` | Create a change with ALL artifacts |
| `/opsx:apply <name>` | Implement tasks (TDD first) |
| `/opsx:archive <name>` | Archive completed change |

### Key Rules

1. **TDD first during apply/development** — write failing test, watch it fail, write minimal code to pass. Use `vitest` or `pytest`.
2. **Change names**: kebab-case (`add-qc-subcommand`, `implement-loess-normalization`)
3. **`context` and `rules`** from `openspec instructions` are AI constraints only — never copy into artifact files.
4. **Archive naming**: `YYYY-MM-DD-<change-name>`

## Core Connection

This node connects to the IRE core via the core's `registry.yaml`:

```yaml
nodes:
  geo-microarray-processing:
    source: git
    url: https://github.com/ire-product/geo-microarray-processing.git
    versions:
      - version: "1.0.0"
        commit: "<git-sha>"
        sha256: "<package-sha256>"
```

The core's `ire sync` clones this repo, checks out the pinned commit, and verifies integrity. For local testing, place this directory at `nodes/<name>@<version>/` in a core checkout.
