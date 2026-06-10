# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with this node package.

## Project Identity

This is a node package in the IRE agentic bioinformatics workflow framework. It is a standalone git repository — one repo per node. The framework core (`ire-core`) links to this repo via its `registry.yaml`.

**This node fetches and processes GEO microarray expression data** with a 5-tier fallback strategy (local cache → series matrix → supplementary matrix → raw CEL files → metadata-only). It supports Affymetrix, Agilent, and Illumina platforms, producing probe-level and gene-level expression matrices. The original implementation was named `geo-microarray-fetch`; the name `geo-microarray-processing` broadens the scope to include downstream cleaning in addition to data retrieval.

### Current Status: Pre-implementation

The node package files (`SKILL.md`, `env.yaml`, `scripts/`, `references/`) have **not yet been created**. The reference implementation lives in `original/geo-microarray-fetch.zip` — a complete, working R implementation that must be adapted to the node-package v2 format.

## Reference Implementation

The `original/geo-microarray-fetch.zip` contains the working R code this node is based on:

```
geo-microarray-fetch/
├── SKILL.md                   # Original skill definition (NOT node-package v2 format)
├── scripts/
│   ├── fetch_geo.R            # Main script: 5-tier fallback, GEO data retrieval (~26 KB)
│   ├── platform_detect.R      # Platform detection: Affymetrix/Agilent/Illumina/methylation
│   └── validate.R             # Validation: CEL integrity, expression matrix QC
└── references/
    ├── ERROR_CODES.md         # E001–E007 error codes with severity and recovery
    └── PLATFORMS.md           # File pattern detection rules for each platform
```

Unzip with `unzip original/geo-microarray-fetch.zip -d /tmp/geo-ref` when you need to consult the original implementation.

**GEO reference**: [GEO FTP README](https://ftp.ncbi.nlm.nih.gov/geo/README.txt) — GEO data formats, file conventions, and submission guidelines.

### Architecture: 5-Tier Fallback Strategy

The original `fetch_geo.R` implements this decision flow:

```
1. Local data?     → use cached files, skip download
2. Series matrix?  → GEOquery::getGEO(GSEMatrix=TRUE), normalized, fast
3. Suppl matrix?   → author-provided *.txt.gz from suppl/ directory
4. Raw CEL files?  → oligo::rma() normalization for Affymetrix, slower
5. Metadata only?  → when all expression data retrieval fails
```

### Platform Detection Rules (from original)

| File Pattern | Platform | Processing |
|---|---|---|
| `*.CEL(.gz)?` | Affymetrix | `oligo::rma()` |
| `*.GPR(.gz)?` | Agilent | `limma::read.maimages()` |
| `*.idat` only | Illumina | `illuminaio::readIDAT()` |
| `*.idat` + `*.BPM` | Methylation | **SKIP** — not expression data |
| `series_matrix.txt.gz` | TXT/CSV | Direct import |

### Original R Dependencies

- **Bioconductor**: `GEOquery`, `Biobase`, `oligo` (optional, for raw CEL)
- **CRAN**: `dplyr`, `tidyr`, `stringr`, `tibble`

These inform what goes into `env.yaml`.

## Language

**English is the working language.** SKILL.md, env.yaml, code comments, commit messages, error messages, NDJSON reporting, and OpenSpec artifacts are all in English.

## Development Environment

This node has its own conda environment, declared in `env.yaml`. The framework's Execution Engine creates it automatically — during development, the node author creates and tests it manually.

**Dual R version testing**: `env.yaml` targets R 4.3 / Bioconductor 3.18 (production). `env-4.5.yaml` targets R 4.5 / Bioconductor 3.20 (forward-compat). CI tests both.

```bash
# Create and activate the node-specific environment
conda env create -f env.yaml
conda activate geo-microarray-processing      # name from env.yaml

# Forward-compat testing
conda env create -f env-4.5.yaml
conda activate geo-microarray-processing-4.5

# Run tests in either env
Rscript -e 'testthat::test_dir("tests/testthat/")'

# Python: use uv for package management inside the conda env
conda install -c conda-forge uv
uv pip install <package>

# R: use R installed via conda (not system R)
conda install -c conda-forge r-base r-essentials
# R packages: prefer conda-forge, fall back to install.packages()
```

### Conda Mirror Configuration (Development Only)

For faster package installation during development, configure conda to use Chinese mirror sources. **Never commit mirror configs to env.yaml** — they are development-local only.

```bash
# Add mirror channels to conda config (one-time setup)
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge
conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda

# Or set via environment variable for one-off use
export CONDA_CHANNELS="https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge,https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda"

# Show current channel priority
conda config --show channels
```

Channel priority: mirrors first (fast), then upstream conda-forge/bioconda (authoritative). The `env.yaml` file always uses upstream channels only — mirrors are a local dev convenience.

### Proxy for R Package Installation

Use proxy `http://localhost:2999` to speed up CRAN package installation:

```bash
# Set proxy before installing R packages
export http_proxy=http://localhost:2999
export https_proxy=http://localhost:2999

# Or set within R
Rscript -e 'Sys.setenv(http_proxy="http://localhost:2999", https_proxy="http://localhost:2999")'
```

### AnnoProbe (CRAN Only)

`r-annoprobe` is NOT available on conda channels. Install from CRAN after creating the conda environment:

```bash
conda activate geo-microarray-processing
Rscript -e 'install.packages("AnnoProbe", repos = "https://cloud.r-project.org")'
```

### Rules

1. **Use the node's own env**: create from `env.yaml`, activate it for all development work.
2. **Python packages**: use `uv pip install` inside the node's conda env.
3. **R packages**: prefer `conda install -c conda-forge r-<package>`. Fall back to `install.packages()` for CRAN-only packages (e.g., AnnoProbe).
4. **All languages**: install via conda-forge first.
5. **No hardcoded secrets**: node packages are open-source. No API keys, tokens, or credentials.

### GitHub Token

The GitHub token lives in `.claude/settings.local.json` (never committed). Available as `$GITHUB_TOKEN`.

### R Development

This node is primarily R. Key commands once the package is created:

```bash
# Run the node
Rscript scripts/main.R fetch --gse-id GSE100155 --outdir ./output

# Run with subcommands (when multi-action)
Rscript scripts/main.R qc --input ./output/expr_gene_GSE100155.csv
Rscript scripts/main.R clean --input ./output/expr_gene_GSE100155.csv

# Run a single test
Rscript -e 'testthat::test_file("tests/testthat/test-fetch.R")'

# Run all tests
Rscript -e 'testthat::test_dir("tests/testthat/")'

# Install R packages for the node env
conda activate <env-name>
conda install -c conda-forge r-testthat r-lintr
Rscript -e 'BiocManager::install(c("GEOquery", "Biobase", "oligo"))'
```

**Testing**: Use `testthat` for R. TDD workflow: write test → `Rscript -e 'testthat::test_file("tests/testthat/test-<name>.R")'` → write minimal code → repeat.

**Linting**: `Rscript -e 'lintr::lint_dir("scripts/")'`

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

1. **TDD first during apply/development** — write failing test, watch it fail, write minimal code to pass. Use `testthat` for R (this node) or `pytest` for Python nodes.
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

## Migration: What Needs to Be Built

The original `geo-microarray-fetch` reference implementation must be adapted to the node-package v2 format. Key differences from the original:

| Original | Node Package v2 |
|---|---|
| SKILL.md is a skill definition | SKILL.md frontmatter is machine-readable agent contract with `inputs`, `outputs`, `parameters`, `exceptions`, `hardware` |
| Multi-file scripts with `source()` | Single-entry `scripts/main.R` with subcommand dispatch |
| No env.yaml | Declarative `env.yaml` (conda/mamba) |
| Returns R list | Reports NDJSON to stdout (`{"level":"info",...}`, `{"level":"result",...}`) |
| Hardcoded proxy `localhost:1086` | No hardcoded config; proxy via `parameters` with `bind: config` |
| Status codes as strings | Structured `exceptions` with `pattern`, `nature`, `action` |

**Build order** (use `/opsx:new` for each):
1. `init-node-package` — create `SKILL.md` frontmatter, `env.yaml`, `scripts/main.R` skeleton with fetch subcommand
2. `port-fetch-geo` — adapt `fetch_geo.R` logic into `main.R`, convert to NDJSON reporting
3. `add-qc-subcommand` — add QC subcommand for expression matrix validation
4. `add-clean-subcommand` — downstream cleaning/normalization (the "processing" scope)
