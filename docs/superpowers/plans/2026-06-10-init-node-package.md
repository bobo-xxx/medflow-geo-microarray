# init-node-package Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create the node package skeleton — SKILL.md frontmatter + body, env.yaml (R 4.3 + R 4.5), scripts/main.R subcommand dispatch, references/, and test infrastructure. All downstream changes (`port-fetch-geo`, `add-qc-subcommand`, `add-clean-subcommand`) depend on this.

**Architecture:** `scripts/main.R` is the single entry point per node-package v2 spec. It parses CLI args, dispatches to subcommand handlers, and delegates NDJSON output to `scripts/report.R`. The subcommand handlers (`fetch`, `qc`, `clean`) are stubs in this change — they echo their name and exit. SKILL.md frontmatter is the machine-readable agent contract; the body is human/LLM narrative following the skill-like pattern (11 sections).

**Tech Stack:** R >= 4.3, YAML (SKILL.md frontmatter), testthat, conda/mamba, AnnoProbe (CRAN) for probe annotation

**Test Datasets:** GSE318047 (human GPL570), GSE156508 (human GPL16686 no-gene-symbol), GSE11381 (mouse GPL339), GSE4105 (rat GPL85), GSE84422 (multi-platform GPL96+97+570), GSE42861 (methylation skip)

---

## File Structure

| File | Purpose |
|---|---|
| `SKILL.md` | Frontmatter (YAML agent contract) + body (11-section narrative) |
| `env.yaml` | Conda environment: R 4.3, Bioconductor 3.18 |
| `env-4.5.yaml` | Conda environment: R 4.5, Bioconductor 3.20 |
| `scripts/main.R` | Single entry point: arg parsing, subcommand dispatch, NDJSON reporting |
| `scripts/report.R` | NDJSON helper: `report_info()`, `report_result()`, `report_error()` |
| `references/ERROR_CODES.md` | E001–E007 error codes with severity and recovery (from original) |
| `references/PLATFORMS.md` | Platform detection rules (from original) |
| `tests/testthat/test-main.R` | Integration tests for CLI dispatch and NDJSON output |
| `tests/testthat/test-skillmd.R` | Structural validation of SKILL.md frontmatter |
| `tests/testthat/helpers.R` | Shared test utilities |
| `tests/fixtures/` | Empty directory placeholder (fixtures added in port-fetch-geo) |

---

### Task 1: Create references/ directory from original

**Files:**
- Create: `references/ERROR_CODES.md`
- Create: `references/PLATFORMS.md`

- [ ] **Step 1: Create references/ERROR_CODES.md**

```bash
mkdir -p references
```

Copy from original zip's `references/ERROR_CODES.md` (the entire file as-is, since it's a static reference document):

```bash
unzip -p original/geo-microarray-fetch.zip geo-microarray-fetch/references/ERROR_CODES.md > references/ERROR_CODES.md
```

- [ ] **Step 2: Create references/PLATFORMS.md**

```bash
unzip -p original/geo-microarray-fetch.zip geo-microarray-fetch/references/PLATFORMS.md > references/PLATFORMS.md
```

- [ ] **Step 3: Verify files exist and have content**

```bash
wc -l references/ERROR_CODES.md references/PLATFORMS.md
```

Expected: ERROR_CODES.md ~67 lines, PLATFORMS.md ~37 lines

- [ ] **Step 4: Commit**

```bash
git add references/ERROR_CODES.md references/PLATFORMS.md
git commit -m "feat: add ERROR_CODES.md and PLATFORMS.md references from original

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Create testthat infrastructure and test-skillmd.R

**Files:**
- Create: `tests/testthat/test-skillmd.R`
- Create: `tests/testthat/helpers.R`

- [ ] **Step 1: Create directory structure**

```bash
mkdir -p tests/testthat
```

- [ ] **Step 2: Write the failing test for SKILL.md frontmatter validation**

Create `tests/testthat/test-skillmd.R`:

```r
library(testthat)

# helpers
source("tests/testthat/helpers.R")

describe("SKILL.md", {

  it("exists at the package root", {
    expect_true(file.exists("SKILL.md"))
  })

  it("has valid YAML frontmatter", {
    content <- readLines("SKILL.md", warn = FALSE)

    # Find YAML frontmatter delimiters
    dashes <- which(content == "---")
    expect_gte(length(dashes), 2,
      "SKILL.md must have opening and closing '---' YAML delimiters")

    yaml_lines <- content[(dashes[1] + 1):(dashes[2] - 1)]
    expect_gt(length(yaml_lines), 0,
      "YAML frontmatter must not be empty")

    # Write YAML to temp file and parse with yaml package
    tmp <- tempfile(fileext = ".yaml")
    writeLines(yaml_lines, tmp)
    frontmatter <- yaml::read_yaml(tmp)

    expect_type(frontmatter, "list")
  })

  it("declares required frontmatter fields", {
    content <- readLines("SKILL.md", warn = FALSE)
    dashes <- which(content == "---")
    yaml_lines <- content[(dashes[1] + 1):(dashes[2] - 1)]
    tmp <- tempfile(fileext = ".yaml")
    writeLines(yaml_lines, tmp)
    fm <- yaml::read_yaml(tmp)

    required <- c("name", "description", "type", "inputs", "outputs",
                  "entry", "parameters", "exceptions", "hardware")
    for (field in required) {
      expect_true(field %in% names(fm),
        sprintf("SKILL.md frontmatter must have '%s' field", field))
    }
  })

  it("has type 'standard'", {
    fm <- parse_skillmd_frontmatter("SKILL.md")
    expect_equal(fm$type, "standard")
  })

  it("has entry pointing to scripts/main.R", {
    fm <- parse_skillmd_frontmatter("SKILL.md")
    expect_equal(fm$entry, "scripts/main.R")
  })

  it("has parameters for subcommand, gse-id, outdir, proxy, api-key", {
    fm <- parse_skillmd_frontmatter("SKILL.md")
    param_names <- vapply(fm$parameters, `[[`, "", "name")
    expected <- c("subcommand", "--gse-id", "--outdir", "--input",
                  "--proxy", "--api-key")
    for (p in expected) {
      expect_true(p %in% param_names,
        sprintf("SKILL.md must declare parameter '%s'", p))
    }
  })

  it("has 6 exception patterns", {
    fm <- parse_skillmd_frontmatter("SKILL.md")
    expect_gte(length(fm$exceptions), 6)
  })
})
```

- [ ] **Step 3: Write helpers.R**

Create `tests/testthat/helpers.R`:

```r
#' Parse SKILL.md frontmatter YAML
#' @param path Path to SKILL.md
#' @return List of frontmatter fields
parse_skillmd_frontmatter <- function(path = "SKILL.md") {
  content <- readLines(path, warn = FALSE)
  dashes <- which(content == "---")
  yaml_lines <- content[(dashes[1] + 1):(dashes[2] - 1)]
  tmp <- tempfile(fileext = ".yaml")
  writeLines(yaml_lines, tmp)
  yaml::read_yaml(tmp)
}
```

- [ ] **Step 4: Run test, verify it FAILS (no SKILL.md yet)**

```bash
conda activate geo-microarray-processing 2>/dev/null || true
Rscript -e 'testthat::test_file("tests/testthat/test-skillmd.R")'
```

Expected: FAIL — "SKILL.md does not exist" or "file.exists('SKILL.md') is not TRUE"

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/test-skillmd.R tests/testthat/helpers.R
git commit -m "test: add SKILL.md frontmatter validation tests

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: Create SKILL.md

**Files:**
- Create: `SKILL.md`

- [ ] **Step 1: Write SKILL.md with full frontmatter + body**

Create `SKILL.md`:

```markdown
---
name: geo-microarray-processing
description: >
  Fetch and process GEO microarray expression data with 5-tier fallback
  (local cache → series matrix → supplementary → raw CEL → metadata-only).
  Supports Affymetrix, Agilent, and Illumina platforms with automatic
  platform detection. Produces probe-level and gene-level expression
  matrices with normalization, probe-to-gene aggregation, and quality
  validation. Primary species: human (9606), mouse (10090), rat (10116).
type: standard

inputs: []

outputs:
  - name: expr_probe_{gse_id}_{gpl}.csv
    format: csv
    semantic_type: probe_expression_matrix
    description: >
      Probe-level expression matrix (probes as rows, samples as columns).
      One file per platform when multiple GPLs are present.

  - name: expr_gene_{gse_id}_{gpl}.csv
    format: csv
    semantic_type: gene_expression_matrix
    description: >
      Gene-level expression matrix aggregated from probes via mean.
      Gene symbols as rows, samples as columns. One file per platform.

entry: scripts/main.R

parameters:
  - name: subcommand
    type: choice
    choices: [fetch, qc, clean]
    required: true
    bind: upstream
    description: Operation to perform (fetch data, quality check, or clean/normalize)

  - name: --gse-id
    type: string
    required: true
    bind: upstream
    description: GEO Series identifier (e.g., GSE100155)

  - name: --outdir
    type: file_out
    required: false
    default: .
    bind: framework
    description: Output directory for expression matrix CSV files

  - name: --input
    type: file
    required: false
    bind: upstream
    description: Input expression matrix CSV (required for qc and clean subcommands)

  - name: --proxy
    type: string
    required: false
    bind: config
    description: HTTP/HTTPS proxy URL for GEO network access

  - name: --api-key
    type: string
    required: false
    bind: config
    description: NCBI API key for higher rate limits (10 req/s vs 3 req/s)

exceptions:
  - exit_code: 1
    pattern: "No suppl directory"
    nature: data_insufficient
    action: skip_with_warning

  - exit_code: 1
    pattern: "Methylation BPM present"
    nature: data_mismatch
    action: skip_with_warning

  - exit_code: 1
    pattern: "No series matrix"
    nature: data_insufficient
    action: skip_with_warning

  - exit_code: 1
    pattern: "No raw files"
    nature: data_insufficient
    action: skip_with_warning

  - exit_code: 1
    pattern: "Corrupted CEL file"
    nature: data_corrupt
    action: skip_with_warning

  - exit_code: 2
    pattern: "All data retrieval methods failed"
    nature: data_insufficient
    action: halt

hardware:
  memory_gb: 4
  cpu: 2
  gpu: false
  runtime: "~5–30 minutes depending on dataset size and fallback tier"
---

# Overview

This node fetches and processes microarray expression data from the
Gene Expression Omnibus (GEO) database using R/Bioconductor. It
implements a smart five-tier fallback strategy to maximize data
retrieval success and automatically handles platform detection,
normalization, probe-to-gene aggregation, and quality validation.

## When to Use

Use this node when you need to:
- Download gene expression microarray data from GEO
- Process Affymetrix, Agilent, or Illumina arrays
- Convert probe-level to gene-level expression
- Handle missing data with automatic fallback
- Batch download multiple GSE datasets
- Validate expression matrix quality

**Trigger phrases**: "fetch GEO", "download microarray", "get GSE",
"GEO expression", "series matrix"

## Key Features

- **Five-tier fallback**: local cache → processed series matrix →
  supplementary matrix → raw CEL files → metadata-only
- **Automatic platform detection**: Affymetrix, Agilent, Illumina,
  TXT/CSV, methylation arrays
- **Methylation array skip**: BPM + IDAT files detected and excluded
- **Multi-GPL support**: per-platform output files with GPL suffix
- **Probe-to-gene aggregation**: mean-based aggregation with `///`
  separator handling
- **Smart normalization**: auto-detect raw/centered/log expression
  type and apply appropriate transform
- **Gene annotation**: multi-tier fallback (fData → GPL annotation →
  probe IDs)

## Supported Platforms

See [PLATFORMS.md](references/PLATFORMS.md) for detailed detection rules.

| Platform | File Pattern | Processing |
|---|---|---|
| Affymetrix | `*.CEL(.gz)?` | `oligo::rma()` |
| Agilent | `*.GPR(.gz)?` | `limma::read.maimages()` |
| Agilent FE | TXT with ProbeName/GeneName | `limma::read.maimages(source="agilent")` |
| Illumina | `*.idat` (no BPM) | `illuminaio::readIDAT()` |
| TXT/CSV | `series_matrix.txt.gz` | Direct import |
| Methylation | `*.idat` + `*.BPM` | **SKIP** |

## Species Support

The node detects organism from GEO metadata
(`experimentData(eset)@other$sample_taxid`) and routes to the
appropriate annotation database:

| TaxID | Species | Annotation | Tier |
|---|---|---|---|
| 9606 | *Homo sapiens* | `org.Hs.eg.db` | 1 — validated |
| 10090 | *Mus musculus* | `org.Mm.eg.db` | 1 — validated |
| 10116 | *Rattus norvegicus* | `org.Rn.eg.db` | 1 — validated |
| Other | Any | GPL table only | 2 — pass-through |

Tier-1 species get validated gene symbol mapping against the species
annotation database. Tier-2 species rely on GPL annotation table gene
symbols with a warning; if unavailable, probe IDs are used as gene symbols.

## Data Priority

```
1. Local data?     → use cached files, skip download
2. Series matrix?  → GEOquery::getGEO(GSEMatrix=TRUE), normalized, fast
3. Suppl matrix?   → author-provided *.txt.gz from suppl/ directory
4. Raw CEL files?  → oligo::rma() normalization for Affymetrix, slower
5. Metadata only?  → when all expression data retrieval fails
```

## Usage

```bash
# Fetch GEO data
Rscript scripts/main.R fetch --gse-id GSE100155 --outdir ./output

# Quality check on expression matrix
Rscript scripts/main.R qc --input ./output/GSE100155/expr_gene_GSE100155_GPL570.csv

# Clean/normalize expression matrix
Rscript scripts/main.R clean --input ./output/GSE100155/expr_gene_GSE100155_GPL570.csv

# With proxy and API key
Rscript scripts/main.R fetch --gse-id GSE100155 --proxy http://proxy:8080 --api-key abc123
```

## Parameters

### Required

- `subcommand` — Operation: `fetch`, `qc`, or `clean`. Set by the
  orchestrator via `bind: upstream`.
- `--gse-id` — GEO Series identifier. Required for `fetch`.

### Optional

- `--outdir` — Output directory. Defaults to current directory.
  Set by the framework via `bind: framework`.
- `--input` — Input expression matrix CSV. Required for `qc` and
  `clean` subcommands.
- `--proxy` — HTTP/HTTPS proxy URL for GEO network access.
  Read from config via `bind: config`.
- `--api-key` — NCBI API key for higher rate limits. Falls back to
  `NCBI_API_KEY` environment variable. Via `bind: config`.

## Output

The fetch subcommand produces per-platform files in the output directory:

```
output/
└── GSE12345/
    ├── expr_probe_GSE12345_GPL570.csv   # Probe-level (probes × samples)
    ├── expr_gene_GSE12345_GPL570.csv    # Gene-level (genes × samples)
    └── suppl/                            # Downloaded supplementary files
```

**File formats**:
- **Probe-level CSV**: Probes as rows, samples as columns
- **Gene-level CSV**: Genes as rows, samples as columns (mean-aggregated)

**NDJSON reporting** (to stdout):

Progress messages:
```json
{"level":"info","msg":"Downloading suppl files for GSE100155..."}
{"level":"info","msg":"Platform detected: Affymetrix (GPL570)"}
```

Result on success:
```json
{"level":"result","status":"success_matrix","files":[
  {"path":"output/GSE100155/expr_probe_GSE100155_GPL570.csv","rows":54675,"cols":24},
  {"path":"output/GSE100155/expr_gene_GSE100155_GPL570.csv","rows":20838,"cols":24}
],"metadata":{"platform":"GPL570","organism":"Homo sapiens","n_samples":24}}
```

Result on skip/methylation:
```json
{"level":"result","status":"skipped_methylation",
 "reason":"Methylation BPM present — not an expression array",
 "decision":"skip_with_warning"}
```

### Status Codes

| Code | Description |
|---|---|
| `success_local` | Used existing local data, skipped download |
| `success_matrix` | Successfully retrieved GEO processed series matrix |
| `success_suppl` | Author-provided processed matrix from suppl directory |
| `success_raw` | Raw CEL files processed with RMA normalization |
| `metadata_only` | Metadata only, no expression data available |
| `skipped_methylation` | Methylation array detected, excluded |
| `error` | All data retrieval methods failed |

## Error Handling

See [ERROR_CODES.md](references/ERROR_CODES.md) for detailed error
code specifications and recovery strategies.

| Exit Code | Pattern | Nature | Action |
|---|---|---|---|
| 1 | No suppl directory | `data_insufficient` | `skip_with_warning` |
| 1 | Methylation BPM present | `data_mismatch` | `skip_with_warning` |
| 1 | No series matrix | `data_insufficient` | `skip_with_warning` |
| 1 | No raw files | `data_insufficient` | `skip_with_warning` |
| 1 | Corrupted CEL file | `data_corrupt` | `skip_with_warning` |
| 2 | All data retrieval methods failed | `data_insufficient` | `halt` |

The node handles transient errors (network timeouts, rate limiting)
internally — no `retry` action is declared. Rate limiting is automatic:
3 req/s without API key, 10 req/s with key.

## Integration

This node works with:
- Downstream QC node — validate expression matrix quality
- Downstream clean node — normalization, batch correction
- Downstream analysis nodes — differential expression, enrichment

For multi-platform datasets, downstream nodes consume per-platform
CSV files individually. Cross-platform merging is deferred to a
dedicated merge node or meta-analysis pipeline.
```

- [ ] **Step 2: Run test-skillmd.R, verify it PASSES**

```bash
conda activate geo-microarray-processing 2>/dev/null || Rscript -e 'install.packages("yaml", repos="https://cloud.r-project.org")'
Rscript -e 'testthat::test_file("tests/testthat/test-skillmd.R")'
```

Expected: all 6 tests PASS

- [ ] **Step 3: Commit**

```bash
git add SKILL.md
git commit -m "feat: add SKILL.md with v2 node-package frontmatter and skill-like body

Frontmatter: 6 parameters (subcommand, gse-id, outdir, input, proxy, api-key),
6 exceptions (E001-E005 skip_with_warning, all-failed halt),
hardware profile (4GB/2CPU).
Body: 11-section narrative following skill pattern.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: Create env.yaml and env-4.5.yaml

**Files:**
- Create: `env.yaml`
- Create: `env-4.5.yaml`

- [ ] **Step 1: Write env.yaml (R 4.3 / Bioconductor 3.18)**

Create `env.yaml`:

```yaml
name: geo-microarray-processing

channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  # R base
  - r-base >=4.3,<4.4
  - r-essentials

  # Bioconductor — data access
  - bioconductor-geoquery
  - bioconductor-biobase

  # Bioconductor — raw data processing
  - bioconductor-oligo
  - bioconductor-limma

  # CRAN — data manipulation
  - r-dplyr
  - r-tidyr
  - r-stringr
  - r-tibble

  # CRAN — probe annotation (handles GB_ACC-only platforms like GPL16686)
  - r-annoprobe

  # YAML parsing (SKILL.md validation)
  - r-yaml

  # Development & testing
  - r-testthat
  - r-lintr
```

- [ ] **Step 2: Write env-4.5.yaml (R 4.5 / Bioconductor 3.20)**

Create `env-4.5.yaml`:

```yaml
name: geo-microarray-processing-4.5

channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  # R base
  - r-base >=4.5,<4.6
  - r-essentials

  # Bioconductor — data access
  - bioconductor-geoquery
  - bioconductor-biobase

  # Bioconductor — raw data processing
  - bioconductor-oligo
  - bioconductor-limma

  # CRAN — data manipulation
  - r-dplyr
  - r-tidyr
  - r-stringr
  - r-tibble

  # CRAN — probe annotation (handles GB_ACC-only platforms like GPL16686)
  - r-annoprobe

  # YAML parsing (SKILL.md validation)
  - r-yaml

  # Development & testing
  - r-testthat
  - r-lintr
```

- [ ] **Step 3: Verify YAML syntax**

```bash
python3 -c "import yaml; yaml.safe_load(open('env.yaml')); print('env.yaml OK')"
python3 -c "import yaml; yaml.safe_load(open('env-4.5.yaml')); print('env-4.5.yaml OK')"
```

Expected: both print "OK"

- [ ] **Step 4: Commit**

```bash
git add env.yaml env-4.5.yaml
git commit -m "feat: add conda environments for R 4.3 and R 4.5

env.yaml: R >=4.3,<4.4 with Bioconductor 3.18 (production)
env-4.5.yaml: R >=4.5,<4.6 with Bioconductor 3.20 (forward-compat)
Channels: conda-forge, bioconda, defaults (no mirrors committed)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: Create test-main.R (integration tests for CLI dispatch)

**Files:**
- Create: `tests/testthat/test-main.R`

- [ ] **Step 1: Write the failing integration test**

Create `tests/testthat/test-main.R`:

```r
library(testthat)

source("tests/testthat/helpers.R")

describe("main.R CLI dispatch", {

  it("prints help when no arguments given", {
    result <- system2("Rscript", c("scripts/main.R"),
                      stdout = TRUE, stderr = TRUE)
    expect_true(any(grepl("Usage|usage|fetch|subcommand",
                          result, ignore.case = TRUE)),
      "Should print usage information when no args")
  })

  it("exits with error for unknown subcommand", {
    exit_code <- system2("Rscript", c("scripts/main.R", "unknown"),
                         stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 1)
  })

  it("exits with error when fetch missing --gse-id", {
    exit_code <- system2("Rscript", c("scripts/main.R", "fetch"),
                         stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 1)
  })

  it("accepts fetch subcommand with --gse-id", {
    # Stub implementation — just validates it doesn't crash
    exit_code <- system2("Rscript",
      c("scripts/main.R", "fetch", "--gse-id", "GSE100155"),
      stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 0)
  })

  it("accepts qc subcommand with --input", {
    exit_code <- system2("Rscript",
      c("scripts/main.R", "qc", "--input", "test.csv"),
      stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 0)
  })

  it("accepts clean subcommand with --input", {
    exit_code <- system2("Rscript",
      c("scripts/main.R", "clean", "--input", "test.csv"),
      stdout = FALSE, stderr = FALSE)
    expect_equal(exit_code, 0)
  })
})

describe("NDJSON output", {

  it("produces valid JSON lines on stdout", {
    result <- system2("Rscript",
      c("scripts/main.R", "fetch", "--gse-id", "GSE100155"),
      stdout = TRUE, stderr = FALSE)

    # Each line should be parseable JSON
    for (line in result) {
      parsed <- tryCatch(jsonlite::fromJSON(line),
                         error = function(e) NULL)
      expect_false(is.null(parsed),
        sprintf("Line is not valid JSON: '%s'", substr(line, 1, 80)))
    }
  })

  it("has 'level' field in every NDJSON line", {
    result <- system2("Rscript",
      c("scripts/main.R", "fetch", "--gse-id", "GSE100155"),
      stdout = TRUE, stderr = FALSE)

    for (line in result) {
      parsed <- jsonlite::fromJSON(line)
      expect_true("level" %in% names(parsed),
        sprintf("NDJSON line missing 'level' field: '%s'", substr(line, 1, 80)))
    }
  })
})
```

- [ ] **Step 2: Run test, verify FAILS (main.R doesn't exist yet)**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-main.R")'
```

Expected: FAIL — "cannot open the connection" or "No such file or directory"

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-main.R
git commit -m "test: add CLI dispatch and NDJSON integration tests

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 6: Create scripts/main.R and scripts/report.R

**Files:**
- Create: `scripts/main.R`
- Create: `scripts/report.R`

- [ ] **Step 1: Write scripts/report.R (NDJSON helpers)**

Create directory and file:

```bash
mkdir -p scripts
```

Create `scripts/report.R`:

```r
# report.R — NDJSON reporting helpers for geo-microarray-processing
#
# All output to stdout is valid NDJSON. Each function writes one line.
# Use report_info() for progress, report_result() for final output,
# and report_error() for terminal failures.

#' Write an info-level NDJSON message to stdout
#'
#' @param msg Character string with progress message
#' @param ... Additional named fields to include in the JSON object
report_info <- function(msg, ...) {
  extra <- list(...)
  obj <- c(list(level = "info", msg = msg), extra)
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
}

#' Write a result-level NDJSON message to stdout
#'
#' @param status Character status code (success_matrix, error, etc.)
#' @param files List of file info (each with path, rows, cols)
#' @param metadata List of metadata (platform, organism, n_samples)
#' @param ... Additional named fields
report_result <- function(status, files = list(), metadata = list(), ...) {
  extra <- list(...)
  obj <- c(list(level = "result", status = status,
                files = files, metadata = metadata), extra)
  # Remove empty lists for cleaner output
  if (length(files) == 0) obj$files <- NULL
  if (length(metadata) == 0) obj$metadata <- NULL
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
}

#' Write an error-level NDJSON message and exit
#'
#' @param msg Error message
#' @param exit_code Integer exit code (default: 1)
report_error <- function(msg, exit_code = 1) {
  obj <- list(level = "error", msg = msg)
  cat(jsonlite::toJSON(obj, auto_unbox = TRUE), "\n", sep = "")
  quit(status = exit_code)
}
```

- [ ] **Step 2: Write scripts/main.R**

Create `scripts/main.R`:

```r
#!/usr/bin/env Rscript
#
# main.R — Single entry point for geo-microarray-processing node
#
# Usage:
#   Rscript scripts/main.R fetch --gse-id GSE100155 --outdir ./output
#   Rscript scripts/main.R qc    --input ./output/expr_gene.csv
#   Rscript scripts/main.R clean --input ./output/expr_gene.csv
#
# The first positional argument is the subcommand.
# All parameters declared in SKILL.md frontmatter are accepted.
# Output is NDJSON to stdout.

source("scripts/report.R")

# -------------------------------------------------------------------
# Argument parsing
# -------------------------------------------------------------------

#' Parse command-line arguments
#'
#' @param args Character vector of CLI args (default: commandArgs)
#' @return Named list of parsed values
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {

  if (length(args) == 0) {
    cat("Usage: Rscript scripts/main.R <subcommand> [options]\n")
    cat("\nSubcommands:\n")
    cat("  fetch   Download and process GEO microarray data\n")
    cat("  qc      Quality check on expression matrix\n")
    cat("  clean   Clean and normalize expression matrix\n")
    cat("\nOptions:\n")
    cat("  --gse-id ID     GEO Series identifier (fetch)\n")
    cat("  --outdir DIR    Output directory (default: .)\n")
    cat("  --input FILE    Input expression matrix (qc, clean)\n")
    cat("  --proxy URL     HTTP/HTTPS proxy for GEO access\n")
    cat("  --api-key KEY   NCBI API key for higher rate limits\n")
    quit(status = 1)
  }

  subcommand <- args[1]

  valid_subcommands <- c("fetch", "qc", "clean")
  if (!subcommand %in% valid_subcommands) {
    report_error(sprintf("Unknown subcommand '%s'. Valid: %s",
      subcommand, paste(valid_subcommands, collapse = ", ")))
  }

  # Parse --key=value and --key value pairs
  opts <- list(
    subcommand      = subcommand,
    gse_id          = NULL,
    outdir          = ".",
    input           = NULL,
    proxy           = NULL,
    api_key         = Sys.getenv("NCBI_API_KEY", unset = NULL)
  )

  remaining <- args[-1]
  i <- 1
  while (i <= length(remaining)) {
    key <- remaining[i]

    if (key == "--gse-id") {
      i <- i + 1
      if (i > length(remaining)) {
        report_error("--gse-id requires a value")
      }
      opts$gse_id <- remaining[i]
    } else if (key == "--outdir") {
      i <- i + 1
      if (i <= length(remaining)) opts$outdir <- remaining[i]
    } else if (key == "--input") {
      i <- i + 1
      if (i <= length(remaining)) opts$input <- remaining[i]
    } else if (key == "--proxy") {
      i <- i + 1
      if (i <= length(remaining)) opts$proxy <- remaining[i]
    } else if (key == "--api-key") {
      i <- i + 1
      if (i <= length(remaining)) opts$api_key <- remaining[i]
    } else if (startsWith(key, "--gse-id=")) {
      opts$gse_id <- sub("^--gse-id=", "", key)
    } else if (startsWith(key, "--outdir=")) {
      opts$outdir <- sub("^--outdir=", "", key)
    } else if (startsWith(key, "--input=")) {
      opts$input <- sub("^--input=", "", key)
    } else if (startsWith(key, "--proxy=")) {
      opts$proxy <- sub("^--proxy=", "", key)
    } else if (startsWith(key, "--api-key=")) {
      opts$api_key <- sub("^--api-key=", "", key)
    } else {
      report_error(sprintf("Unknown option: %s", key))
    }
    i <- i + 1
  }

  # Validate required args per subcommand
  if (opts$subcommand == "fetch" && is.null(opts$gse_id)) {
    report_error("fetch subcommand requires --gse-id")
  }
  if (opts$subcommand %in% c("qc", "clean") && is.null(opts$input)) {
    report_error(sprintf("%s subcommand requires --input", opts$subcommand))
  }

  return(opts)
}

# -------------------------------------------------------------------
# Subcommand stubs (filled in by port-fetch-geo, add-qc, add-clean)
# -------------------------------------------------------------------

#' Fetch subcommand — stub
do_fetch <- function(opts) {
  report_info(sprintf("Fetching GEO data for %s...", opts$gse_id))
  report_info("Subcommand 'fetch' not yet implemented")
  report_result("success_matrix", files = list(), metadata = list())
}

#' QC subcommand — stub
do_qc <- function(opts) {
  report_info(sprintf("Running QC on %s...", opts$input))
  report_info("Subcommand 'qc' not yet implemented")
  report_result("pass", files = list(), metadata = list())
}

#' Clean subcommand — stub
do_clean <- function(opts) {
  report_info(sprintf("Cleaning %s...", opts$input))
  report_info("Subcommand 'clean' not yet implemented")
  report_result("success", files = list(), metadata = list())
}

# -------------------------------------------------------------------
# Main dispatch
# -------------------------------------------------------------------

main <- function() {
  opts <- parse_args()

  switch(opts$subcommand,
    fetch = do_fetch(opts),
    qc    = do_qc(opts),
    clean = do_clean(opts)
  )
}

if (sys.nframe() == 0) {
  main()
}
```

- [ ] **Step 3: Run test-main.R, verify all tests PASS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-main.R")'
```

Expected: all 8 tests PASS (6 CLI dispatch + 2 NDJSON format)

- [ ] **Step 4: Commit**

```bash
git add scripts/main.R scripts/report.R
git commit -m "feat: add main.R entry point with subcommand dispatch and NDJSON reporting

main.R: CLI arg parsing, subcommand dispatch (fetch/qc/clean),
NCBI_API_KEY env var fallback, stub handlers for all subcommands.
report.R: report_info(), report_result(), report_error() NDJSON helpers.
Tested: 8 integration tests pass for CLI dispatch and NDJSON output.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 7: Create tests/fixtures/ directory placeholder

**Files:**
- Create: `tests/fixtures/.gitkeep`

- [ ] **Step 1: Create fixtures directory with placeholder**

```bash
mkdir -p tests/fixtures
touch tests/fixtures/.gitkeep
```

- [ ] **Step 2: Add fixture file listing to CLAUDE.md**

Append this note to CLAUDE.md after the "R Development" section:

```bash
# Fixture files (tests/fixtures/) are not tracked in git.
# Generate them from GEO test datasets before running module-level tests.
# See tests/fixtures/README.md for generation instructions.
```

Actually, let me create a `tests/fixtures/README.md` instead:

Create `tests/fixtures/README.md`:

```markdown
# Test Fixtures for geo-microarray-processing

Test fixtures are pre-saved RDS files used by module-level tests.
They eliminate live GEO network calls during testing.

## Required Fixtures

| File | Source | How to Generate |
|---|---|---|
| `GSE100155_eset.rds` | GSE100155 | `getGEO("GSE100155", GSEMatrix=TRUE)[[1]]` then `saveRDS()` |
| `GSE12345_eset_list.rds` | Multi-platform GSE | `getGEO("GSE12345", GSEMatrix=TRUE)` then `saveRDS()` |
| `GSE_methylation_meta.rds` | Methylation GSE | GEO metadata with BPM+IDAT detected |
| `cel_valid.rds` | Any valid CEL | `oligo::read.celfiles("valid.CEL")` then `saveRDS()` |
| `cel_corrupted.rds` | Corrupted CEL | Empty/truncated file fixture |
| `gpl570_annotation.rds` | GPL570 | `getGEO("GPL570")` then `Table()` then `saveRDS()` |
| `expr_raw.rds` | Raw matrix | Matrix with 99pct > 100 |
| `expr_centered.rds` | Centered matrix | Matrix with mean ≈ 0, 99pct < 10 |
| `expr_log.rds` | Log matrix | Log2-transformed matrix |

## Generation Script

```r
# Run once to generate fixtures — requires active conda env with GEOquery
library(GEOquery)
library(Biobase)

# Single-platform ExpressionSet
gse <- getGEO("GSE100155", GSEMatrix = TRUE)[[1]]
saveRDS(gse, "tests/fixtures/GSE100155_eset.rds")

# GPL570 annotation
gpl <- getGEO("GPL570")
gpl_table <- Table(gpl)
saveRDS(gpl_table, "tests/fixtures/gpl570_annotation.rds")

# Raw expression (99pct > 100)
expr_raw <- matrix(runif(100 * 5, 0, 500), nrow = 100, ncol = 5)
rownames(expr_raw) <- paste0("probe_", 1:100)
colnames(expr_raw) <- paste0("sample_", 1:5)
saveRDS(expr_raw, "tests/fixtures/expr_raw.rds")

# Centered expression (mean ~0, range narrow)
set.seed(42)
expr_centered <- matrix(rnorm(100 * 5, 0, 2), nrow = 100, ncol = 5)
rownames(expr_centered) <- paste0("probe_", 1:100)
colnames(expr_centered) <- paste0("sample_", 1:5)
saveRDS(expr_centered, "tests/fixtures/expr_centered.rds")

# Log expression (typical microarray range)
expr_log <- matrix(runif(100 * 5, 2, 14), nrow = 100, ncol = 5)
rownames(expr_log) <- paste0("probe_", 1:100)
colnames(expr_log) <- paste0("sample_", 1:5)
saveRDS(expr_log, "tests/fixtures/expr_log.rds")

message("Fixtures generated in tests/fixtures/")
```
```

- [ ] **Step 3: Commit**

```bash
git add tests/fixtures/.gitkeep tests/fixtures/README.md
git commit -m "feat: add tests/fixtures/ directory and generation guide

Create placeholder for fixture RDS files (generated on-demand, not tracked).
README.md documents all required fixtures with generation scripts.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 8: Final verification — run all tests

- [ ] **Step 1: Run the full test suite**

```bash
Rscript -e 'testthat::test_dir("tests/testthat/")'
```

Expected: all tests PASS (6 test-skillmd.R + 8 test-main.R = 14 tests)

- [ ] **Step 2: Verify SKILL.md frontmatter is valid YAML**

```bash
Rscript -e '
  content <- readLines("SKILL.md", warn = FALSE)
  dashes <- which(content == "---")
  yaml_lines <- content[(dashes[1] + 1):(dashes[2] - 1)]
  tmp <- tempfile(fileext = ".yaml")
  writeLines(yaml_lines, tmp)
  fm <- yaml::read_yaml(tmp)
  cat("name:", fm$name, "\n")
  cat("type:", fm$type, "\n")
  cat("parameters:", length(fm$parameters), "\n")
  cat("exceptions:", length(fm$exceptions), "\n")
  cat("OK\n")
'
```

Expected: prints `name: geo-microarray-processing`, `type: standard`, `parameters: 6`, `exceptions: 6`, `OK`

- [ ] **Step 3: Verify env YAML files are valid**

```bash
python3 -c "
import yaml
for f in ['env.yaml', 'env-4.5.yaml']:
    d = yaml.safe_load(open(f))
    print(f'{f}: name={d[\"name\"]}, channels={d[\"channels\"]}, deps={len(d[\"dependencies\"])}')
print('OK')
"
```

Expected: prints summary for both env files, `OK`

- [ ] **Step 4: Commit final state**

```bash
git add -A
git commit -m "test: final verification — all 14 tests pass, SKILL.md and env files valid

Test results:
- test-skillmd.R: 6/6 pass
- test-main.R: 8/8 pass
- SKILL.md frontmatter: valid YAML, all required fields present
- env.yaml + env-4.5.yaml: valid YAML, correct channel/dep structure

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```
