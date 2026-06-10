# geo-microarray-processing — Design Doc

**Date:** 2026-06-10  
**Status:** Approved  
**Pipeline:** superpowers brainstorming → writing-plans → OpenSpec (/opsx:apply with TDD)

## Overview

This node fetches and processes GEO microarray expression data, adapted from the reference implementation in `original/geo-microarray-fetch.zip` to the IRE node-package v2 format. The original is a 779-line R script that already covers fetch, normalization, probe-to-gene aggregation, and validation. The port is primarily structural: single-entry `main.R` with subcommand dispatch, NDJSON reporting, and a declarative `env.yaml`.

## Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Scope | Fetch + QC + clean | Matches original scope; QC/clean added as incremental changes |
| Language | R + Bioconductor | Deep dependency on `oligo::rma()`, `limma::read.maimages()`, `GEOquery`, `Biobase::ExpressionSet` — no Python equivalents |
| Architecture | Functional modules (Option B) | Minimal abstraction, close to original code, fits Bioconductor idioms. `main.R` dispatcher + `fetch.R`, `normalize.R`, `annotate.R`, `validate.R`, `report.R` |
| Subcommands | `fetch`, `qc`, `clean` (Option A) | Each independently callable; orchestrator wires them. Agents get per-stage decision points for error handling |
| Multi-platform output | Per-platform files (Option A) | Established pattern (both `virtualArray` and `crossmeta` expect per-platform input). Merging deferred to downstream node or future `add-platform-merge` change |
| R versions | 4.3 / BioC 3.18 (prod) + 4.5 / BioC 3.20 (forward-compat) | Bioconductor is tightly coupled to R version; test both |
| Conda channels | Upstream only in `env.yaml` | Mirrors configured locally via `conda config`, never committed |
| Testing | Hybrid (Option C) | Pure functions → unit tests; stateful modules → module-level with fixture RDS; `main.R` dispatch → integration tests |
| Dev approach | TDD per OpenSpec change | Sequential pipeline-shaped work; CLAUDE.md mandates TDD |
| API key | `NCBI_API_KEY` env var, fallback to `--api-key` config | No hardcoded secrets; runtime resolution |
| Proxy | `--proxy` with `bind: config` | Original hardcoded `localhost:1086` removed |

## Architecture

### Package Structure

```
geo-microarray-processing@1.0.0/
├── SKILL.md                   # Agent contract (frontmatter) + narrative
├── env.yaml                   # Conda env: R 4.3 + Bioconductor 3.18
├── env-4.5.yaml               # Forward-compat: R 4.5 + Bioconductor 3.20
├── scripts/
│   ├── main.R                 # Single entry point, subcommand dispatch
│   ├── fetch.R                # GEO download, 5-tier fallback, platform detection
│   ├── normalize.R            # detect_expr_type, normalize_expr_matrix, RMA
│   ├── annotate.R             # GPL annotation download, probe-to-gene aggregation
│   ├── validate.R             # CEL integrity, matrix QC, gene expression validation
│   └── report.R               # NDJSON output helpers
├── tests/
│   ├── testthat/
│   │   ├── test-fetch.R
│   │   ├── test-normalize.R
│   │   ├── test-annotate.R
│   │   ├── test-validate.R
│   │   ├── test-report.R
│   │   ├── test-main.R
│   │   └── helpers.R
│   └── fixtures/
│       ├── GSE100155_eset.rds
│       ├── GSE12345_eset_list.rds
│       ├── GSE_methylation_meta.rds
│       ├── cel_valid.rds
│       ├── cel_corrupted.rds
│       ├── gpl570_annotation.rds
│       ├── expr_raw.rds
│       ├── expr_centered.rds
│       └── expr_log.rds
└── references/
    ├── ERROR_CODES.md
    └── PLATFORMS.md
```

### Data Flow (fetch subcommand)

```
GEO Database
    │
    ▼
fetch.R        5-tier fallback (local → series matrix → suppl → raw → meta)
               Platform detection (Affy/Agilent/Illumina), methylation skip
    │ expr_matrix (probes × samples)
    ▼
normalize.R    detect_expr_type() → raw|centered|log
               normalize_expr_matrix() → log2(x+1), shift ≥ 0
    │
    ▼
annotate.R     get_gpl_annotation() → probe-to-gene mapping
               aggregate_probe_to_gene() → gene-level (mean aggregation, /// split)
    │
    ▼
validate.R     validate_expr_matrix(), validate_gene_expression()
    │
    ▼
Output         Per-platform CSV files:
               ├── expr_probe_{gse_id}_{gpl}.csv
               └── expr_gene_{gse_id}_{gpl}.csv
               
               NDJSON to stdout:
               {"level":"info","msg":"..."}
               {"level":"result","status":"...",...}
```

## Subcommands

```
Rscript scripts/main.R fetch --gse-id GSE100155 --outdir ./output
Rscript scripts/main.R qc    --input ./output/GSE100155/expr_gene_GSE100155_GPL570.csv
Rscript scripts/main.R clean --input ./output/GSE100155/expr_gene_GSE100155_GPL570.csv
```

## SKILL.md Frontmatter (Key Fields)

- **inputs:** `[]` — fetch pulls from GEO, not upstream nodes
- **outputs:** Two CSV patterns with `{gse_id}` and `{gpl}` variables
- **parameters:** `subcommand` (upstream), `--gse-id` (upstream), `--outdir` (framework), `--input` (upstream), `--proxy` (config), `--api-key` (config)
- **exceptions:** 6 patterns (E001–E005 errors → `skip_with_warning`, all-failed → `halt`)
- **hardware:** 4 GB / 2 CPU, no GPU

## NDJSON Report Format

```json
{"level":"info","msg":"Downloading suppl files for GSE100155..."}
{"level":"info","msg":"Platform detected: Affymetrix (GPL570)"}
{"level":"result","status":"success_matrix","files":[
  {"path":"output/GSE100155/expr_probe_GSE100155_GPL570.csv","rows":54675,"cols":24},
  {"path":"output/GSE100155/expr_gene_GSE100155_GPL570.csv","rows":20838,"cols":24}
],"metadata":{"platform":"GPL570","organism":"Homo sapiens","n_samples":24}}
```

## Testing Strategy

**Unit tests** for pure functions:
- `detect_expr_type()` — quantile-based classification
- `normalize_expr_matrix()` — transform correctness
- `aggregate_probe_to_gene()` — `///` split, mean aggregation
- `validate_cel_integrity()` — file corruption detection

**Module-level tests** for stateful modules (fixture RDS, no live GEO):
- `do_fetch()` — fallback logic, platform detection, error handling
- `do_qc()` — QC thresholds, outlier detection
- `do_clean()` — normalization pipeline

**Integration tests** for `main.R`:
- CLI arg parsing, subcommand dispatch
- NDJSON output format validation
- `NCBI_API_KEY` env var fallback

**CI:** Matrix across R 4.3 and R 4.5

## OpenSpec Build Order

1. **`init-node-package`** — `SKILL.md` frontmatter, `env.yaml`, `env-4.5.yaml`, `scripts/main.R` skeleton, `references/`, test fixtures
2. **`port-fetch-geo`** — `fetch.R`, `normalize.R`, `annotate.R`, `validate.R`, `report.R` ported from original, NDJSON reporting
3. **`add-qc-subcommand`** — standalone QC subcommand, expression matrix validation
4. **`add-clean-subcommand`** — downstream cleaning/normalization
5. *(future)* **`add-platform-merge`** — optional union merge for cross-platform needs

## Reference

- Original implementation: `original/geo-microarray-fetch.zip`
- Node-package v2 spec: `openspec/specs/node-package/spec.md`
- Core connection protocol: `openspec/specs/core-connection/spec.md`
- GEO FTP README: https://ftp.ncbi.nlm.nih.gov/geo/README.txt
- Bioconductor: https://bioconductor.org/
- `crossmeta` (meta-analysis pattern): https://bioconductor.org/packages/crossmeta
- `virtualArray` (merge pattern): Heider & Alt, BMC Bioinformatics 2013
