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

  - name: phenotype_{gse_id}_{gpl}.csv
    format: csv
    semantic_type: sample_phenotype_table
    description: >
      Sample-level phenotype metadata (pData from ExpressionSet). Contains
      clinical covariates, group labels, and experimental factors for
      downstream grouping and co-factor selection. One file per platform.

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
- **5-tier gene annotation**: fData direct column → gene_assignment
  parsing → GPL annotation table → AnnoProbe pipe → probe IDs

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
| Other | Any | GPL table + AnnoProbe | 2 — pass-through |

Tier-1 species get validated gene symbol mapping against the species
annotation database. Tier-2 species use the 5-tier annotation fallback
chain with a warning; if all tiers fail, probe IDs are used as gene symbols.

### Gene Annotation Fallback

The node tries five sources in priority order:

```
1. fData() direct column          — Gene Symbol / GENE_SYMBOL / Symbol
2. fData() gene_assignment        — parse "ACC // SYMBOL // desc // ..." → extract field 2
3. GPL annotation table           — GEOquery::Table(getGEO(GPL))
4. AnnoProbe pipe alignment       — probe FASTA → genome → GENCODE (147 platforms)
5. Probe IDs as gene symbols      — last resort, with warning
```

Tier 2 handles platforms like GPL17586 (HTA-2.0) where annotation
uses `gene_assignment` format instead of a direct gene symbol column.
Tier 4 handles platforms like GPL16686 (HuGene-2_0-st) where the GPL
table only provides GB_ACC accessions. See
[AnnoProbe on CRAN](https://cran.r-project.org/web/packages/AnnoProbe/)
for supported platforms.

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
    ├── expr_probe_GSE12345_GPL570.csv    # Probe-level (probes × samples)
    ├── expr_gene_GSE12345_GPL570.csv     # Gene-level (genes × samples)
    ├── phenotype_GSE12345_GPL570.csv     # Sample metadata (pData)
    └── suppl/                             # Downloaded supplementary files
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
