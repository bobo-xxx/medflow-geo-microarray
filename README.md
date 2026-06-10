# medflow-geo-microarray

GEO microarray expression data processing node for the IRE agentic bioinformatics workflow framework.

## What it does

Fetches, quality-checks, and normalizes microarray expression data from NCBI GEO with a 5-tier fallback strategy.

```
fetch → QC → clean
```

| Subcommand | Purpose |
|---|---|
| `fetch` | Download GEO data via 5-tier fallback (cache → series matrix → supplementary → raw CEL → metadata-only) |
| `qc` | Quality check with 7 metrics, returns pass / caution / rerun / veto |
| `clean` | Normalize expression values to log2 scale |

All output is NDJSON to stdout — machine-readable by downstream agents and nodes.

## Quick start

```bash
# Create environment
conda env create -f node/env.yaml
conda activate geo-microarray-processing

# Install AnnoProbe (CRAN-only, not on conda)
Rscript -e 'install.packages("AnnoProbe", repos = "https://cloud.r-project.org")'

# Fetch GEO data
Rscript node/scripts/main.R fetch --gse-id GSE318047 --outdir ./output

# Quality check
Rscript node/scripts/main.R qc --input ./output/GSE318047/expr_gene_GSE318047.csv

# Normalize
Rscript node/scripts/main.R clean --input ./output/GSE318047/expr_probe_GSE318047.csv
```

## Supported platforms

| Platform | File Pattern | Processing |
|---|---|---|
| Affymetrix | `*.CEL(.gz)?` | `oligo::rma()` |
| Agilent | `*.GPR(.gz)?` | `limma::read.maimages()` |
| Agilent FE | TXT with ProbeName/GeneName | `limma::read.maimages(source="agilent")` |
| Illumina | `*.idat` (no BPM) | `illuminaio::readIDAT()` |
| TXT/CSV | `series_matrix.txt.gz` | Direct import |
| Methylation | `*.idat` + `*.BPM` | **SKIP** |

## Species

Human (taxId 9606), mouse (10090), rat (10116). Other species pass through with GPL table + AnnoProbe annotation.

## Gene annotation

5-tier fallback: `fData()` direct column → `gene_assignment` parsing → GPL annotation table → AnnoProbe pipe → probe IDs.

## QC metrics

| Metric | Warning | Hard | Hard Decision |
|---|---|---|---|
| Missing values | > 10% | > 50% | veto |
| Zero-variance genes | > 20% | — | caution |
| Extreme values | > 1e100 | — | caution |
| Sample count | < 6 | — | caution |
| Gene count | < 5000 | < 100 | rerun |
| Sample correlation | < 0.8 | < 0.5 | rerun |
| Expression scale | not log2 | — | caution |

## Development

```bash
# R 4.3 (production)
conda env create -f node/env.yaml
conda activate geo-microarray-processing

# R 4.5 (forward-compat)
conda env create -f node/env-4.5.yaml
conda activate geo-microarray-processing-4.5

# Run tests
Rscript -e 'testthat::test_dir("tests/testthat/")'
```

133 tests, 0 failures across both R 4.3 and R 4.5.

## Architecture

```
node/
├── SKILL.md              # Agent contract
├── env.yaml              # Conda environment
├── scripts/
│   ├── main.R            # Entry point + subcommand dispatch
│   ├── report.R          # NDJSON helpers
│   ├── fetch.R           # 5-tier fallback + platform detection
│   ├── qc.R              # 7 metrics + decision engine
│   ├── clean.R           # log2 normalization
│   ├── normalize.R       # detect_expr_type + normalize_expr_matrix
│   ├── annotate.R        # 5-tier gene annotation
│   ├── validate.R        # Expression matrix + CEL validation
│   └── species.R         # taxId → species name
└── references/
    ├── ERROR_CODES.md
    └── PLATFORMS.md
```

## License

Research use. Data fetched from GEO is subject to GEO's terms of use.
