# Platform Detection & Raw Data Processing

This document defines detection rules for microarray platform file types
and the standard preprocessing pipeline for each raw data type.

## File Pattern Detection

| File Pattern | Platform | Raw Processing | Output Scale |
|---|---|---|---|
| `*.CEL(.gz)?` | Affymetrix | `oligo::rma()` (bg correct → QN → median polish) | Log2 |
| `*.idat` + `*.BPM` | Methylation | **SKIP** — not expression data | — |
| `*.idat(.gz)?` only | Illumina BeadChip | `limma::neqc()` (normexp bg → offset 16 → QN) | Log2 |
| `*.GPR(.gz)?` | Agilent 2-color | `limma::read.maimages()` → `backgroundCorrect(normexp)` → `normalizeWithinArrays(loess)` → `normalizeBetweenArrays(quantile)` | Log2 |
| `*.txt` (Agilent FE) | Agilent single-color | `limma::read.maimages(source="agilent", green.only=TRUE)` → `normalizeBetweenArrays(quantile)` | Log2 |
| `*.PAIR(.gz)?` | NimbleGen | `oligo::read.xys()` → `rma()` | Log2 |
| `series_matrix.txt.gz` | TXT/CSV | Pipeline detection from `data_processing` metadata | Varies |

## Raw Data Processing Pipelines

### Affymetrix CEL

```
read.celfiles(files)          # oligo
    │
    ▼
rma()                        # oligo
  ├─ Background correct       (convolution model, signal+noise)
  ├─ Quantile normalize       (force identical distributions)
  ├─ Median polish summarization (probe set → expression value)
  └─ Log2 transform
    │
    ▼
Probe-level log2 expression matrix
```

Standard Bioconductor: `oligo::rma(celfiles)` or `affy::rma(celfiles)` (3' IVT arrays).

### Illumina BeadChip IDAT

```
read.idat(files, bgxfile)     # limma
    │
    ▼
neqc()                       # limma
  ├─ Normexp background       (negative control probes)
  ├─ Offset +16               (stabilize low intensities)
  ├─ Quantile normalize       (force identical distributions)
  └─ Log2 transform
    │
    ▼
Probe-level log2 expression matrix
```

Standard Bioconductor: `limma::neqc(limma::read.idat(files, bgxfile))`.

### Agilent Two-Color GPR

```
read.maimages(files, source="genepix")  # limma
    │
    ▼
backgroundCorrect(method="normexp")     # limma
    │
    ▼
normalizeWithinArrays(method="loess")   # limma (dye bias)
    │
    ▼
normalizeBetweenArrays(method="quantile") # limma
    │
    ▼
Probe-level log2 expression matrix (M-values or A-values)
```

### Agilent Single-Color FE TXT

```
read.maimages(files, source="agilent", green.only=TRUE)  # limma
    │
    ▼
normalizeBetweenArrays(method="quantile")                # limma
    │
    ▼
Probe-level log2 expression matrix
```

Detection: TXT files in RAW.tar with header containing `ProbeName`, `GeneName`,
`gTotalGeneSignal`, or `gProcessedSignal`.

### NimbleGen PAIR

```
read.xys(files)              # oligo
    │
    ▼
rma()                        # oligo
    │
    ▼
Probe-level log2 expression matrix
```

## GPR Header Detection

GPR files require header analysis to distinguish variants:

- `^ATF` → GenePix format (use `source="genepix"`)
- `^TYPE` → Agilent format (use `source="agilent"`)

## Agilent FE TXT Column Detection

TXT files in RAW.tar require column header analysis:

| Columns Present | Channel Type | limma call |
|---|---|---|
| `gProcessedSignal`, `gMeanSignal` | Single-channel (green only) | `read.maimages(source="agilent", green.only=TRUE)` |
| `rProcessedSignal`, `gProcessedSignal` | Two-channel (red/green) | `read.maimages(source="agilent")` |

## Platform → Router Mapping

The Tier 4 raw data router inspects extracted files and dispatches:

```
RAW.tar extracted
  │
  ├─ filelist.txt parsed → files vector
  │
  ├─ *.CEL(.gz)? found     → Affymetrix  → oliRma::rma()
  ├─ BPM + IDAT found      → Methylation → SKIP (with warning)
  ├─ *.idat only           → Illumina    → limma::neqc()
  ├─ *.GPR(.gz)?           → Agilent 2C  → limma::read.maimages() + loess + QN
  ├─ *.txt (Agilent FE)    → Agilent 1C  → limma::read.maimages(source="agilent")
  ├─ *.PAIR(.gz)?          → NimbleGen   → oligo::read.xys() → rma()
  └─ none                  → Tier 5 (metadata only)
```

## Gene Symbol Column Priorities

When extracting gene identifiers from processed data:

1. `Gene Symbol` / `GENE_SYMBOL` / `Symbol` — direct column
2. `gene_assignment` — parse `"ACC // SYMBOL // desc // ..."` → field 2
3. GPL annotation table — `GEOquery::Table(getGEO(GPL))`
4. AnnoProbe pipe — probe FASTA → genome alignment → GENCODE
5. Probe IDs — last resort (with warning)
