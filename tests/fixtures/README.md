# Test Fixtures for geo-microarray-processing

Test fixtures are pre-saved RDS files used by module-level tests.
They eliminate live GEO network calls during testing.

## Test Datasets

| # | GSE | Species | Platform | Samples | Tests |
|---|---|---|---|---|---|
| 1 | GSE318047 | Human | GPL570 | 12 | Standard annotation (Gene Symbol in fData) |
| 2 | GSE156508 | Human | GPL16686 | 12 | No gene symbol — GB_ACC only, tests AnnoProbe fallback |
| 3 | GSE11381 | Mouse | GPL339 | 12 | Mouse single-platform |
| 4 | GSE4105 | Rat | GPL85 | 6 | Rat single-platform |
| 5 | GSE84422 | Human | GPL96+97+570 | ~51 | Multi-platform (3 GPLs in one series) |
| 6 | GSE42861 | Human | GPL13534 (450k) | ~12 | Methylation skip (BPM+IDAT) |

## Required Fixture Files

| File | Source | How to Generate |
|---|---|---|
| `GSE318047_eset.rds` | GSE318047 | `getGEO("GSE318047", GSEMatrix=TRUE)[[1]]` then `saveRDS()` |
| `GSE156508_eset.rds` | GSE156508 | `getGEO("GSE156508", GSEMatrix=TRUE)[[1]]` then `saveRDS()` |
| `GSE11381_eset.rds` | GSE11381 | `getGEO("GSE11381", GSEMatrix=TRUE)[[1]]` then `saveRDS()` |
| `GSE4105_eset.rds` | GSE4105 | `getGEO("GSE4105", GSEMatrix=TRUE)[[1]]` then `saveRDS()` |
| `GSE84422_eset_list.rds` | GSE84422 | `getGEO("GSE84422", GSEMatrix=TRUE)` then `saveRDS()` |
| `GSE42861_meta.rds` | GSE42861 | GEO metadata with BPM+IDAT detected |
| `GPL570_annotation.rds` | GPL570 | `getGEO("GPL570")` then `Table()` then `saveRDS()` |
| `GPL16686_annotation.rds` | GPL16686 | `getGEO("GPL16686")` then `Table()` then `saveRDS()` |
| `cel_valid.rds` | Any valid CEL | `oligo::read.celfiles("valid.CEL")` then `saveRDS()` |
| `cel_corrupted.rds` | Corrupted CEL | Empty/truncated file fixture |
| `expr_raw.rds` | Generated | Matrix with 99pct > 100 |
| `expr_centered.rds` | Generated | Matrix with mean ≈ 0, 99pct < 10 |
| `expr_log.rds` | Generated | Log2-transformed matrix |

## Generation Script

Run once to generate fixtures — requires the conda env with GEOquery and Biobase:

```r
library(GEOquery)
library(Biobase)

# Human GPL570 — standard annotation
gse318047 <- getGEO("GSE318047", GSEMatrix = TRUE)[[1]]
saveRDS(gse318047, "tests/fixtures/GSE318047_eset.rds")

# Human GPL16686 — GB_ACC only, no gene symbol
gse156508 <- getGEO("GSE156508", GSEMatrix = TRUE)[[1]]
saveRDS(gse156508, "tests/fixtures/GSE156508_eset.rds")

# Mouse GPL339
gse11381 <- getGEO("GSE11381", GSEMatrix = TRUE)[[1]]
saveRDS(gse11381, "tests/fixtures/GSE11381_eset.rds")

# Rat GPL85
gse4105 <- getGEO("GSE4105", GSEMatrix = TRUE)[[1]]
saveRDS(gse4105, "tests/fixtures/GSE4105_eset.rds")

# Multi-platform (3 GPLs)
gse84422 <- getGEO("GSE84422", GSEMatrix = TRUE)
saveRDS(gse84422, "tests/fixtures/GSE84422_eset_list.rds")

# GPL570 annotation
gpl570 <- getGEO("GPL570")
gpl570_table <- Table(gpl570)
saveRDS(gpl570_table, "tests/fixtures/GPL570_annotation.rds")

# GPL16686 annotation (GB_ACC only)
gpl16686 <- getGEO("GPL16686")
gpl16686_table <- Table(gpl16686)
saveRDS(gpl16686_table, "tests/fixtures/GPL16686_annotation.rds")

# Raw expression (99pct > 100)
set.seed(123)
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
set.seed(456)
expr_log <- matrix(runif(100 * 5, 2, 14), nrow = 100, ncol = 5)
rownames(expr_log) <- paste0("probe_", 1:100)
colnames(expr_log) <- paste0("sample_", 1:5)
saveRDS(expr_log, "tests/fixtures/expr_log.rds")

message("Fixtures generated in tests/fixtures/")
```
