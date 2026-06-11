## Architecture

### Pipeline Detection (new)

`node/scripts/pipeline.R` — replaces `detect_expr_type()` as the primary decision engine:

```
fetch receives expr_matrix + pData
  │
  ├─ detect_pipeline(pData)
  │   Parse data_processing text → identify submitter's preprocessing
  │
  ├─ Pipeline matched?
  │   ├─ YES → apply pipeline-specific transform
  │   │
  │   └─ NO → fall back to raw file processing
  │       ├─ Affymetrix CEL  → oligo::rma()          (RMA = bg + QN + log2)
  │       ├─ Illumina IDAT    → limma::neqc()         (normexp + QN + log2)
  │       └─ Agilent GPR      → limma backgroundCorrect → loess → quantile
  │
  ├─ Post-normalization validation
  └─ Report pipeline + QN status in NDJSON metadata
```

### Pipeline → Transform Mapping

| Pipeline | Keywords | Transform |
|---|---|---|
| Affy RMA | `RMA`, `robust multi-array` | pass-through (already log2 + QN) |
| Affy MAS5/GCOS | `GCOS`, `MAS5`, `target intensity` | `log2(x + 1e-6)` |
| Affy GCRMA | `GCRMA`, `gc content` | pass-through (already log2 + QN) |
| Affy SST-RMA/TAC | `TAC`, `SST-RMA`, `signal space` | pass-through (already log2, check gene_assignment) |
| Illumina GenomeStudio avg | `Genome Studio`, `average norm` (without `neqc`/`normexp`) | shift + `log2(x + 1e-6)` |
| Illumina neqc | `neqc`, `normexp`, `offset 16` | pass-through (already log2 + QN) |
| Illumina lumi | `lumi`, `VST`, `variance stabil` | pass-through (VST-transformed) |
| Agilent FE | `Feature Extraction`, `Agilent`, `loess` | pass-through (already normalized) |
| None identified | — | fall back to raw file processing |

### QN Detection

5-percentile check: `cv(colQuantiles) < 0.002` at p25, p50, p75, p90, p95.

After quantile normalization, all sample distributions are identical — every percentile is the same across samples. `cv(colQuantile)` at any percentile must approach zero.

```r
is_quantile_normalized <- function(expr, tol = 0.002) {
  pcts <- c(0.25, 0.50, 0.75, 0.90, 0.95)
  vals <- apply(expr, 2, quantile, probs = pcts, na.rm = TRUE)
  cvs <- apply(vals, 1, function(r) sd(r) / abs(mean(r)))
  all(cvs < tol)
}
```

Validated on 8 datasets across 6 platforms (GPL570, GPL96, GPL339, GPL85, GPL6884, GPL16686, GPL13497, GPL341). After explicit QN, all pass (max CV = 0.00158). Before QN, only GSE19804 passes (3.99e-05) — correctly, because it was already QN'd by the submitter.

### Metadata → Data Consistency Check

After detecting pipeline, verify that data matches claimed preprocessing:

| Claim | Expected data | Check |
|---|---|---|
| RMA/log2 | min > 0, max < 20, 99pct < 100 | `detect_expr_type() == "log"` |
| GCOS/MAS5 | min ≥ 0, max > 100, 99pct > 100 | `detect_expr_type() == "raw"` |
| GenomeStudio avg | min < 0, 99pct > 100 | `min < 0 && 99pct > 100` |
| QN | cv(colQuantiles) < 0.002 | 5-percentile CV check |

Mismatch → WARNING in NDJSON, fall back to raw processing.

### File changes

| File | Change |
|---|---|
| `node/scripts/pipeline.R` | NEW: `detect_pipeline()`, `apply_pipeline_transform()`, `is_quantile_normalized()` |
| `node/scripts/normalize.R` | MODIFIED: `detect_expr_type()` retained for reporting only; `normalize_expr_matrix()` calls pipeline transforms |
| `node/scripts/fetch.R` | MODIFIED: call pipeline detection after pData extraction |
| `node/scripts/qc.R` | MODIFIED: report QN status in QC metrics |
| `node/SKILL.md` | MODIFIED: update normalization description |
| `tests/testthat/test-pipeline.R` | NEW: pipeline detection tests |
| `tests/testthat/test-normalize.R` | MODIFIED: update for pipeline-driven transforms |
