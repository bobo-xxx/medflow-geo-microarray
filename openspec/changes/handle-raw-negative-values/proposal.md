## Why

GEO returns expression data with inconsistent preprocessing — RMA log2, GCOS-scaled linear, GenomeStudio background-subtracted, TAC/SST-RMA log2, or raw CEL/IDAT with no preprocessing. The current `detect_expr_type()` heuristic (99pct > 100 → raw) cannot distinguish MAS5 from quantile-normalized linear data, and the `log2(x+1)` transform produces NaN for negative values common in Illumina data. GEO's `data_processing` metadata contains the submitter's preprocessing description, which we now know is reliable (5/5 datasets tested showed perfect consistency between claim and data). This change replaces heuristic detection with metadata-driven pipeline identification and platform-appropriate transforms.

## What Changes

- **New**: `detect_pipeline()` — parse `data_processing` column(s) for preprocessing keywords (RMA, GCOS, MAS5, GenomeStudio, TAC, neqc, lumi, quantile normalization, etc.)
- **New**: Pipeline-specific normalization transforms:
  - Affy RMA/GCRMA/neqc/lumi → pass-through (already log2 + normalized)
  - Affy MAS5/GCOS → `log2(x + 1e-6)` (linear, non-negative)
  - Illumina GenomeStudio avg norm → shift + `log2(x + 1e-6)` (linear, negatives possible)
  - SST-RMA/TAC → pass-through (log2, check gene_assignment column)
  - Agilent FE → pass-through (already normalized)
- **New**: QN detection via CV(colMedians) < 0.002 (scale-invariant)
- **New**: metadata-pipeline consistency check — warn on mismatch
- **Modified**: `normalize_expr_matrix()` replaced by pipeline-specific transforms
- **Modified**: `post-normalization validation` in fetch.R (already added)
- **Removed**: heuristic `detect_expr_type()` for transformation decisions (kept for QC reporting only)

## Capabilities

### New Capabilities
- `pipeline-detection`: Parse GEO `data_processing` text to identify the preprocessing pipeline applied by the submitter
- `pipeline-normalization`: Apply platform-appropriate transform based on detected pipeline (log2, shift+log2, or pass-through)

### Modified Capabilities
- `normalization`: `detect_expr_type()` demoted from decision-making to reporting-only; transforms now driven by pipeline detection

## Impact

- `node/scripts/normalize.R`: replaced by `node/scripts/pipeline.R` with `detect_pipeline()`, `apply_pipeline_transform()`
- `node/scripts/fetch.R`: call pipeline detection after extracting pData, apply pipeline-appropriate transform
- `node/scripts/qc.R`: report QN status from CV(colMedians)
- `node/scripts/species.R`: no change
- `node/scripts/annotate.R`: no change
- `node/SKILL.md`: update normalize description to reference pipeline detection
- `tests/testthat/test-normalize.R`: update to test pipeline detection
