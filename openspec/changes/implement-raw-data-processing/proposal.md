## Why

The current fetch pipeline only implements Tier 2 (series matrix). When no series matrix exists, or when metadata-data consistency fails, the node must fall back to raw file processing (Tiers 3-4) or metadata-only (Tier 5). Without this, the 5-tier fallback is incomplete and many GEO datasets are unreachable.

## What Changes

- **Tier 3**: Download supplementary files, check for author-provided processed matrix (*.txt.gz), apply pipeline detection from metadata
- **Tier 4**: Download RAW.tar, extract raw files, detect platform from file patterns, route to platform-specific processor, produce log2-scale expression matrix
- **Tier 5**: Return metadata only when all data retrieval fails
- **New**: `node/scripts/raw.R` — `detect_raw_type()`, `process_raw_files()`, plus platform-specific processors: `process_cel()`, `process_idat()`, `process_gpr()`, `process_agilent_fe()`, `process_pair()`
- **Modified**: `node/scripts/fetch.R` — wire Tiers 3-5 into `fetch_geo_data()`

## Capabilities

### New Capabilities
- `raw-data-processing`: Detect raw file type from file patterns, route to platform-specific Bioconductor pipeline, produce log2-scale probe-level matrix
- `platform-processors`: Six processors (CEL→RMA, IDAT→neqc, GPR→limma+loess, Agilent FE→limma, PAIR→RMA, Methylation→skip)

## Impact

- `node/scripts/raw.R` — NEW
- `node/scripts/fetch.R` — MODIFIED (Tiers 3-5)
- `node/references/PLATFORMS.md` — already updated ✅
- `tests/testthat/test-raw.R` — NEW
