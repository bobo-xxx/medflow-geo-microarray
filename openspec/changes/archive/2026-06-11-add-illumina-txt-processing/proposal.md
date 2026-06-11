## Why

~80% of Illumina expression submissions on GEO provide non-normalized TXT exports from GenomeStudio, not binary IDAT files. Our current `process_illumina()` only handles IDAT, leaving most Illumina data unprocessable through the raw file path. These TXT files contain `AVG_Signal` + `Detection Pval` columns per probe and can be processed by `limma::read.ilmn()` → `neqc()` — the standard Bioconductor pipeline documented in the limma User's Guide.

## What Changes

- **New**: `process_illumina_txt()` — read non-normalized TXT via `limma::read.ilmn()`, process via `neqc()`
- **Modified**: `detect_raw_type()` — recognize Illumina non-normalized TXT files (have `TargetID` and `AVG_Signal` columns)
- **Modified**: `process_raw_files()` routing — dispatch to `process_illumina_txt()` for Illumina TXT

## Capabilities

### New Capabilities
- `illumina-txt-processing`: Process Illumina GenomeStudio non-normalized TXT exports via limma::read.ilmn() → neqc()

### Modified Capabilities  
- None

## Impact

- `node/scripts/raw.R` — NEW function, modified detect + router
- `tests/testthat/test-raw.R` — NEW tests
