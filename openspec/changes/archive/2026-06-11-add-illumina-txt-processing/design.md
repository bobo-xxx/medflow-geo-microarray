## Architecture

`process_illumina_txt()` handles Illumina GenomeStudio non-normalized TXT exports — the format used by ~80% of Illumina GEO submissions.

### Detection

`detect_raw_type()` checks TXT headers for `TargetID` and `AVG_Signal` columns (after Agilent FE check, before falling to `unknown`).

### Processing

1. `read.table()` — parse tab-delimited TXT with `ID_REF` row names
2. Split interleaved columns: odd = expression values, even = detection p-values
3. Build `EListRaw` with `$E` (expression) and `$other$Detection` (p-values)
4. `limma::neqc()` — normexp bg → offset +16 → QN → log2
5. Extract `$E` matrix

### Files

| File | Change |
|---|---|
| `node/scripts/raw.R` | NEW function + modified detect/routing |
| `tests/testthat/test-raw.R` | NEW tests |
