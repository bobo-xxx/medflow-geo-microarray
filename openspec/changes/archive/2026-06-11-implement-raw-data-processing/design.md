## Architecture

### Router logic

```
fetch_geo_data(gse_id, ...)
  │
  ├─ Tier 2: series matrix? → pipeline detection from metadata → done ✅
  │
  ├─ Tier 3: supplementary matrix? (*.txt.gz in suppl/) → pipeline detection → done
  │
  ├─ Tier 4: raw files in RAW.tar?
  │   │
  │   ├─ download + extract RAW.tar
  │   ├─ detect_raw_type(files)        ← new
  │   ├─ route to processor             ← new
  │   └─ shared downstream path         ← existing (validate, annotate, aggregate, save)
  │
  └─ Tier 5: metadata only
```

### `detect_raw_type(files)` — platform detection from filenames

```r
detect_raw_type <- function(files) {
  if (any(grepl("[.]CEL([.]gz)?$", files, ignore.case = TRUE))) return("affymetrix")
  if (any(grepl("[.]BPM$", files, ignore.case = TRUE))) return("methylation")
  if (any(grepl("[.]idat([.]gz)?$", files, ignore.case = TRUE))) return("illumina")
  if (any(grepl("[.]GPR([.]gz)?$", files, ignore.case = TRUE))) return("agilent_2c")
  if (any(grepl("[.]PAIR([.]gz)?$", files, ignore.case = TRUE))) return("nimblegen")
  txt_files <- grep("[.]txt$", files, value = TRUE, ignore.case = TRUE)
  if (length(txt_files) > 0) {
    header <- readLines(txt_files[1], n = 5)
    if (any(grepl("ProbeName|GeneName|gTotalGeneSignal|gProcessedSignal", header)))
      return("agilent_1c")
  }
  "unknown"
}
```

### Platform processors

| Type | Processor | Pipeline |
|---|---|---|
| affymetrix | `process_affy(files, out_dir)` | `oligo::read.celfiles()` → `rma()` |
| illumina | `process_illumina(files, out_dir)` | `limma::read.idat()` → `neqc()` |
| agilent_2c | `process_agilent_2c(files, out_dir)` | `read.maimages()` → header detect → normexp → loess → quantile |
| agilent_1c | `process_agilent_1c(files, out_dir)` | `read.maimages(source="agilent", green.only=TRUE)` → quantile |
| nimblegen | `process_nimblegen(files, out_dir)` | `oligo::read.xys()` → `rma()` |
| methylation | SKIP | warning + skip |
| unknown | SKIP | warning + fallback to Tier 5 |

### Shared downstream path (after any processor)

```
log2 probe-level matrix
  → validate_expr_matrix()
  → normalize_expr_matrix()    (pass-through, already log2)
  → 5-tier gene annotation
  → aggregate to gene-level
  → post-normalization validation
  → save CSV + metadata + NDJSON
```

### Files

| File | Change |
|---|---|
| `node/scripts/raw.R` | NEW: `detect_raw_type()`, `process_raw_files()`, 5 processors |
| `node/scripts/fetch.R` | MODIFIED: wire Tiers 3-5 |
| `node/scripts/main.R` | MODIFIED: source raw.R |
| `tests/testthat/test-raw.R` | NEW: platform detection + processor tests |
