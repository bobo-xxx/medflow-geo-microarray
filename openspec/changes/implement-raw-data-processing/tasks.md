### Task 1: Create node/scripts/raw.R — detection and routing

- [x] `detect_raw_type()` — 7 file patterns + Agilent FE header detection
- [x] `detect_gpr_source()` — ATF/TYPE header detection
- [x] `process_raw_files()` — router dispatching to platform processors
- [x] `tests/testthat/test-raw.R` — 16 tests

### Task 2: Create platform processors

- [x] `process_affy()` — `oligo::read.celfiles()` → `rma()`
- [x] `process_illumina()` — `limma::read.idat()` → `neqc()`
- [x] `process_agilent_2c()` — `read.maimages()` → normexp → loess → quantile
- [x] `process_agilent_1c()` — `read.maimages(source="agilent", green.only=TRUE)` → quantile
- [x] `process_nimblegen()` — `oligo::read.xys()` → `rma()`
- [x] Methylation skip — BPM+IDAT detection before processing

### Task 3: Wire Tiers 3-5 into fetch.R

- [x] Refactor Tier 2 per-platform loop into `process_expression_set()`
- [x] Add `process_raw_matrix()` for raw-file-derived matrices
- [x] Tier 3: download suppl, scan for processed matrix
- [x] Tier 4: download RAW.tar, extract, `process_raw_files()`, shared downstream
- [x] Tier 5: `getGEO(GSEMatrix=FALSE)` → metadata only
- [x] Source raw.R in main.R

### Task 4: Tests and verification

- [x] 174 tests pass (0 FAIL) across both R 4.3 and R 4.5
- [x] Tier 2 end-to-end verified on GSE318047
- [ ] Real data verification: test processors against actual raw files (CEL/IDAT/GPR)
- [ ] Mock-based unit tests for `process_illumina()`, `process_agilent_2c()`, `process_agilent_1c()`, `process_nimblegen()`
