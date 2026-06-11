### Task 1: Add process_illumina_txt() to raw.R

- [x] `detect_raw_type()` — recognize TargetID/AVG_Signal in TXT headers
- [x] `process_illumina_txt()` — parse interleaved columns, build EListRaw, neqc
- [x] `process_raw_files()` routing — dispatch `illumina_txt`

### Task 2: Add tests

- [x] Detection test for `illumina_txt` type
- [x] Unit tests: interleaved column parsing, empty file list, requireNamespace error
- [x] Processor existence check

### Task 3: Real data verification

- [x] GSE100155 non-normalized TXT: 48803 × 124, log2, CV=7.1e-07, 0 NaN
