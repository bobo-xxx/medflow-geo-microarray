### Task 1: Fix entg| prefix on gene symbols

- [x] `clean_gene_symbols()` strips entg|/ref|/gb| prefixes
- [x] Applied at Tier 3 (gpl_table) and Tier 4 (BioC DB) output
- [x] 4 TDD tests

### Task 2: GPL supplementary file fallback

- [x] `parse_gpl_suppl_soft()` extracts ID→GENE_SYMBOL from SOFT format
- [x] `try_get_gpl_suppl()` downloads suppl file with 1GB limit, caches results
- [x] FTP download fix: check curl status_code, fallback to download.file with mode="wb"
- [x] gunzip fallback: system gunzip when R gzfile() fails
- [x] NAME column preference for Agilent probe IDs
- [x] 3 TDD tests

### Task 3: Suppress redundant gene matrix at Tier 5

- [x] Gene aggregation guarded by `anno_tier < 5`
- [x] SKILL.md output contract: `condition: annotation_tier <= 4`
- [x] 2 TDD tests

### Task 4: Verify with real data

- [x] 325 tests pass (0 FAIL)
- [x] 26/27 datasets successful in agentic test run
- [x] entg| prefix eliminated on GSE102541
- [x] GPL19072 correctly falls to Tier 5 (manufacturer never populated annotation)
