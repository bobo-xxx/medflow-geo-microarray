# Error Codes Reference

This document defines all error codes used in the geo-microarray-fetch skill, along with their meanings and recommended recovery actions.

## Error Codes

| Code | Message | Action |
|------|---------|--------|
| E001 | No suppl directory | Skip GSE - No supplementary files available, cannot proceed with this dataset |
| E002 | Methylation BPM present | Skip GSE - Methylation data not supported, focus on expression data only |
| E003 | No series matrix | Fallback to raw - Try downloading raw CEL files instead |
| E004 | No raw files | Return metadata only - No data files available, return GSE metadata |
| E005 | Corrupted CEL file | Log and skip file - Continue processing other files, log the corrupted file |
| E006 | GPL not found | Auto-download - Retrieve GPL platform annotation from GEO database |
| E007 | Rate limit exceeded | Retry with backoff - Implement exponential backoff before retrying request |

## Recovery Strategies

### E001: No suppl directory
- **Severity:** Warning
- **Recovery:** Skip current GSE, continue to next dataset
- **Impact:** Complete data loss for this GSE
- **Prevention:** None - some GEO datasets simply lack supplementary files

### E002: Methylation BPM present
- **Severity:** Warning
- **Recovery:** Skip current GSE, continue to next dataset
- **Impact:** Methylation data excluded from processing
- **Prevention:** None - by design, this skill focuses on expression data

### E003: No series matrix
- **Severity:** Recoverable
- **Recovery:** Attempt fallback to raw CEL file download
- **Impact:** May increase processing time and storage requirements
- **Prevention:** Check series matrix availability before attempting download

### E004: No raw files
- **Severity:** Warning
- **Recovery:** Return GSE metadata without expression data
- **Impact:** Incomplete data for this GSE
- **Prevention:** Check for both series matrix and raw file availability before processing

### E005: Corrupted CEL file
- **Severity:** Recoverable
- **Recovery:** Log file path and error, skip to next file
- **Impact:** Partial data loss for this GSE
- **Prevention:** Validate file integrity after download, consider checksum verification

### E006: GPL not found
- **Severity:** Recoverable
- **Recovery:** Auto-download GPL from GEO database
- **Impact:** Additional network request, slight delay
- **Prevention:** Cache commonly used GPLs locally

### E007: Rate limit exceeded
- **Severity:** Recoverable
- **Recovery:** Implement exponential backoff (start with 5s, double up to 60s max)
- **Impact:** Increased processing time during rate-limiting
## Validation Thresholds

### EXTREME_VALUE_THRESHOLD
- **Value:** 1e50
- **Purpose:** Upper bound for detecting extreme/outlier values in gene expression data
- **Rationale:** Expression data typically ranges from -20 to +20 on log2 scale. A threshold of 1e50 provides a reasonable upper bound to detect corrupted data while avoiding false positives on legitimate high-expression values.
- **Usage:** Applied in `validate_gene_expression()` and Tier 2/Tier 3 matrix validation
- **Recovery:** Values exceeding this threshold trigger rejection of the dataset, forcing fallback to lower-tier data sources.
