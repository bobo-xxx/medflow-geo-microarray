## Why

The result NDJSON line reports `status: "success_matrix"` with no indication of which gene annotation tier succeeded or why earlier tiers failed. Agents consuming this output cannot determine data quality — they can't tell whether gene symbols came from fData (Tier 1, high quality) or probe IDs (Tier 5, low quality).

## What Changes

- **Modified**: `process_expression_set()` in fetch.R — capture `annotation_tier` (1-5), `annotation_method` (source name), and `annotation_warning` (why lower tiers failed) in the per-platform metadata
- **Modified**: Result NDJSON lines — include annotation metadata

## Capabilities

### Modified Capabilities
- `gene-annotation`: Report annotation tier, method, and failure reasons in result metadata
