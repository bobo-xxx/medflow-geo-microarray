### Task 1: Track annotation tier/method/warning

- [x] `process_expression_set`: add anno_tier, anno_method, anno_reasons tracking
- [x] Each tier sets tier+method on success, appends reason on failure
- [x] Store annotation_tier, annotation_method, annotation_warning in result$metadata
- [x] `process_raw_matrix`: set NA tier with deferred warning

### Task 2: Verify NDJSON output

- [x] Full test suite passes
- [ ] End-to-end test on GSE318047 (Tier 1) and GSE156508 (Tier 3+)
