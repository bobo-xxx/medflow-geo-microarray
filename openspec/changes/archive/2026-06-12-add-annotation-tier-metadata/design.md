## Design

Add `annotation_tier` (1-5), `annotation_method` (source name), and `annotation_warning` (failure chain) to per-platform result metadata.

### Tracking

Each tier sets `anno_tier` and `anno_method` when it succeeds. Failed tiers append to `anno_reasons`. If Tier 5 fires, `anno_warning` contains the full failure chain.

### Result NDJSON

```json
{"level":"result","status":"success_matrix",
 "metadata":{"GPL570":{
   "platform":"GPL570","organism":"Homo sapiens",
   "n_samples":12,"n_probes":54675,
   "annotation_tier":1,"annotation_method":"fData:Gene Symbol",
   "annotation_warning":null}}}

{"level":"result","status":"success_matrix",
 "metadata":{"GPL16686":{
   "annotation_tier":5,"annotation_method":"probe_ids",
   "annotation_warning":"Tier 1: no Gene Symbol column; Tier 2: no gene_assignment column; Tier 3: GPL table unavailable; Tier 4: AnnoProbe not installed"}}}
```

### Files

| File | Change |
|---|---|
| `node/scripts/fetch.R` | Track annotation tier/method/warning in process_expression_set and process_raw_matrix |
