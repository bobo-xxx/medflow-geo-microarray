## Design

Three annotation quality fixes from 27-dataset agentic test run.

### Issue 1: entg| prefix on gene symbols

`clean_gene_symbols()` strips `entg|`, `ref|`, `gb|` prefixes from BioC annotation DB and GPL table output. Applied at Tier 3 (gpl_table) and Tier 4 (annotate_with_bioc_db).

### Issue 2: Agilent GPL19072 supplementary file

`try_get_gpl_suppl()` downloads the GPL supplementary annotation file when Table() lacks gene symbols. Uses `download.file` with `mode="wb"`, falls back to gunzip system command. Uses NAME column for probe IDs (Agilent suppl files use row numbers as ID). 1GB size limit.

### Issue 3: Tier 5 redundant gene matrix

Gene CSV output suppressed when `anno_tier == 5` (probe IDs as gene names). SKILL.md output contract updated with `condition: annotation_tier <= 4`.
