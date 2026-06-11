# 8-GEO: GEO Data Cleaning Workflow

**Purpose:** Clean and normalize microarray expression data from GEO datasets.

## WORKFLOW PIPELINE

```
1-预处理清洗.R → 2-判断数据类型.R → 3-数据清洗*.R → 4-*.R
     ↓                  ↓                  ↓              ↓
   Filter GSE      Detect platform     Clean data    Post-process
   (skip BPM)      (CEL/IDAT/GPR/TXT)
```

## SCRIPT REFERENCE

| Script | Input | Output | Function |
|--------|-------|--------|----------|
| `1-预处理清洗.R` | GSE directories | `GSE_suppl_judgement_summary.csv` | Filter out methylation BPM files |
| `2-判断数据类型.R` | Filtered GSE | `GSE_keep_file_composition_summary_final.csv` | Infer array platform type |
| `3-数据清洗1.R` | CEL files | `output/GSE*/expr_*.csv` | Affymetrix RMA + probe2gene |
| `3-数据清洗IDAT.R` | IDAT files | `output/GSE*/expr_*.csv` | Illumina normalization |
| `3-数据清洗GPR.R` | GPR files | `outputT/GSE*/expr_*.csv` | Agilent/Genepix processing |
| `3-数据清洗TXT.R` | Series matrix | `output/GSE*/expr_*.csv` | TXT/CSV expression files |

## DATA STRUCTURE

```
data/
├── GSE*/suppl/           # Raw supplementary files
│   ├── filelist.txt      # GEO file manifest
│   ├── RAW/              # Extracted raw files
│   └── *_RAW.tar         # Archived raw data
├── platforms/GPL*/       # Platform annotation (SOFT files)
└── output/GSE*/          # Cleaned expression matrices
    ├── expr_probe_*.csv  # Probe-level expression
    └── expr_gene_*.csv   # Gene-level expression
```

## KEY PATTERNS

**Expression type detection:**
```r
detect_expr_type <- function(x) {
  q <- quantile(x, probs = c(0.01, 0.5, 0.99))
  if (q[3] > 100) return("raw")      # Needs log2
  if (abs(mean(x)) < 0.5) return("centered")  # Already processed
  return("log")                      # Log-transformed
}
```

**Probe-to-gene aggregation:**
```r
probe2gene <- gpl_table %>%
  select(probe_id = ID, gene_symbol = `Gene Symbol`) %>%
  mutate(gene_symbol = str_split(gene_symbol, " /// ")) %>%
  unnest(gene_symbol)

expr_gene <- expr_df %>%
  pivot_longer(-probe_id, names_to = "sample", values_to = "expr") %>%
  inner_join(probe2gene, by = "probe_id") %>%
  group_by(gene_symbol, sample) %>%
  summarise(expr = mean(expr, na.rm = TRUE))
```

## PLATFORM DETECTION

| File Pattern | Inferred Type |
|--------------|---------------|
| `*.CEL(.gz)?` | Affymetrix_microarray |
| `*.GPR(.gz)?` | Agilent_microarray |
| `*.idat` + `*.bpm` | Methylation_array (SKIP) |
| `*.idat` only | Illumina_expression_array |
| `*.txt` only | TXT_array |

## ANTI-PATTERNS

- **SKIP methylation data**: BPM files indicate methylation arrays, not expression
- **Handle negative values**: Shift to non-negative before log2 transform
- **Filter bad CEL files**: Some CEL files may be corrupted, validate before RMA
