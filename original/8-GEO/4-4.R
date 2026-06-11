#模糊匹配ID"GSE28894"

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(GEOquery)
library(limma)

# ===================================================
## 0. 参数
# ===================================================
data_dir   <- "/work/run/projects/bio-30/projects/8-GEO/data"
output_dir <- file.path(data_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gse_id <- "GSE28894"

expr_file <- file.path(
  data_dir, gse_id, "suppl",
  paste0(gse_id, "_non-normalized.txt.gz")
)

# ===================================================
## 1. 读取 RAW
# ===================================================
if (!file.exists(expr_file)) stop("RAW file not found")

expr <- fread(
  cmd = paste("zcat", shQuote(expr_file)),
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

# ===================================================
## 2. 提取 probe + 表达量（去 Detection_Pval）
# ===================================================
probe_col <- intersect(colnames(expr), c("ID_REF", "PROBE_ID"))[1]
if (is.na(probe_col)) stop("No probe column found")

expr_cols <- colnames(expr)
expr_cols <- expr_cols[!grepl("Detection_Pval|Detection Pval", expr_cols)]
expr_cols <- setdiff(expr_cols, probe_col)

expr_raw <- expr %>%
  select(all_of(probe_col), all_of(expr_cols))

expr_raw[expr_cols] <- lapply(expr_raw[expr_cols], as.numeric)

# ===================================================
## 3. probe-level 表达矩阵（log2 + quantile）
# ===================================================
expr_matrix <- as.matrix(expr_raw[expr_cols])

message(
  "Raw range: ",
  paste(range(expr_matrix, na.rm = TRUE), collapse = " ~ ")
)

expr_log  <- log2(expr_matrix + 1)
expr_norm <- normalizeBetweenArrays(expr_log, method = "quantile")
expr_norm <- as.data.frame(expr_norm)

rownames(expr_norm) <- expr_raw[[probe_col]]
expr_probe <- expr_norm

# ===================================================
## 4. ⭐ 正确替换列名为 GSM（按顺序）
# ===================================================
gse_obj <- getGEO(gse_id, GSEMatrix = FALSE)
gsm_ids <- names(GSMList(gse_obj))

if (length(gsm_ids) != ncol(expr_probe)) {
  stop("GSM number != expression matrix columns")
}

colnames(expr_probe) <- gsm_ids

# ===================================================
## 5. 读取 GPL 注释（probe → gene）
# ===================================================
gpl_id <- names(GPLList(gse_obj))[1]
message("Using platform: ", gpl_id)

gpl <- getGEO(gpl_id)
gpl_table <- Table(gpl)

gene_col <- intersect(
  colnames(gpl_table),
  c("Gene Symbol", "GENE_SYMBOL", "Symbol",
    "gene_assignment", "GENE_NAME")
)[1]

if (is.na(gene_col)) stop("No gene symbol column in GPL")

probe2gene <- gpl_table %>%
  select(
    probe_id    = ID,
    gene_symbol = all_of(gene_col)
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  mutate(gene_symbol = str_split(gene_symbol, " /// ")) %>%
  unnest(gene_symbol)

# ===================================================
## 6. gene-level 聚合
# ===================================================
expr_probe$probe_id <- rownames(expr_probe)

expr_gene <- expr_probe %>%
  pivot_longer(-probe_id, names_to = "sample", values_to = "expr") %>%
  inner_join(probe2gene, by = "probe_id") %>%
  group_by(gene_symbol, sample) %>%
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = expr) %>%
  as.data.frame()

rownames(expr_gene) <- expr_gene$gene_symbol
expr_gene$gene_symbol <- NULL

# ===================================================
## 7. 保存结果
# ===================================================
out_dir <- file.path(output_dir, gse_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(expr_probe,
          file.path(out_dir, paste0("expr_probe_", gse_id, ".csv")),
          row.names = TRUE)

write.csv(expr_gene,
          file.path(out_dir, paste0("expr_gene_", gse_id, ".csv")),
          row.names = TRUE)

message("✅ Finished processing ", gse_id)
message("Probe-level rows: ", nrow(expr_probe))
message("Gene-level rows : ", nrow(expr_gene))
message("Samples         : ", ncol(expr_gene))
