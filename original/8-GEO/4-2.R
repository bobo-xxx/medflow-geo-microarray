
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

gse_id <- "GSE37816"

expr_file <- file.path(
  data_dir, gse_id, "suppl",
  paste0(gse_id, "_signal_intensities.txt.gz")
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
## 2. 修正 Illumina RAW（多行说明 + 真正 header）
# ===================================================

# 去掉前 3 行说明
#expr <- expr[-c(1:4), ]

# 第 1 行是真正列名
real_colnames <- as.character(unlist(expr[1, ]))
expr <- expr[-1, ]
colnames(expr) <- real_colnames

# probe 列统一命名
colnames(expr)[1] <- "PROBE_ID"

# ===================================================
## 3. 提取 PROBE + AVG_Signal（去 Detection Pval）
# ===================================================


signal_cols <- colnames(expr)[
  !grepl("Detection", colnames(expr), ignore.case = TRUE)
]

expr_raw <- expr %>%
  select(PROBE_ID, all_of(signal_cols[-1]))

# 转数值（非常关键）
expr_raw[-1] <- lapply(expr_raw[-1], as.numeric)

# 
# # 只保留 AVG_Signal
# signal_cols <- grep("\\.AVG_Signal$", colnames(expr), value = TRUE)
# 
# stopifnot(length(signal_cols) > 0)
# 
# expr_raw <- expr %>%
#   select(PROBE_ID, all_of(signal_cols))
# 
# # 转为数值（非常关键）
# expr_raw[signal_cols] <- lapply(expr_raw[signal_cols], as.numeric)
# ===================================================
## 4. probe-level 表达矩阵 + 标准化
# ===================================================

expr_matrix <- as.matrix(expr_raw[-1])

message(
  "Raw range: ",
  paste(range(expr_matrix, na.rm = TRUE), collapse = " ~ ")
)

# Illumina RAW 常有负值 → shift
min_val <- min(expr_matrix, na.rm = TRUE)
if (min_val <= 0) {
  expr_matrix <- expr_matrix - min_val + 1
}

expr_log  <- log2(expr_matrix+1)
expr_norm <- normalizeBetweenArrays(expr_log, method = "quantile")
expr_norm <- as.data.frame(expr_norm)

rownames(expr_norm) <- expr_raw$PROBE_ID
expr_probe <- expr_norm

# ===================================================
## 5. GSM 列名修正（按顺序，最稳）
# ===================================================
gse_obj <- getGEO(gse_id, GSEMatrix = FALSE)
gsm_ids <- names(GSMList(gse_obj))

stopifnot(length(gsm_ids) == ncol(expr_probe))
colnames(expr_probe) <- gsm_ids

# ===================================================
## 6. 读取 GPL 注释（probe → gene）
# ===================================================
gpl_id <- names(GPLList(gse_obj))[1]
message("Using platform: ", gpl_id)

gpl <- getGEO(gpl_id)
gpl_tab <- Table(gpl)

gene_col <- intersect(
  colnames(gpl_tab),
  c("Gene Symbol", "GENE_SYMBOL", "Symbol",
    "gene_assignment", "GENE_NAME")
)[1]

if (is.na(gene_col)) stop("No gene symbol column in GPL")

probe2gene <- gpl_tab %>%
  select(
    probe_id    = ID,
    gene_symbol = all_of(gene_col)
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  mutate(gene_symbol = str_split(gene_symbol, " /// ")) %>%
  unnest(gene_symbol)

# ===================================================
## 7. gene-level 聚合
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
## 8. 保存结果
# ===================================================
out_dir <- file.path(output_dir, gse_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  expr_probe,
  file.path(out_dir, paste0("expr_probe_", gse_id, ".csv")),
  row.names = TRUE
)

write.csv(
  expr_gene,
  file.path(out_dir, paste0("expr_gene_", gse_id, ".csv")),
  row.names = TRUE
)

message("✅ Finished processing ", gse_id)
message("Probe-level rows: ", nrow(expr_probe))
message("Gene-level rows : ", nrow(expr_gene))
message("Samples         : ", ncol(expr_gene))

