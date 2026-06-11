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
  paste0(gse_id, "_non_normalized.txt.gz")
)

# ===================================================
## 1. 读取 RAW
# ===================================================
expr <- fread(
  cmd = paste("zcat", shQuote(expr_file)),
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "null")
)



# ===================================================
## 2. 修正 Illumina RAW header
# ===================================================
real_colnames <- as.character(unlist(expr[1, ]))
expr <- expr[-1, ]
colnames(expr) <- real_colnames

colnames(expr)[1] <- "PROBE_ID"
# 
# # 自动识别 gene symbol 列
# gene_col <- intersect(
#   colnames(expr),
#   c("TargetID", "SYMBOL", "Gene Symbol", "GENE_SYMBOL")
# )[1]
# if (is.na(gene_col)) stop("❌ No gene symbol column")



# 先找所有 Signal 列
signal_cols <- grep("\\.Signal$", colnames(expr), value = TRUE)

# 按样本名分组，去掉 A/B 后缀
sample_names <- unique(sub("\\.[AB]$", "", signal_cols))

# 对每个样本的 A/B 列取平均
expr_avg <- sapply(sample_names, function(s) {
  cols <- grep(paste0("^", s, "\\.[AB]$"), colnames(expr), value = TRUE)
  rowMeans(expr[, cols], na.rm = TRUE)
})

expr_avg <- as.data.frame(expr_avg)
rownames(expr_avg) <- expr$PROBE_ID


# ===================================================
## 3. 提取 PROBE + AVG_Signal
# ===================================================
expr_raw <- expr %>%
  select(PROBE_ID, all_of(gene_col), ends_with("AVG_Signal"))

# 去掉 .AVG_Signal 后缀
colnames(expr_raw) <- sub("\\.AVG_Signal$", "", colnames(expr_raw))

# 转 numeric（除 probe/gene 列）
expr_raw[-c(1,2)] <- lapply(expr_raw[-c(1,2)], as.numeric)

# ===================================================
## 4. probe-level log2 + quantile 标准化
# ===================================================
expr_matrix <- as.matrix(expr_raw[-c(1,2)])


message("Raw range: ",
        paste(range(expr_matrix, na.rm = TRUE), collapse = " ~ "))

# Illumina RAW 允许负值 → shift
min_val <- min(expr_matrix, na.rm = TRUE)
if (min_val <= 0) {
  expr_matrix <- expr_matrix - min_val + 1
}

expr_log  <- log2(expr_matrix+1)
expr_norm <- normalizeBetweenArrays(expr_log, method = "quantile")
expr_norm <- as.data.frame(expr_norm)

rownames(expr_norm) <- expr_raw$PROBE_ID[1:nrow(expr_norm)]

expr_probe <- expr_norm

# ===================================================
## 4. 替换列名为 GSM ID
# ===================================================
gse_obj <- getGEO(gse_id, GSEMatrix = FALSE)

gsm_map <- lapply(GSMList(gse_obj), function(x) {
  data.frame(
    raw_name = Meta(x)$title,
    gsm      = Meta(x)$geo_accession,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

new_names <- gsm_map$gsm[
  match(colnames(expr_norm), gsm_map$raw_name)
]

colnames(expr_norm) <- ifelse(
  is.na(new_names),
  colnames(expr_norm),
  new_names
)

# ===================================================
## 5. 保存 probe-level 矩阵
# ===================================================
expr_probe <- expr_norm
rownames(expr_probe) <- expr_raw[[probe_col]]


# ===================================================
## 6. gene-level 聚合（用 gene symbol 列）
# ===================================================
expr_probe$SYMBOL <- expr_raw[[gene_col]]  # 添加基因列

expr_gene <- expr_probe %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%        # 去掉没有基因符号的探针
  group_by(SYMBOL) %>%                           # 按基因聚合
  summarise(across(where(is.numeric), mean, na.rm = TRUE),  # 求平均
            .groups = "drop") %>%
  as.data.frame()

rownames(expr_gene) <- expr_gene$SYMBOL
expr_gene$SYMBOL <- NULL

# ===================================================
## 7. 保存结果
# ===================================================
out_dir <- file.path(output_dir, gse_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 保存 probe-level矩阵
write.csv(expr_probe,
          file.path(out_dir, paste0("expr_probe_", gse_id, ".csv")),
          row.names = TRUE)

# 保存 gene-level矩阵
write.csv(expr_gene,
          file.path(out_dir, paste0("expr_gene_", gse_id, ".csv")),
          row.names = TRUE)

# ===================================================
## 8. 信息输出
# ===================================================
message("✅ Finished processing ", gse_id)
message("Probe-level rows: ", nrow(expr_probe))
message("Gene-level rows : ", nrow(expr_gene))
message("Samples         : ", ncol(expr_gene))

