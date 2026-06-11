# ===================================================
## Series_matrix 清洗处理（单数据集 + 背景矫正）
# ===================================================

library(dplyr)
library(tidyr)
library(stringr)
library(GEOquery)
library(Biobase)
library(limma)
library(tibble)

# ===================================================
## 0. 路径与参数
# ===================================================
data_dir   <- "/work/run/projects/bio-30/projects/8-GEO/data"
output_dir <- file.path(data_dir, "output")
platform_dir <- file.path(data_dir, "platforms")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gse_id <- "GSE19274"

expr_file <- file.path(data_dir, gse_id, "suppl", paste0(gse_id, "_non-normalized.txt.gz"))
note_file <- file.path(output_dir, "SKIPPED_REASON.txt")

# ===================================================
## 1. 读取 RAW 数据
# ===================================================
if (!file.exists(expr_file)) {
  writeLines(paste0("File not found: ", expr_file), note_file)
  stop("RAW file not found")
}

library(data.table)

expr <- fread(
  cmd = paste("zcat", shQuote(expr_file)),
  sep = "\t",
  header = TRUE,
  data.table = FALSE
)

# ===================================================
## 2. probe + AVG_Signal 提取（这是核心）
# ===================================================
probe_col <- intersect(colnames(expr), c("ID_REF", "PROBE_ID"))[1]
if (is.na(probe_col)) stop("No probe ID column found")

expr_raw <- expr %>%
  select(all_of(probe_col), ends_with("AVG_Signal"))

rownames(expr_raw) <- expr_raw[[probe_col]]
expr_raw[[probe_col]] <- NULL

colnames(expr_raw) <- sub("\\.AVG_Signal$", "", colnames(expr_raw))

# ===================================================
## 3. log2 转换 + 分位数标准化
# ===================================================
expr_matrix <- as.matrix(expr_raw)

# 建议先看一眼
message("Raw range: ", paste(range(expr_matrix, na.rm = TRUE), collapse = " ~ "))
expr_matrix <- as.matrix(expr_raw)

min_val <- min(expr_matrix, na.rm = TRUE)

if (min_val <= 0) {
  expr_matrix <- expr_matrix - min_val + 1
}

message("Shifted range: ",
        paste(range(expr_matrix, na.rm = TRUE), collapse = " ~ "))

# log2
expr_log <- log2(expr_matrix)

# 分位数标准化
expr_norm <- normalizeBetweenArrays(expr_log, method = "quantile")

expr_norm <- as.data.frame(expr_norm)

# ===================================================
## 3.1 替换 RAW 列名为 GSM ID
# ===================================================
Sys.setenv(http_proxy="http://localhost:1086",
           https_proxy="http://localhost:1086")
gse_obj <- getGEO(gse_id, GSEMatrix = FALSE)
Sys.unsetenv(c("http_proxy", "https_proxy"))

gsm_map <- lapply(GSMList(gse_obj), function(x) {
  data.frame(
    raw_name = Meta(x)$title,   # RAW 文件列名
    gsm = Meta(x)$geo_accession # GSM ID
  )
}) %>% bind_rows()

# 替换列名
colnames(expr) <- gsm_map$gsm[match(colnames(expr), gsm_map$raw_name)]

# ===================================================
## 4. 获取 GPL 注释
# ===================================================
gpl_id <- names(GPLList(gse_obj))
if (length(gpl_id) != 1) {
  writeLines(paste("Multiple or zero GPL detected:", paste(gpl_id, collapse = ", ")), note_file)
  stop("Platform detection issue")
}
gpl_id <- gpl_id[1]
message("Using platform: ", gpl_id)

gpl_soft_file <- file.path(platform_dir, gpl_id, paste0(gpl_id, ".soft.gz"))
if (!file.exists(gpl_soft_file)) {
  writeLines(paste("GPL SOFT file not found:", gpl_soft_file), note_file)
  stop("Missing GPL annotation")
}
gpl <- getGEO(filename = gpl_soft_file)
gpl_table <- Table(gpl)

# ===================================================
## 5. probe → gene 映射
# ===================================================
gene_col <- intersect(colnames(gpl_table),
                      c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment", "GENE_NAME"))[1]

if (is.na(gene_col)) {
  message("No gene symbol column found, using probe_id as gene symbol")
  probe2gene <- gpl_table %>%
    select(probe_id = ID) %>%
    mutate(gene_symbol = probe_id)
} else {
  probe2gene <- gpl_table %>%
    select(probe_id = ID, gene_symbol = all_of(gene_col)) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    mutate(gene_symbol = str_split(gene_symbol, " /// ")) %>%
    unnest(gene_symbol)
}

# ===================================================
## 6. probe → gene 聚合
# ===================================================
expr_df <- expr
expr_df$probe_id <- rownames(expr_df)

expr_gene <- expr_df %>%
  pivot_longer(-probe_id, names_to = "sample", values_to = "expr") %>%
  inner_join(probe2gene, by = "probe_id") %>%
  group_by(gene_symbol, sample) %>%
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = expr) %>%
  as.data.frame()

rownames(expr_gene) <- expr_gene$gene_symbol
expr_gene$gene_symbol <- NULL

# ===================================================
## 7. 保存 probe-level + gene-level
# ===================================================
out_gse <- file.path(output_dir, gse_id)
dir.create(out_gse, recursive = TRUE, showWarnings = FALSE)

expr_probe_file <- file.path(out_gse, paste0("expr_probe_", gse_id, ".csv"))
expr_gene_file  <- file.path(out_gse, paste0("expr_gene_", gse_id, ".csv"))

write.csv(expr_df, expr_probe_file, row.names = TRUE)
write.csv(expr_gene, expr_gene_file,  row.names = TRUE)

message("✅ Finished processing ", gse_id)
message("Probe-level saved to: ", expr_probe_file)
message("Gene-level saved to: ", expr_gene_file)

