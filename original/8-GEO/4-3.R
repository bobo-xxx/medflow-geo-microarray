
# ===================================================
## GSE19274 表达矩阵清洗（probe + gene）
# ===================================================

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(GEOquery)
library(limma)

# ===================================================
## 0. 路径与参数
# ===================================================
data_dir   <- "/work/run/projects/bio-30/projects/8-GEO/data"
output_dir <- file.path(data_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gse_id <- "GSE29221"

expr_file <- file.path(
  data_dir, gse_id, "suppl",
  paste0(gse_id, "_non-normalized.txt.gz")
)

# ===================================================
## 1. 快速读取 RAW
# ===================================================
if (!file.exists(expr_file)) stop("RAW file not found")

expr <- fread(
  cmd = paste("zcat", shQuote(expr_file)),
  sep = "\t",
  header = TRUE,
  data.table = FALSE
)

# ===================================================
## 2. 提取 PROBE + SYMBOL + AVG_Signal
# ===================================================
probe_col <- intersect(colnames(expr), c("ID_REF", "PROBE_ID"))[1]
if (is.na(probe_col)) stop("No probe column")
if (!"TargetID" %in% colnames(expr)) stop("No SYMBOL column")

expr_raw <- expr %>%
  select(all_of(probe_col), TargetID, ends_with("AVG_Signal"))

colnames(expr_raw) <- sub("\\.AVG_Signal$", "", colnames(expr_raw))

# ===================================================
## 3. log2 + quantile（probe-level）
# ===================================================
expr_matrix <- expr_raw %>%
  select(-all_of(probe_col), -TargetID) %>%
  as.matrix()

expr_matrix <- expr_matrix[-((nrow(expr_matrix)-2):nrow(expr_matrix)), ]

message("Raw range: ",
        paste(range(expr_matrix, na.rm = TRUE), collapse = " ~ "))

min_val <- min(expr_matrix, na.rm = TRUE)
if (min_val <= 0) {
  expr_matrix <- expr_matrix - min_val + 1
}

expr_log  <- log2(expr_matrix)
expr_norm <- normalizeBetweenArrays(expr_log, method = "quantile")
expr_norm <- as.data.frame(expr_norm)

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
## 6. gene 层面聚合（直接用 SYMBOL）
# ===================================================
expr_norm$SYMBOL <- expr_raw$SYMBOL

expr_gene <- expr_norm %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(expr_gene) <- expr_gene$SYMBOL
expr_gene$SYMBOL <- NULL

# ===================================================
## 7. 保存结果
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
message("Probes: ", nrow(expr_probe))
message("Genes : ", nrow(expr_gene))
message("Samples: ", ncol(expr_gene))

