## ===================================================
## 0. 基础设置
## ===================================================

library(dplyr)
library(tidyr)
library(stringr)
library(oligo)
library(R.utils)
library(Biobase)
library(GEOquery)
library(tibble)
library(pd.primeview)


## ===================================================
## 1. 路径与参数
## ===================================================
data_dir     <- "/work/run/projects/bio-30/projects/8-GEO/data"
platform_dir <- file.path(data_dir, "platforms")
output_dir   <- file.path(data_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Affymetrix microarray 数据集
gse_affy <- gse_keep_summary %>%
  filter(inferred_type == "Affymetrix_microarray") %>%
  pull(GSE)

## ===================================================
## 2. 批量循环处理 GSE
## ===================================================
for (gse_id in gse_affy) {
  
  out_gse <- file.path(output_dir, gse_id)
  expr_gene_file  <- file.path(out_gse, paste0("expr_gene_", gse_id, ".csv"))
  expr_probe_file <- file.path(out_gse, paste0("expr_probe_", gse_id, ".csv"))
  
  # ---- 判断是否已经存在结果，跳过 ----
  if (file.exists(expr_gene_file) && file.exists(expr_probe_file)) {
    message("✅ Result files already exist, skipping ", gse_id)
    next
  }
  
  message("=== Processing ", gse_id, " ===")
  
  suppl_dir <- file.path(data_dir, gse_id, "suppl")
  raw_dir   <- file.path(suppl_dir, "RAW")
  dir.create(out_gse, recursive = TRUE, showWarnings = FALSE)
  
  ## ===================================================
  ## 3. 解压 RAW.tar
  ## ===================================================
  tar_files <- list.files(suppl_dir, pattern = "_RAW\\.tar$", full.names = TRUE)
  if (length(tar_files) == 0) {
    message("❌ No RAW tar file found for ", gse_id)
    next
  }
  
  if (!dir.exists(raw_dir)) dir.create(raw_dir)
  if (length(list.files(raw_dir)) == 0) untar(tar_files[1], exdir = raw_dir)
  
  ## ===================================================
  ## 4. 查找 CEL 文件
  ## ===================================================
  cel_files <- list.files(raw_dir, pattern = "\\.CEL(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  if (length(cel_files) == 0) {
    message("❌ No CEL files found for ", gse_id)
    next
  }
  
  ## ===================================================
  ## 5. 过滤损坏 CEL 文件
  ## ===================================================
  good_flags <- sapply(cel_files, function(f) {
    tryCatch({ read.celfiles(f); TRUE }, 
             error = function(e) { message("❌ Bad CEL: ", basename(f)); FALSE })
  })
  
  good_cel_files <- cel_files[good_flags]
  bad_cel_files  <- basename(cel_files[!good_flags])
  
  if (length(bad_cel_files) > 0) {
    writeLines(bad_cel_files, file.path(out_gse, paste0("bad_cel_", gse_id, ".txt")))
  }
  
  if (length(good_cel_files) == 0) {
    message("❌ All CEL files are bad for ", gse_id)
    next
  }
  
  ## ===================================================
  ## 6. 读取 CEL + RMA
  ## ===================================================
  raw_data <- read.celfiles(good_cel_files)
  eset     <- rma(raw_data)
  expr_matrix <- exprs(eset)
  colnames(expr_matrix) <- basename(colnames(expr_matrix)) |> str_remove("\\.CEL(\\.gz)?$")
  
  ## ===================================================
  ## 7. 获取平台 GPL
  ## ===================================================
  Sys.setenv(http_proxy="http://localhost:1086", https_proxy="http://localhost:1086")
  gse_obj <- getGEO(gse_id, GSEMatrix = FALSE)
  gpl_id  <- names(GPLList(gse_obj))
  Sys.unsetenv(c("http_proxy", "https_proxy"))
  
  if (length(gpl_id) != 1) {
    message("❌ Multiple or zero GPL found for ", gse_id)
    next
  }
  
  gpl_id <- gpl_id[1]
  message("Using platform: ", gpl_id)
  
  ## ===================================================
  ## 8. 读取本地 GPL SOFT
  ## ===================================================
  gpl_soft_file <- file.path(platform_dir, gpl_id, paste0(gpl_id, ".soft.gz"))
  if (!file.exists(gpl_soft_file)) {
    message("❌ Platform SOFT not found: ", gpl_soft_file)
    next
  }
  
  gpl <- getGEO(filename = gpl_soft_file)
  gpl_table <- Table(gpl)
  
  ## ===================================================
  ## 9. 构建 probe → gene 映射
  ## ===================================================
  gene_col <- intersect(colnames(gpl_table), c("Gene Symbol","GENE_SYMBOL","Symbol","gene_assignment"))[1]
  if (is.na(gene_col)) {
    message("❌ No gene symbol column found for ", gse_id)
    next
  }
  
  probe2gene <- gpl_table %>%
    select(probe_id=ID, gene_symbol=all_of(gene_col)) %>%
    filter(!is.na(gene_symbol), gene_symbol!="") %>%
    mutate(gene_symbol=str_split(gene_symbol," /// ")) %>%
    unnest(gene_symbol)
  
  ## ===================================================
  ## 10. probe → gene（取 mean）
  ## ===================================================
  expr_df <- as.data.frame(expr_matrix)
  expr_df$probe_id <- rownames(expr_df)
  
  expr_long <- expr_df %>%
    pivot_longer(cols=-probe_id, names_to="sample", values_to="expr")
  
  expr_annotated <- expr_long %>% inner_join(probe2gene, by="probe_id")
  
  expr_gene <- expr_annotated %>%
    group_by(gene_symbol, sample) %>%
    summarise(expr=mean(expr, na.rm=TRUE), .groups="drop")
  
  expr_gene_matrix <- expr_gene %>%
    pivot_wider(names_from=sample, values_from=expr) %>%
    as.data.frame()
  rownames(expr_gene_matrix) <- expr_gene_matrix$gene_symbol
  expr_gene_matrix$gene_symbol <- NULL
  
  ## ===================================================
  ## 11. 保存结果
  ## ===================================================
  write.csv(expr_matrix, file=expr_probe_file, row.names=TRUE)
  write.csv(expr_gene_matrix, file=expr_gene_file, row.names=TRUE)
  
  message("✅ Finished processing ", gse_id)
}
