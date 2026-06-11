# ===================================================
# 稳健版循环处理 Illumina IDAT 数据集（列名改为 GSM）
# 自动 log2 + 分位数标准化 + probe -> gene 汇总
# ===================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(R.utils)
  library(Biobase)
  library(GEOquery)
  library(tibble)
  library(illuminaio)
  library(limma)  # normalizeBetweenArrays
})

options(stringsAsFactors = FALSE)

# -------------------------
# 判断表达矩阵类型
# -------------------------
detect_expr_type <- function(x) {
  q <- quantile(x, probs = c(0.01, 0.5, 0.99), na.rm = TRUE)
  mean_x <- mean(x, na.rm = TRUE)
  if (q[3] > 100) return("raw")
  if (abs(mean_x) < 0.5 && q[3] < 10 && q[1] > -10) return("centered")
  return("log")
}

# -------------------------
# 背景矫正 + 标准化
# -------------------------
normalize_expr_matrix <- function(x) {
  type <- detect_expr_type(x)
  message("Detected expression type: ", type)
  
  # 负值置零
  x[x < 0] <- 0
  
  # raw 数据 log2
  if (type == "raw") x <- log2(x + 1)
  
  # 分位数标准化
  x <- normalizeBetweenArrays(x, method = "quantile")
  return(x)
}

# ===================================================
# 路径设置
# ===================================================
data_dir     <- "/work/run/projects/bio-30/projects/8-GEO/data"
platform_dir <- file.path(data_dir, "platforms")
output_dir   <- file.path(data_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ===================================================
# 待处理 GSE（假设 gse_keep_summary 已经定义）
# ===================================================
gse_illumina <- gse_keep_summary %>%
  filter(inferred_type == "Illumina_expression_array") %>%
  pull(GSE)

# ===================================================
# 循环处理每个 GSE
# ===================================================
for (gse_id in gse_illumina) {
  
  message("=== Processing ", gse_id, " ===")
  
  out_gse <- file.path(output_dir, gse_id)
  dir.create(out_gse, recursive = TRUE, showWarnings = FALSE)
  
  expr_gene_file  <- file.path(out_gse, paste0("expr_gene_", gse_id, ".csv"))
  expr_probe_file <- file.path(out_gse, paste0("expr_probe_", gse_id, ".csv"))
  note_file       <- file.path(out_gse, "SKIPPED_REASON.txt")
  
  if (file.exists(expr_gene_file) && file.exists(expr_probe_file)) {
    message("✅ Already processed, skipping ", gse_id)
    next
  }
  
  # -----------------------------------------
  # RAW 文件准备
  # -----------------------------------------
  suppl_dir <- file.path(data_dir, gse_id, "suppl")
  raw_dir   <- file.path(suppl_dir, "RAW")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 解压 .tar / .zip
  archives <- list.files(suppl_dir, pattern = "\\.(tar|zip)$", full.names = TRUE)
  if (length(archives) > 0 && length(list.files(raw_dir)) == 0) {
    for (a in archives) {
      tryCatch({
        if (grepl("\\.tar$", a, ignore.case = TRUE)) untar(a, exdir = raw_dir)
        if (grepl("\\.zip$", a, ignore.case = TRUE)) unzip(a, exdir = raw_dir)
      }, error = function(e) message("❌ Archive extract failed: ", basename(a)))
    }
  }
  
  # 解压 .idat.gz
  idat_gz_files <- list.files(raw_dir, pattern = "\\.idat\\.gz$", full.names = TRUE, recursive = TRUE)
  for (gz in idat_gz_files) {
    tryCatch(gunzip(gz, remove = FALSE, overwrite = TRUE), error = function(e){})
  }
  
  # IDAT 文件
  idat_files <- list.files(raw_dir, pattern = "\\.idat$", full.names = TRUE, recursive = TRUE)
  if (length(idat_files) == 0) {
    writeLines("❌ No IDAT files found.", note_file)
    message("❌ No IDAT files found for ", gse_id)
    next
  }
  
  # -----------------------------------------
  # 读取 IDAT 文件
  # -----------------------------------------
  idatlist <- tryCatch(lapply(idat_files, readIDAT), error = function(e) NULL)
  if (is.null(idatlist)) {
    writeLines("❌ Failed to read IDAT files.", note_file)
    next
  }
  
  exprs <- sapply(idatlist, function(x) x$Quants$MeanBinData)
  rownames(exprs) <- idatlist[[1]]$Quants$CodesBinData
  colnames(exprs) <- sapply(idatlist, function(x) paste0(x$Barcode, "_", x$Section))
  
  # -----------------------------------------
  # 替换列名为 GSM
  # -----------------------------------------
  gsm_vec <- str_extract(basename(idat_files), "^GSM\\d+")
  if (any(is.na(gsm_vec)) || length(gsm_vec) != ncol(exprs)) {
    writeLines("❌ GSM extraction/count mismatch.", note_file)
    next
  }
  colnames(exprs) <- gsm_vec
  if (any(duplicated(colnames(exprs)))) {
    writeLines("❌ Duplicated GSM IDs detected.", note_file)
    next
  }
  
  # -----------------------------------------
  # 背景矫正 + 分位数标准化
  # -----------------------------------------
  exprs <- normalize_expr_matrix(exprs)
  
  # -----------------------------------------
  # 读取 GPL 注释并映射 probe -> gene
  # -----------------------------------------
  gpl_id <- tryCatch({
    gse_obj <- getGEO(gse_id, GSEMatrix = FALSE)
    names(GPLList(gse_obj))[1]
  }, error = function(e) NA)
  
  gpl_soft_file <- file.path(platform_dir, gpl_id, paste0(gpl_id, ".soft.gz"))
  if (!file.exists(gpl_soft_file)) {
    writeLines(paste("GPL SOFT file not found:", gpl_soft_file), note_file)
    next
  }
  
  gpl <- getGEO(filename = gpl_soft_file)
  gpl_table <- Table(gpl)
  
  gene_col <- intersect(colnames(gpl_table), c("Gene Symbol","GENE_SYMBOL","Symbol","gene_assignment"))[1]
  if (is.na(gene_col)) {
    writeLines("No gene symbol column found in GPL annotation table.", note_file)
    next
  }
  
  probe2gene <- gpl_table %>%
    select(probe_id = ID, gene_symbol = all_of(gene_col)) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    mutate(gene_symbol = str_split(gene_symbol, " /// ")) %>%
    tidyr::unnest(gene_symbol)
  
  # probe -> gene 汇总
  expr_df <- as.data.frame(exprs)
  expr_df$numeric_id <- as.numeric(rownames(expr_df))
  
  probe_map <- gpl_table %>%
    select(probe_id = ID, numeric_id = Array_Address_Id) %>%
    filter(!is.na(numeric_id))
  
  expr_df <- expr_df %>%
    left_join(probe_map, by = "numeric_id") %>%
    filter(!is.na(probe_id)) %>%
    select(probe_id, everything(), -numeric_id)
  
  expr_gene <- expr_df %>%
    pivot_longer(-probe_id, names_to = "sample", values_to = "expr") %>%
    inner_join(probe2gene, by = "probe_id") %>%
    group_by(gene_symbol, sample) %>%
    summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = expr) %>%
    as.data.frame()
  
  rownames(expr_gene) <- expr_gene$gene_symbol
  expr_gene$gene_symbol <- NULL
  
  # -----------------------------------------
  # 保存结果
  # -----------------------------------------
  write.csv(exprs, expr_probe_file, row.names = TRUE)
  write.csv(expr_gene, expr_gene_file, row.names = TRUE)
  
  message("✅ Finished processing ", gse_id)
}






