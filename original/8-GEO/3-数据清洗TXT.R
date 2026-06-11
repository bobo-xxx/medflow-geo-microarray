# ===================================================
## TXT_array 循环版本（多平台适配 + 非负矩阵输出）
# ===================================================

library(dplyr)
library(tidyr)
library(stringr)
library(GEOquery)
library(Biobase)
library(tibble)

# ===================================================
## 0. 检测类型 + 统一处理函数
# ===================================================
detect_expr_type <- function(x) {
  q <- quantile(x, probs = c(0.01, 0.5, 0.99), na.rm = TRUE)
  mean_x <- mean(x, na.rm = TRUE)
  
  if (q[3] > 100) return("raw")
  if (abs(mean_x) < 0.5 && q[3] < 10 && q[1] > -10) return("centered")
  return("log")
}

normalize_expr_matrix <- function(x) {
  type <- detect_expr_type(x)
  message("Detected expression type: ", type)
  
  if (type == "raw") x <- log2(x + 1)
  
  # 保证矩阵非负
  min_val <- min(x, na.rm = TRUE)
  if (min_val < 0) x <- x - min_val
  return(x)
}

# ===================================================
## 1. 路径与参数
# ===================================================
data_dir     <- "/work/run/projects/bio-30/projects/8-GEO/data"
platform_dir <- file.path(data_dir, "platforms")
output_dir   <- file.path(data_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## 仅处理 TXT_array
gse_txt <- gse_keep_summary %>%
  filter(inferred_type == "TXT_array") %>%
  pull(GSE)

# ===================================================
## 2. 批量循环处理 GSE（支持多平台）
# ===================================================
for (gse_id in gse_txt) {
  
  out_gse <- file.path(output_dir, gse_id)
  dir.create(out_gse, recursive = TRUE, showWarnings = FALSE)
  
  note_file <- file.path(out_gse, "SKIPPED_REASON.txt")
  
  ## ---------- 已完成则跳过 ---------- 
  if (file.exists(expr_gene_file) && file.exists(expr_probe_file)) {
    message("✅ Already processed, skipping ", gse_id)
    next
  }
  
  message("=== Processing ", gse_id, " ===")
  
  
  # ===================================================
  ## 3. 读取 Series Matrix
  # ===================================================
  Sys.setenv(http_proxy="http://localhost:1086",
             https_proxy="http://localhost:1086")
  
  gse_matrix <- tryCatch(getGEO(gse_id, GSEMatrix = TRUE), error = function(e) NULL)
  Sys.unsetenv(c("http_proxy", "https_proxy"))
  
  if (is.null(gse_matrix)) {
    writeLines("Failed to load GSE Matrix.", note_file)
    next
  }
  
  # ===================================================
  ## 4. 遍历每个平台 ExpressionSet
  # ===================================================
  for (i in seq_along(gse_matrix)) {
    
    eset <- gse_matrix[[i]]
    gpl_id <- annotation(eset)
    
    # 文件名加上平台后缀
    expr_gene_file  <- file.path(out_gse, paste0("expr_gene_", gse_id, "_", gpl_id, ".csv"))
    expr_probe_file <- file.path(out_gse, paste0("expr_probe_", gse_id, "_", gpl_id, ".csv"))
    
    if (file.exists(expr_gene_file) && file.exists(expr_probe_file)) {
      message("✅ Already processed, skipping ", gse_id, "_", gpl_id)
      next
    }
    
    expr_matrix <- exprs(eset)
    if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
      writeLines("Expression matrix is empty.", note_file)
      next
    }
    
    # 统一非负矩阵处理
    expr_matrix <- normalize_expr_matrix(expr_matrix)
    colnames(expr_matrix) <- make.names(colnames(expr_matrix), unique = TRUE)
    
    message("Using platform: ", gpl_id)
    
    # ===================================================
    ## 5. 读取 GPL SOFT（只从 platform_dir）
    # ===================================================
    gpl_soft_file <- file.path(platform_dir, gpl_id, paste0(gpl_id, ".soft.gz"))
    if (!file.exists(gpl_soft_file)) {
      writeLines(paste("GPL SOFT file not found:", gpl_soft_file), note_file)
      next
    }
    
    gpl <- getGEO(filename = gpl_soft_file)
    gpl_table <- Table(gpl)
    if (!"ID" %in% colnames(gpl_table)) {
      writeLines("GPL table has no ID column.", note_file)
      next
    }
    
    # ===================================================
    ## 6. probe → gene 映射（兼容无 gene symbol 平台）
    # ===================================================
    gene_col <- intersect(
      colnames(gpl_table),
      c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment", "GENE_NAME")
    )[1]
    
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
    
    if (nrow(probe2gene) == 0) {
      writeLines("Probe-to-gene mapping is empty.", note_file)
      next
    }
    
    # ===================================================
    ## 7. probe → gene（mean 聚合）
    # ===================================================
    expr_df <- as.data.frame(expr_matrix)
    expr_df$probe_id <- rownames(expr_df)
    
    expr_gene <- expr_df %>%
      pivot_longer(-probe_id, names_to = "sample", values_to = "expr") %>%
      inner_join(probe2gene, by = "probe_id") %>%
      group_by(gene_symbol, sample) %>%
      summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = sample, values_from = expr) %>%
      as.data.frame()
    
    if (nrow(expr_gene) == 0) {
      writeLines("Gene-level expression matrix is empty.", note_file)
      next
    }
    
    rownames(expr_gene) <- expr_gene$gene_symbol
    expr_gene$gene_symbol <- NULL
    
    # ===================================================
    ## 8. 保存
    # ===================================================
    write.csv(expr_matrix, expr_probe_file, row.names = TRUE)
    write.csv(expr_gene,   expr_gene_file,  row.names = TRUE)
    
    message("✅ Finished processing ", gse_id, "_", gpl_id)
  }
}
