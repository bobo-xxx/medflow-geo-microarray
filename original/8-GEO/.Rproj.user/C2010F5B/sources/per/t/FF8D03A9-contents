# ===================================================
## 0. 基础设置
# ===================================================
library(dplyr)
library(tidyr)
library(stringr)
library(oligo)
library(R.utils)
library(Biobase)
library(GEOquery)
library(tibble)

# ===================================================
## 1. 路径与参数
# ===================================================
data_dir     <- "/work/run/projects/bio-30/projects/8-GEO/data"  #数据所在文件夹
platform_dir <- file.path(data_dir, "platforms")   #平台所在文件夹
output_dir   <- file.path(data_dir, "output")      #输出文件
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gse_affy <- gse_keep_summary %>%
  filter(inferred_type == "Affymetrix_microarray") %>%
  pull(GSE)

# ===================================================
## 2. 批量循环处理 GSE
# ===================================================

for (gse_id in gse_affy) {
  
  out_gse <- file.path(output_dir, gse_id)
  dir.create(out_gse, recursive = TRUE, showWarnings = FALSE)
  
  expr_gene_file  <- file.path(out_gse, paste0("expr_gene_", gse_id, ".csv"))
  expr_probe_file <- file.path(out_gse, paste0("expr_probe_", gse_id, ".csv"))
  note_file       <- file.path(out_gse, "SKIPPED_REASON.txt")
  
  ## ---------- 已完成则跳过 ----------
  if (file.exists(expr_gene_file) && file.exists(expr_probe_file)) {
    message("✅ Already processed, skipping ", gse_id)
    next
  }
  
  message("=== Processing ", gse_id, " ===")
  
  ## ===================================================
  ## 3. 准备 RAW
  ## ===================================================
  
  suppl_dir <- file.path(data_dir, gse_id, "suppl")
  raw_dir   <- file.path(suppl_dir, "RAW")
  
  tar_files <- list.files(suppl_dir, pattern = "_RAW\\.tar$", full.names = TRUE)
  if (length(tar_files) == 0) {
    writeLines("No RAW.tar file found.", note_file)
    next
  }
  
  if (!dir.exists(raw_dir)) dir.create(raw_dir)
  if (length(list.files(raw_dir)) == 0) {
    untar(tar_files[1], exdir = raw_dir)
  }
  
  ## ===================================================
  ## 4. 查找 CEL
  ## ===================================================
  
  cel_files <- list.files(
    raw_dir,
    pattern = "\\.cel(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(cel_files) == 0) {
    writeLines("No CEL files found after extracting RAW.", note_file)
    next
  }
  
  ## ===================================================
  ## 5. 获取数据集的平台GPL
  ## ===================================================
  
  Sys.setenv(http_proxy="http://localhost:1086",
             https_proxy="http://localhost:1086")
  gse_obj <- getGEO(gse_id, GSEMatrix = FALSE)
  Sys.unsetenv(c("http_proxy", "https_proxy"))
  
  gpl_id <- names(GPLList(gse_obj))
  
  if (length(gpl_id) != 1) {
    writeLines(
      paste("Multiple or zero GPL detected:", paste(gpl_id, collapse = ", ")),
      note_file
    )
    next
  }
  
  gpl_id <- gpl_id[1]
  message("Using platform: ", gpl_id)

  ## ---------- 明确跳过 PrimeView ----------
  if (gpl_id %in% c("GPL15207")) {
    writeLines(
      paste0(
        "Platform ", gpl_id, " (PrimeView) requires pd.primeview, ",
        "which has been removed from Bioconductor and cannot be installed. ",
        "Dataset skipped."
      ),
      note_file
    )
    next
  }
  
  ## ===================================================
  ## 6. 过滤 CEL（真正的）
  ## ===================================================
  
  good_flags <- sapply(cel_files, function(f) {
    tryCatch({
      suppressMessages(read.celfiles(f))
      TRUE
    }, error = function(e) {
      FALSE
    })
  })
  
  good_cel_files <- cel_files[good_flags]
  bad_cel_files  <- basename(cel_files[!good_flags])
  
  if (length(bad_cel_files) > 0) {
    writeLines(
      bad_cel_files,
      file.path(out_gse, paste0("bad_cel_", gse_id, ".txt"))
    )
  }
  
  if (length(good_cel_files) == 0) {
    writeLines(
      "All CEL files failed to load (not due to platform annotation).",
      note_file
    )
    next
  }
  
  ## ===================================================
  ## 7. RMA
  ## ===================================================
  
  raw_data <- read.celfiles(good_cel_files)
  eset     <- rma(raw_data)
  
  expr_matrix <- exprs(eset)
  colnames(expr_matrix) <- basename(colnames(expr_matrix)) |>
    str_remove("\\.CEL(\\.gz)?$")
  
  ## ===================================================
  ## 8. 读取 GPL SOFT
  ## ===================================================
  
  gpl_soft_file <- file.path(platform_dir, gpl_id, paste0(gpl_id, ".soft.gz"))
  if (!file.exists(gpl_soft_file)) {
    writeLines(
      paste("GPL SOFT file not found:", gpl_soft_file),
      note_file
    )
    next
  }
  
  gpl <- getGEO(filename = gpl_soft_file)
  gpl_table <- Table(gpl)
  
  ## ===================================================
  ## 9. probe → gene
  ## ===================================================
  gene_col <- intersect(
    colnames(gpl_table),
    c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment")
  )[1]
  
  if (is.na(gene_col)) {
    writeLines(
      "No gene symbol column found in GPL annotation table.",
      note_file
    )
    next
  }
  
  probe2gene <- gpl_table %>%
    select(probe_id = ID, gene_symbol = all_of(gene_col)) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    mutate(gene_symbol = str_split(gene_symbol, " /// ")) %>%
    unnest(gene_symbol)
  
  ## ===================================================
  ## 10. probe → gene（mean）
  ## ===================================================
  expr_df <- as.data.frame(expr_matrix)
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
  
  ## ===================================================
  ## 11. 保存
  ## ===================================================
  
  write.csv(expr_matrix, expr_probe_file, row.names = TRUE)
  write.csv(expr_gene,   expr_gene_file,  row.names = TRUE)
  
  message("✅ Finished processing ", gse_id)
}
