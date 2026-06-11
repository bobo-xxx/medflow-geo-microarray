# ===================================================
## 0. 基础环境
# ===================================================
suppressPackageStartupMessages({
  library(limma)
  library(GEOquery)
  library(dplyr)
  library(R.utils)
})

options(stringsAsFactors = FALSE)

# ===================================================
## 1. 路径设置
# ===================================================
data_dir   <- "/work/run/projects/bio-30/projects/8-GEO/data"
output_dir <- file.path(data_dir, "outputT")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===================================================
## 2. 待处理 GSE
# ===================================================
gse_gpr <- gse_keep_summary %>%
  filter(inferred_type == "Agilent_microarray") %>%  
  pull(GSE)

# ===================================================
## 3. 自动检测文件类型函数
# ===================================================
detect_gpr_type <- function(f) {
  x <- try(readLines(f, n = 1), silent = TRUE)
  if (inherits(x, "try-error")) return("unknown")
  if (grepl("^ATF", x)) return("genepix")     # 双通道 GenePix
  if (grepl("^TYPE", x)) return("agilent")    # Agilent Feature Extraction
  return("unknown")
}

# ===================================================
## 4. 主循环
# ===================================================
for (i in seq_along(gse_gpr)) {
  
  gse_id <- gse_gpr[i]
  message(sprintf("\n=== [%d/%d] Processing %s ===",
                  i, length(gse_gpr), gse_id))
  
  out_gse <- file.path(output_dir, gse_id)
  dir.create(out_gse, showWarnings = FALSE, recursive = TRUE)
  
  expr_file <- file.path(out_gse, paste0("expr_probe_", gse_id, ".csv"))
  note_file <- file.path(out_gse, "SKIPPED_REASON.txt")
  
  if (file.exists(expr_file)) {
    message("✅ Already processed, skipping ", gse_id)
    next
  }
  
  # ---------------------------------------------------
  # 4.1 定位 RAW 目录
  # ---------------------------------------------------
  suppl_dir <- file.path(data_dir, gse_id, "suppl")
  raw_dir   <- file.path(suppl_dir, "RAW")
  dir.create(raw_dir, showWarnings = FALSE)
  
  # ---------------------------------------------------
  # 4.2 解压 tar / zip
  # ---------------------------------------------------
  archives <- list.files(
    suppl_dir,
    pattern = "\\.(tar|zip)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(archives) > 0 && length(list.files(raw_dir)) == 0) {
    message("Extracting archive(s)...")
    for (a in archives) {
      tryCatch({
        if (grepl("\\.tar$", a, ignore.case = TRUE)) {
          untar(a, exdir = raw_dir)
        } else {
          unzip(a, exdir = raw_dir)
        }
      }, error = function(e) {
        message("❌ Archive extract failed: ", basename(a))
      })
    }
  }
  
  # ---------------------------------------------------
  # 4.3 解压 .gpr.gz
  # ---------------------------------------------------
  gz_files <- list.files(
    raw_dir,
    pattern = "\\.gpr\\.gz$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(gz_files) > 0) {
    message("Decompressing .gpr.gz files...")
    for (gz in gz_files) {
      tryCatch({
        gunzip(gz, remove = FALSE, overwrite = TRUE)
      }, error = function(e) {})
    }
  }
  
  # ---------------------------------------------------
  # 4.4 收集 GPR 文件
  # ---------------------------------------------------
  gpr_files <- list.files(
    raw_dir,
    pattern = "\\.gpr$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(gpr_files) == 0) {
    writeLines("❌ No GPR files found.", note_file)
    next
  }
  
  message("GPR files found: ", length(gpr_files))
  
  # ---------------------------------------------------
  # 4.5 检测文件类型
  # ---------------------------------------------------
  gpr_type <- detect_gpr_type(gpr_files[1])
  message("Detected GPR type: ", gpr_type)
  
  if (gpr_type == "unknown") {
    writeLines("❌ Unknown GPR file type.", note_file)
    next
  }
  
  # ---------------------------------------------------
  # 4.6 读取并处理
  # ---------------------------------------------------
  RG <- tryCatch({
    if (gpr_type == "genepix") {
      message("Reading GenePix (two-color) files...")
      read.maimages(gpr_files, source = "genepix")
    } else if (gpr_type == "agilent") {
      message("Reading Agilent FE (single-color) files...")
      read.maimages(gpr_files, source = "agilent", green.only = TRUE)
    }
  }, error = function(e) {
    writeLines(paste("❌ read.maimages failed:", e$message), note_file)
    return(NULL)
  })
  
  if (is.null(RG)) next
  
  # ---------------------------------------------------
  # 4.7 背景校正 + 归一化
  # ---------------------------------------------------
  message("Background correction + normalization...")
  
  if (gpr_type == "genepix") {
    RG.bg <- backgroundCorrect(RG, method = "normexp", offset = 16)
    MA <- normalizeWithinArrays(RG.bg, method = "loess")
    MA <- normalizeBetweenArrays(MA, method = "quantile")
    expr <- as.data.frame(MA$A)
  } else if (gpr_type == "agilent") {
    RG.bg <- backgroundCorrect(RG, method = "normexp", offset = 16)
    RG.qn <- normalizeBetweenArrays(RG.bg, method = "quantile")
    expr <- as.data.frame(RG.qn$E)
  }
  
  colnames(expr) <- basename(colnames(expr))
  write.csv(expr, expr_file, row.names = TRUE)
  
  message("✅ Finished ", gse_id)
}
