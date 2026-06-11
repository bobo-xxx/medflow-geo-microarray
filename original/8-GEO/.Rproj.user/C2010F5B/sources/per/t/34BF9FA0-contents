## ===============================
## 0. 基础设置
## ===============================
rm(list = ls())
gc()

library(dplyr)
library(stringr)
library(readr)
library(purrr)

data_dir <- "/work/run/projects/bio-30/projects/8-GEO/data"

## ===============================
## 1. 获取所有 GSE 目录
## ===============================
gse_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
gse_dirs <- gse_dirs[grepl("GSE", basename(gse_dirs))]

length(gse_dirs)

## ===============================
## 2. 定义：单个 GSE 的判断函数
## ===============================
inspect_one_gse <- function(gse_path) {
  
  gse_id <- basename(gse_path)
  suppl_dir <- file.path(gse_path, "suppl")
  
  ## ---- 初始化结果 ----
  res <- tibble(
    GSE = gse_id,
    suppl_exists = FALSE,
    suppl_empty = NA,
    has_filelist = FALSE,
    has_BPM = FALSE,
    decision = NA_character_,
    reason = NA_character_
  )
  
  ## ---- Step 1: 检查 suppl 是否存在 ----
  if (!dir.exists(suppl_dir)) {
    res$suppl_exists <- FALSE
    res$suppl_empty  <- TRUE
    res$decision <- "skip"
    res$reason   <- "no_suppl_directory"
    return(res)
  }
  
  res$suppl_exists <- TRUE
  
  suppl_files <- list.files(suppl_dir, full.names = TRUE)
  
  if (length(suppl_files) == 0) {
    res$suppl_empty <- TRUE
    res$decision <- "skip"
    res$reason   <- "suppl_empty"
    return(res)
  }
  
  res$suppl_empty <- FALSE
  
  ## ---- Step 2: 检查 filelist.txt ----
  filelist_path <- file.path(suppl_dir, "filelist.txt")
  
  if (!file.exists(filelist_path)) {
    res$has_filelist <- FALSE
    res$decision <- "undetermined"
    res$reason   <- "no_filelist"
    return(res)
  }
  
  res$has_filelist <- TRUE
  
  ## ---- Step 3: 读取 filelist，判断 BPM ----
  fl <- suppressWarnings(
    read.table(filelist_path, header = FALSE, stringsAsFactors = FALSE)
  )
  
  if (ncol(fl) >= 2) {
    filenames <- fl[[2]]
    res$has_BPM <- any(grepl("\\.bpm(\\.gz)?$", filenames, ignore.case = TRUE))
  }
  
  if (isTRUE(res$has_BPM)) {
    res$decision <- "skip"
    res$reason   <- "methylation_BPM_present"
    return(res)
  }
  
  ## ---- Step 4: 其他情况（暂不跳过） ----
  res$decision <- "keep"
  res$reason   <- "no_BPM_suppl_nonempty"
  
  return(res)
}

## ===============================
## 3. 批量运行所有 GSE
## ===============================
gse_judgement <- map_dfr(gse_dirs, inspect_one_gse)

## ===============================
## 4. 排序 & 查看结果
## ===============================
gse_judgement <- gse_judgement %>%
  arrange(decision, GSE)

print(head(gse_judgement, 10))

## ===============================
## 5. 保存结果
## ===============================
write.csv(
  gse_judgement,
  file = "GSE_suppl_judgement_summary.csv",
  row.names = FALSE
)

