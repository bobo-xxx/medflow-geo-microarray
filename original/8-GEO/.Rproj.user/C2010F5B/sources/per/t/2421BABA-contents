#Step 0：准备环境 & 路径
library(dplyr)
library(stringr)
library(purrr)

data_dir <- "/work/run/projects/bio-30/projects/8-GEO/data"

#Step 1：只取 decision == "keep" 的 GSE
gse_keep <- gse_judgement %>%
  filter(decision == "keep") %>%
  pull(GSE)

#Step 2：构造 keep GSE 对应的 filelist.txt 路径
filelist_paths <- file.path(
  data_dir,
  gse_keep,
  "suppl",
  "filelist.txt"
)

filelist_paths <- filelist_paths[file.exists(filelist_paths)]
length(filelist_paths)


#Step 3：稳健解析单个 filelist.txt（重点）
parse_one_filelist <- function(filelist_path) {
  
  gse_id <- basename(dirname(dirname(filelist_path)))
  
  df <- read.table(
    filelist_path,
    header = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE,
    quote = ""
  )
  
  # GEO 的 filelist.txt 至少 6 列，防御式处理
  if (ncol(df) < 6) return(NULL)
  
  df <- df[, 1:6]
  colnames(df) <- c(
    "entry_type",
    "filename",
    "date",
    "time",
    "size",
    "suffix"
  )
  
  df$GSE <- gse_id
  df
}

#Step 4：合并所有 keep GSE 的 filelist 信息
all_files_df <- map_dfr(filelist_paths, parse_one_filelist)

head(all_files_df)


#Step 5：逐 GSE 判断包含哪些关键文件类型
library(stringr)
library(dplyr)

gse_file_summary <- all_files_df %>%
  group_by(GSE) %>%
  summarise(
    has_RAW_tar = any(str_detect(filename, regex("_RAW\\.tar", ignore_case = TRUE))),
    has_CEL     = any(str_detect(filename, regex("\\.CEL(\\.gz)?$", ignore_case = TRUE))),
    has_GPR     = any(str_detect(filename, regex("\\.GPR(\\.gz)?$", ignore_case = TRUE))),
    has_IDAT    = any(str_detect(filename, regex("\\.idat(\\.gz)?$", ignore_case = TRUE))),
    has_BPM     = any(str_detect(filename, regex("\\.bpm", ignore_case = TRUE))),
    has_TXT     = any(str_detect(filename, regex("\\.txt(\\.gz)?$", ignore_case = TRUE))),
    has_CSV     = any(str_detect(filename, regex("\\.csv(\\.gz)?$", ignore_case = TRUE))),
    has_XLSX    = any(str_detect(filename, regex("\\.xlsx(\\.gz)?$", ignore_case = TRUE))),
    has_matrix  = any(str_detect(filename, regex("series_matrix", ignore_case = TRUE))),
    file_types  = paste(sort(unique(suffix)), collapse = ";"),
    .groups = "drop"
  )


#Step 6：推断「表达数据来源类型」（只做描述，不做清洗决策）
gse_file_summary <- gse_file_summary %>%
  mutate(
    inferred_type = case_when(
      has_CEL ~ "Affymetrix_microarray",
      has_GPR ~ "Agilent_microarray",
      has_IDAT & has_BPM ~ "Methylation_array",
      has_IDAT & !has_BPM ~ "Illumina_expression_array",
      has_TXT & !has_IDAT & !has_CEL & !has_GPR ~ "TXT_array",
      has_matrix ~ "Matrix_only",
      TRUE ~ "Unknown"
    )
  )

#Step 7：和原始判断表合并（形成“总控表”）
gse_keep_summary <- gse_judgement %>%
  filter(decision == "keep") %>%
  left_join(gse_file_summary, by = "GSE")

#Step 8：导出（这一张表非常关键）
write.csv(
  gse_keep_summary,
  file = "GSE_keep_file_composition_summary.csv",
  row.names = FALSE
)


#Step 0：锁定 Unknown 的 GSE
unknown_gse <- gse_keep_summary %>%
  filter(inferred_type == "Unknown", has_RAW_tar) %>%
  pull(GSE)

length(unknown_gse)

#Step 1：定义 RAW.tar 解压 + 扫描函数
inspect_raw_tar <- function(gse_id, data_dir, tmp_base = tempdir()) {
  
  tar_path <- file.path(
    data_dir,
    gse_id,
    "suppl",
    paste0(gse_id, "_RAW.tar")
  )
  
  if (!file.exists(tar_path)) {
    return(tibble(
      GSE = gse_id,
      raw_has_CEL = FALSE,
      raw_has_GPR = FALSE,
      raw_has_IDAT = FALSE,
      raw_has_BPM = FALSE,
      raw_file_types = NA,
      raw_inferred_type = "No_RAW_tar"
    ))
  }
  
  tmp_dir <- file.path(tmp_base, gse_id)
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  
  untar(tar_path, exdir = tmp_dir)
  
  files <- list.files(tmp_dir, recursive = TRUE, full.names = FALSE)
  
  res <- tibble(
    GSE = gse_id,
    raw_has_CEL  = any(str_detect(files, regex("\\.CEL$", ignore_case = TRUE))),
    raw_has_GPR  = any(str_detect(files, regex("\\.GPR$", ignore_case = TRUE))),
    raw_has_IDAT = any(str_detect(files, regex("\\.idat$", ignore_case = TRUE))),
    raw_has_BPM  = any(str_detect(files, regex("\\.bpm", ignore_case = TRUE))),
    raw_file_types = paste(sort(unique(tools::file_ext(files))), collapse = ";")
  )
  
  res <- res %>%
    mutate(
      raw_inferred_type = case_when(
        raw_has_CEL ~ "Affymetrix_microarray",
        raw_has_GPR ~ "Agilent_microarray",
        raw_has_IDAT & raw_has_BPM ~ "Methylation_array",
        raw_has_IDAT & !raw_has_BPM ~ "Illumina_expression_array",
        TRUE ~ "Still_Unknown"
      )
    )
  
  unlink(tmp_dir, recursive = TRUE)
  res
}

#Step 2：批量扫描 Unknown GSE 的 RAW.tar
raw_scan_results <- purrr::map_dfr(
  unknown_gse,
  inspect_raw_tar,
  data_dir = data_dir
)

head(raw_scan_results)

#Step 3：回填修正 gse_keep_summary（关键一步）
gse_keep_summary_updated <- gse_keep_summary %>%
  left_join(raw_scan_results, by = "GSE") %>%
  mutate(
    inferred_type_final = if_else(
      inferred_type == "Unknown" & !is.na(raw_inferred_type),
      raw_inferred_type,
      inferred_type
    )
  )

#Step 4：导出最终「可处理判定表」
write.csv(
  gse_keep_summary_updated,
  file = "GSE_keep_file_composition_summary_final.csv",
  row.names = FALSE
)


#setdiff(table_names_clean, all_subdirs_clean)

#[1] "GSE36002" "GSE37816" "GSE38860" "GSE43256" "GSE59444" "GSE92324"






