# ===================================================
## 0. 基础环境
# ===================================================
suppressPackageStartupMessages({
  library(limma)
  library(GEOquery)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(R.utils)
})

options(stringsAsFactors = FALSE)

# ===================================================
## 1. 路径设置
# ===================================================
data_dir     <- "/work/run/projects/bio-30/projects/8-GEO/data"
output_dir   <- file.path(data_dir, "outputT")
platform_dir <- file.path(data_dir, "platforms")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===================================================
## 2. 待处理 GSE（Agilent GPR）
# ===================================================
gse_gpr <- gse_keep_summary %>%
  filter(inferred_type == "Agilent_microarray") %>%
  pull(GSE)

# ===================================================
## 3. 自动检测 GPR 类型
# ===================================================
detect_gpr_type <- function(f) {
  x <- try(readLines(f, n = 1), silent = TRUE)
  if (inherits(x, "try-error")) return("unknown")
  if (grepl("^ATF", x))  return("genepix")
  if (grepl("^TYPE", x)) return("agilent")
  return("unknown")
}

# ===================================================
## 4. 主循环
# ===================================================
for (gse_id in gse_gpr) {
  
  message("========== ", gse_id, " ==========")
  
  out_gse <- file.path(output_dir, gse_id)
  dir.create(out_gse, showWarnings = FALSE)
  
  note_file <- file.path(out_gse, paste0(gse_id, "_note.txt"))
  
  # -------------------------------------------------
  # 4.1 SUPPL / RAW
  # -------------------------------------------------
  suppl_dir <- file.path(data_dir, gse_id, "suppl")
  raw_dir   <- file.path(suppl_dir, "RAW")
  dir.create(raw_dir, showWarnings = FALSE)
  
  # -------------------------------------------------
  # 4.2 解压 tar / zip
  # -------------------------------------------------
  archives <- list.files(
    suppl_dir,
    pattern = "\\.(tar|zip)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(archives) > 0 && length(list.files(raw_dir)) == 0) {
    for (a in archives) {
      tryCatch({
        if (grepl("\\.tar$", a, ignore.case = TRUE)) {
          untar(a, exdir = raw_dir)
        } else {
          unzip(a, exdir = raw_dir)
        }
      }, error = function(e) {})
    }
  }
  
  # -------------------------------------------------
  # 4.3 解压 .gpr.gz
  # -------------------------------------------------
  gz_files <- list.files(
    raw_dir,
    pattern = "\\.gpr\\.gz$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  for (gz in gz_files) {
    tryCatch(
      gunzip(gz, remove = FALSE, overwrite = TRUE),
      error = function(e) {}
    )
  }
  
  # -------------------------------------------------
  # 4.4 收集 GPR
  # -------------------------------------------------
  gpr_files <- list.files(
    raw_dir,
    pattern = "\\.gpr$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(gpr_files) == 0) {
    writeLines("No GPR files found.", note_file)
    next
  }
  
  # -------------------------------------------------
  # 4.5 检测 GPR 类型
  # -------------------------------------------------
  gpr_type <- detect_gpr_type(gpr_files[1])
  message("Detected type: ", gpr_type)
  
  if (gpr_type == "unknown") {
    writeLines("Unknown GPR type.", note_file)
    next
  }
  
  # -------------------------------------------------
  # 4.6 读取 GPR
  # -------------------------------------------------
  RG <- tryCatch({
    if (gpr_type == "genepix") {
      read.maimages(gpr_files, source = "genepix")
    } else {
      read.maimages(gpr_files, source = "agilent", green.only = TRUE)
    }
  }, error = function(e) {
    writeLines(paste("read.maimages failed:", e$message), note_file)
    return(NULL)
  })
  
  if (is.null(RG)) next
  
  # -------------------------------------------------
  # 4.7 背景校正 + 归一化
  # -------------------------------------------------
  if (gpr_type == "genepix") {
    RG.bg <- backgroundCorrect(RG, method = "normexp", offset = 16)
    MA    <- normalizeWithinArrays(RG.bg, method = "loess")
    MA    <- normalizeBetweenArrays(MA, method = "quantile")
    expr  <- as.data.frame(MA$A)
  } else {
    RG.bg <- backgroundCorrect(RG, method = "normexp", offset = 16)
    RG.qn <- normalizeBetweenArrays(RG.bg, method = "quantile")
    expr  <- as.data.frame(RG.qn$E)
  }
  
  colnames(expr) <- basename(colnames(expr))
  
  # -------------------------------------------------
  # 4.8 GSE → GPL
  # -------------------------------------------------
  gse_obj <- getGEO(gse_id, GSEMatrix = FALSE)
  gpl_id  <- names(GPLList(gse_obj))
  
  if (length(gpl_id) != 1) {
    writeLines(paste("Invalid GPL:", paste(gpl_id, collapse = ",")), note_file)
    next
  }
  
  gpl_id <- gpl_id[1]
  message("Using GPL: ", gpl_id)
  
  # -------------------------------------------------
  # 4.9 读取 GPL SOFT
  # -------------------------------------------------
  gpl_soft <- file.path(platform_dir, gpl_id, paste0(gpl_id, ".soft.gz"))
  if (!file.exists(gpl_soft)) {
    writeLines("GPL SOFT not found.", note_file)
    next
  }
  
  gpl <- getGEO(filename = gpl_soft)
  gpl_table <- Table(gpl)
  
  # -------------------------------------------------
  # 4.10 自动确定 probe ID 列
  # -------------------------------------------------
  probe_col <- intersect(
    colnames(RG$genes),
    c("ID", "ProbeName", "Name", "ID_REF")
  )[1]
  
  if (is.na(probe_col)) {
    writeLines(
      paste("No probe ID column in RG$genes. Columns:",
            paste(colnames(RG$genes), collapse = ", ")),
      note_file
    )
    next
  }
  
  expr$probe_id <- RG$genes[[probe_col]]
  
  # -------------------------------------------------
  # 4.11 probe → gene
  # -------------------------------------------------
  gene_col <- intersect(
    colnames(gpl_table),
    c("Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_assignment")
  )[1]
  
  if (is.na(gene_col)) {
    writeLines("No gene symbol column in GPL.", note_file)
    next
  }
  
  probe2gene <- gpl_table %>%
    select(probe_id = ID, gene_symbol = all_of(gene_col)) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    mutate(gene_symbol = str_split(gene_symbol, " /// ")) %>%
    unnest(gene_symbol) %>%
    filter(probe_id %in% expr$probe_id)
  
  # -------------------------------------------------
  # 4.12 probe → gene（mean）
  # -------------------------------------------------
  expr_gene <- expr %>%
    pivot_longer(-probe_id, names_to = "sample", values_to = "expr") %>%
    inner_join(probe2gene, by = "probe_id") %>%
    group_by(gene_symbol, sample) %>%
    summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = expr) %>%
    as.data.frame()
  
  rownames(expr_gene) <- expr_gene$gene_symbol
  expr_gene$gene_symbol <- NULL
  
  # -------------------------------------------------
  # 4.13 输出
  # -------------------------------------------------
  write.csv(expr,
            file = file.path(out_gse, paste0(gse_id, "_probe.csv")),
            row.names = TRUE)
  
  write.csv(expr_gene,
            file = file.path(out_gse, paste0(gse_id, "_gene.csv")),
            row.names = TRUE)
  
  message("✅ Finished ", gse_id)
}
