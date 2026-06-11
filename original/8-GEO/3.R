library(dplyr)
library(stringr)
library(affy)
library(tibble)
data_dir <- "/work/run/projects/bio-30/projects/8-GEO/data/"

gse_affy <- gse_keep_summary %>%
  filter(inferred_type == "Affymetrix_microarray") %>%
  pull(GSE)

gse_affy[1:10]

gse_id <- "GSE125989"
suppl_dir <- file.path(data_dir, gse_id, "suppl")

tar_files <- list.files(
  suppl_dir,
  pattern = "_RAW\\.tar$",
  full.names = TRUE
)

stopifnot(length(tar_files) >= 1)


raw_dir <- file.path(suppl_dir, "RAW")

if (!dir.exists(raw_dir)) {
  dir.create(raw_dir)
  untar(tar_files[1], exdir = raw_dir)
}

cel_files <- list.files(
  raw_dir,
  pattern = "\\.CEL(\\.gz)?$",
  full.names = TRUE,
  recursive = TRUE
)

length(cel_files)
stopifnot(length(cel_files) > 0)

raw_affy <- ReadAffy(filenames = cel_files)

eset <- rma(raw_affy)

expr_probe <- exprs(eset)
dim(expr_probe)



#测试
library(affy)

test_one_cel <- function(cel_file) {
  tryCatch(
    {
      ReadAffy(filenames = cel_file)
      TRUE
    },
    error = function(e) {
      message("❌ Bad CEL: ", basename(cel_file))
      FALSE
    }
  )
}

good_flags <- sapply(cel_files, test_one_cel)

table(good_flags)


good_flags



good_cel_files <- cel_files[good_flags]

length(good_cel_files)
stopifnot(length(good_cel_files) >= 3)  # 防止样本太少

raw_affy <- ReadAffy(filenames = good_cel_files)
eset <- rma(raw_affy)

expr_probe <- exprs(eset)
dim(expr_probe)


bad_cel <- cel_files[!good_flags]

bad_log <- tibble(
  GSE = gse_id,
  bad_CEL = basename(bad_cel)
)

bad_log




