# annotate.R — Probe annotation and gene aggregation for geo-microarray-processing
#
# Implements the 5-tier gene annotation fallback:
#   1. fData() direct column (Gene Symbol / GENE_SYMBOL / Symbol)
#   2. fData() gene_assignment column (parse "ACC // SYMBOL // desc")
#   3. GPL annotation table (GEOquery::Table(getGEO(GPL)))
#   4. Bioconductor annotation DB (platform-specific .db package)
#   5. Probe IDs as gene symbols (last resort with warning)

library(GEOquery)
library(dplyr)
library(tidyr)
library(stringr)

#' Extract gene symbol from Affymetrix gene_assignment format
#'
#' Parses "Accession // GeneSymbol // Description // ..." format.
#'
#' @param x Character vector of gene_assignment strings
#' @return Character vector of gene symbols (NA if unparseable)
extract_gene_from_assignment <- function(x) {
  parts <- strsplit(as.character(x), " // ", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) >= 2) trimws(p[2]) else NA_character_
  }, "")
}

#' Get GPL annotation table with probe-to-gene mapping
#'
#' Downloads GPL annotation from GEO and extracts ID -> gene symbol mapping.
#' Searches common column names: Gene Symbol, GENE_SYMBOL, Symbol, etc.
#'
#' @param gpl_id GPL identifier (e.g., "GPL570")
#' @param destdir Optional directory for caching
#' @return data.frame with probe_id and gene_symbol columns, or NULL
get_gpl_annotation <- function(gpl_id, destdir = NULL, probe_ids = NULL) {
  if (is.null(gpl_id) || gpl_id == "") return(NULL)

  gpl_id <- toupper(gpl_id)
  if (!grepl("^GPL", gpl_id)) gpl_id <- paste0("GPL", gpl_id)

  # Check cache
  if (!is.null(destdir) && dir.exists(destdir)) {
    cache_file <- file.path(destdir, paste0(gpl_id, "_annotation.rds"))
    if (file.exists(cache_file)) {
      message("Loading cached GPL annotation: ", gpl_id)
      return(readRDS(cache_file))
    }
  }

  # Download
  message("Downloading GPL annotation: ", gpl_id)
  gpl <- tryCatch(getGEO(gpl_id), error = function(e) {
    message("Failed to download GPL ", gpl_id, ": ", e$message)
    NULL
  })
  if (is.null(gpl)) return(NULL)

  gpl_table <- tryCatch(Table(gpl), error = function(e) {
    message("Failed to extract GPL table: ", e$message)
    NULL
  })
  if (is.null(gpl_table)) return(NULL)

  col_names <- colnames(gpl_table)

  # Find ID column
  id_col <- which(toupper(col_names) == "ID")[1]
  if (is.na(id_col)) {
    message("Cannot find ID column in GPL ", gpl_id)
    return(NULL)
  }

  # Find Gene Symbol column - ordered preference
  preferred_cols <- c("GENE_SYMBOL", "SYMBOL", "GENE NAME",
    "GENE SYMBOL (EXTERNAL)", "GENESYMBOL", "GENE_NAME", "GENE SYMBOL")
  gene_col <- NA
  for (pref_col in preferred_cols) {
    idx <- which(toupper(col_names) == toupper(pref_col))[1]
    if (!is.na(idx)) { gene_col <- idx; break }
  }

  if (is.na(gene_col)) {
    message("No gene symbol column in GPL ", gpl_id,
      ". Available: ", paste(col_names[1:min(10, length(col_names))], collapse = ", "))
    # Try supplementary annotation file (Agilent arrays often have richer GPL files)
    suppl_result <- try_get_gpl_suppl(gpl, gpl_table, id_col, destdir)
    if (!is.null(suppl_result)) {
      result <- suppl_result
    } else {
      result <- data.frame(probe_id = as.character(gpl_table[, id_col]),
                           gene_symbol = NA_character_, stringsAsFactors = FALSE)
    }
  } else {
    result <- data.frame(
      probe_id    = as.character(gpl_table[, id_col]),
      gene_symbol = as.character(gpl_table[, gene_col]),
      stringsAsFactors = FALSE
    )
    result <- result[result$gene_symbol != "" & !is.na(result$gene_symbol), ]
  }

  # Cache
  if (!is.null(destdir) && dir.exists(destdir)) {
    saveRDS(result, file.path(destdir, paste0(gpl_id, "_annotation.rds")))
  }

  result
}

#' Aggregate probe-level expression to gene-level
#'
#' Joins expression matrix with probe-to-gene mapping and aggregates
#' multiple probes per gene using mean.
#'
#' @param expr_matrix Numeric matrix (probes x samples) with probe rownames
#' @param gpl_table data.frame with probe_id and gene_symbol columns
#' @return data.frame with gene_symbol rows and sample columns, or NULL
aggregate_probe_to_gene <- function(expr_matrix, gpl_table) {
  if (is.null(expr_matrix) || nrow(expr_matrix) == 0) return(NULL)
  if (is.null(gpl_table) || nrow(gpl_table) == 0) return(NULL)
  if (!all(c("probe_id", "gene_symbol") %in% colnames(gpl_table))) return(NULL)

  # Convert matrix to data.frame with probe_id column
  expr_df <- as.data.frame(expr_matrix)
  expr_df$probe_id <- rownames(expr_matrix)

  # Handle " /// " separated gene symbols (expand one-to-many)
  probe2gene <- gpl_table %>%
    select(probe_id, gene_symbol) %>%
    mutate(gene_symbol = str_split(gene_symbol, " /// ")) %>%
    unnest(gene_symbol) %>%
    mutate(gene_symbol = trimws(gene_symbol)) %>%
    filter(gene_symbol != "" & !is.na(gene_symbol))

  # Melt, join, aggregate by mean
  expr_long <- expr_df %>%
    pivot_longer(cols = -probe_id, names_to = "sample", values_to = "expr")

  expr_gene <- expr_long %>%
    inner_join(probe2gene, by = "probe_id", relationship = "many-to-many") %>%
    group_by(gene_symbol, sample) %>%
    summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
    mutate(expr = ifelse(is.nan(expr), NA_real_, expr)) %>%
    pivot_wider(names_from = sample, values_from = expr) %>%
    as.data.frame()

  rownames(expr_gene) <- expr_gene$gene_symbol
  expr_gene$gene_symbol <- NULL

  expr_gene
}

#' Clean gene symbols from Bioconductor annotation DB output
#'
#' Strips entg| prefix that some BioC annotation packages return.
#'
#' @param symbols Character vector of gene symbols
#' @return Cleaned character vector
clean_gene_symbols <- function(symbols) {
  sub("^entg\\|", "", symbols)
}

#' Try downloading and parsing a GPL supplementary annotation file
#'
#' Checks Meta(gpl)$supplementary_file for a richer annotation table.
#' Downloads if the file is under 1GB.
#'
#' @param gpl GPL object from GEOquery::getGEO
#' @param gpl_table Current (limited) GPL table
#' @param id_col Index of the ID column
#' @param destdir Optional cache directory
#' @return data.frame with probe_id and gene_symbol, or NULL
try_get_gpl_suppl <- function(gpl, gpl_table, id_col, destdir = NULL) {
  MAX_SIZE <- 1e9  # 1 GB limit

  suppl_url <- tryCatch(Meta(gpl)$supplementary_file, error = function(e) NULL)
  if (is.null(suppl_url) || nchar(suppl_url) == 0) return(NULL)

  message("Trying GPL supplementary file: ", basename(suppl_url))

  # Check cache
  cache_file <- if (!is.null(destdir) && dir.exists(destdir)) {
    file.path(destdir, paste0(basename(suppl_url), ".rds"))
  } else NULL

  if (!is.null(cache_file) && file.exists(cache_file)) {
    message("  Loading cached suppl annotation")
    return(readRDS(cache_file))
  }

  # Check size before download
  tmp <- tempfile(fileext = ".gz")
  h <- tryCatch(curl::curl_fetch_disk(suppl_url, tmp), error = function(e) NULL)
  if (is.null(h)) {
    # Fallback to base R download
    tryCatch(download.file(suppl_url, tmp, method = "auto", quiet = TRUE),
             error = function(e) { message("  Download failed"); return(NULL) })
  }

  file_size <- file.info(tmp)$size
  if (is.na(file_size) || file_size == 0) { unlink(tmp); return(NULL) }
  if (file_size > MAX_SIZE) {
    message(sprintf("  Suppl file too large (%.0f MB > 1 GB limit), skipping",
                    file_size / 1e6))
    unlink(tmp); return(NULL)
  }

  # Parse
  result <- tryCatch(parse_gpl_suppl_soft(gzfile(tmp)), error = function(e) {
    message("  Parse failed: ", e$message); NULL
  })
  unlink(tmp)

  if (is.null(result) || nrow(result) == 0) return(NULL)

  message(sprintf("  Extracted %d probe→gene mappings from suppl file", nrow(result)))

  # Cache
  if (!is.null(cache_file)) saveRDS(result, cache_file)

  result
}

#' Parse a GPL supplementary SOFT file for gene symbol annotations
#'
#' Reads the tab-delimited data table from a GPL SOFT file and extracts
#' probe ID → gene symbol mappings. Used when Table(gpl) lacks gene symbols.
#'
#' @param lines Character vector of SOFT file lines (or path to file)
#' @return data.frame with probe_id and gene_symbol, or NULL
parse_gpl_suppl_soft <- function(lines_or_path) {
  # If given a file path, read it
  if (length(lines_or_path) == 1 && file.exists(lines_or_path)) {
    lines_or_path <- readLines(lines_or_path, warn = FALSE)
  }

  # Find the data table section
  table_begin <- grep("^!Platform_table_begin", lines_or_path, ignore.case = TRUE)
  if (length(table_begin) == 0) return(NULL)

  header_line <- table_begin[1] + 1
  if (header_line > length(lines_or_path)) return(NULL)

  # Parse header
  headers <- strsplit(lines_or_path[header_line], "\t")[[1]]

  # Check for gene symbol column
  symbol_col <- grep("GENE_SYMBOL|gene_symbol|Gene Symbol|SYMBOL", headers, ignore.case = TRUE)[1]
  id_col <- grep("^ID$", headers, ignore.case = FALSE)[1]
  if (is.na(symbol_col) || is.na(id_col)) return(NULL)

  # Parse data rows (skip header, stop at table_end or !)
  data_start <- header_line + 1
  data_end <- length(lines_or_path)
  table_end <- grep("^!Platform_table_end", lines_or_path, ignore.case = TRUE)
  if (length(table_end) > 0) data_end <- table_end[1] - 1

  data_lines <- lines_or_path[data_start:data_end]
  data_lines <- data_lines[!grepl("^[!]|^$", data_lines)]
  if (length(data_lines) == 0) return(NULL)

  # Parse each data line as tab-separated
  parsed <- strsplit(data_lines, "\t")
  # Filter to rows with enough columns
  ncols_expected <- length(headers)
  valid <- vapply(parsed, length, 0L) >= ncols_expected
  parsed <- parsed[valid]

  probe_ids <- vapply(parsed, `[`, "", id_col)
  symbols  <- vapply(parsed, `[`, "", symbol_col)

  # Build result, exclude --- entries
  result <- data.frame(
    probe_id    = probe_ids,
    gene_symbol = symbols,
    stringsAsFactors = FALSE
  )
  result <- result[!grepl("^---$", result$gene_symbol) & result$gene_symbol != "", ]
  if (nrow(result) == 0) return(NULL)
  result
}

#' Map GPL platform ID to Bioconductor annotation package name
#'
#' Returns the .db package name for common Affymetrix platforms.
#' These packages provide probe ID → gene symbol mappings via AnnotationDbi.
#'
#' @param gpl_id GPL identifier (e.g., "GPL570")
#' @return Character package name, or NULL if no mapping exists
gpl_to_bioc_package <- function(gpl_id) {
  mapping <- c(
    "GPL570"   = "hgu133plus2.db",
    "GPL96"    = "hgu133a.db",
    "GPL571"   = "hgu133a2.db",
    "GPL16686" = "hugene20sttranscriptcluster.db",
    "GPL17586" = "hta20transcriptcluster.db",
    "GPL6244"  = "hugene10sttranscriptcluster.db",
    "GPL1261"  = "mouse4302.db",
    "GPL339"   = "mouse430a2.db"
  )
  gpl_id <- toupper(gpl_id)
  if (gpl_id %in% names(mapping)) unname(mapping[gpl_id]) else NULL
}

#' Annotate probes using Bioconductor annotation database
#'
#' Maps probe IDs to gene symbols using a platform-specific .db package.
#'
#' @param probe_ids Character vector of probe IDs
#' @param gpl_id GPL platform identifier
#' @return data.frame with probe_id and gene_symbol columns, or NULL
annotate_with_bioc_db <- function(probe_ids, gpl_id) {
  pkg <- gpl_to_bioc_package(gpl_id)
  if (is.null(pkg)) return(NULL)
  if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)

  suppressMessages({
    db <- get(pkg, envir = asNamespace(pkg))
    mapping <- tryCatch(
      AnnotationDbi::select(db, keys = probe_ids,
                            columns = "SYMBOL", keytype = "PROBEID"),
      error = function(e) { message("  BioC DB query failed: ", e$message); NULL }
    )
  })

  if (is.null(mapping) || nrow(mapping) == 0) return(NULL)

  result <- mapping[!is.na(mapping$SYMBOL) & mapping$SYMBOL != "", ]
  if (nrow(result) == 0) return(NULL)

  data.frame(
    probe_id    = as.character(result$PROBEID),
    gene_symbol = clean_gene_symbols(as.character(result$SYMBOL)),
    stringsAsFactors = FALSE
  )
}
