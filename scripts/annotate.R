# annotate.R — Probe annotation and gene aggregation for geo-microarray-processing
#
# Implements the 5-tier gene annotation fallback:
#   1. fData() direct column (Gene Symbol / GENE_SYMBOL / Symbol)
#   2. fData() gene_assignment column (parse "ACC // SYMBOL // desc")
#   3. GPL annotation table (GEOquery::Table(getGEO(GPL)))
#   4. AnnoProbe pipe alignment (probe FASTA -> genome -> GENCODE)
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
get_gpl_annotation <- function(gpl_id, destdir = NULL) {
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
    result <- data.frame(probe_id = as.character(gpl_table[, id_col]),
                         gene_symbol = NA_character_, stringsAsFactors = FALSE)
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
    inner_join(probe2gene, by = "probe_id") %>%
    group_by(gene_symbol, sample) %>%
    summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
    mutate(expr = ifelse(is.nan(expr), NA_real_, expr)) %>%
    pivot_wider(names_from = sample, values_from = expr) %>%
    as.data.frame()

  rownames(expr_gene) <- expr_gene$gene_symbol
  expr_gene$gene_symbol <- NULL

  expr_gene
}
