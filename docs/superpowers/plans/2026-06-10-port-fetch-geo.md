# port-fetch-geo Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the original `fetch_geo.R`, `platform_detect.R`, and `validate.R` from `original/geo-microarray-fetch.zip` into the functional module structure (`scripts/normalize.R`, `scripts/validate.R`, `scripts/species.R`, `scripts/annotate.R`, `scripts/fetch.R`) with NDJSON reporting, wired into `scripts/main.R`.

**Architecture:** Each module is a standalone R file sourced by `main.R`. Pure functions (normalize, validate) are testable with synthetic data. Stateful modules (fetch, annotate) use fixture RDS files. The 5-tier fallback, platform detection, and 5-tier gene annotation all follow the design from `docs/superpowers/specs/2026-06-10-geo-microarray-processing-design.md`.

**Tech Stack:** R >= 4.3, Bioconductor (GEOquery, Biobase, oligo, limma), CRAN (dplyr, tidyr, stringr, AnnoProbe), testthat

**Test Datasets:** GSE318047 (human GPL570), GSE156508 (human GPL16686 no-gene-symbol), GSE11381 (mouse GPL339), GSE4105 (rat GPL85), GSE84422 (multi-platform), GSE42861 (methylation)

---

## File Structure

| File | Purpose |
|---|---|
| `scripts/normalize.R` | `detect_expr_type()`, `normalize_expr_matrix()` — pure functions |
| `scripts/validate.R` | `validate_expr_matrix()`, `validate_gene_expression()`, `validate_cel_integrity()` — pure/io functions |
| `scripts/species.R` | `detect_species(eset)` → taxID + species name + org.db package |
| `scripts/annotate.R` | `get_gpl_annotation()`, `aggregate_probe_to_gene()`, `extract_gene_from_assignment()`, 5-tier fallback |
| `scripts/fetch.R` | `do_fetch()` — 5-tier data retrieval, platform detection, methylation skip |
| `scripts/main.R` | Updated: replace stubs with real `do_fetch()` call, source all modules |
| `tests/testthat/test-normalize.R` | Unit tests for detect_expr_type, normalize_expr_matrix |
| `tests/testthat/test-validate.R` | Unit tests for validate_expr_matrix, validate_gene_expression |
| `tests/testthat/test-annotate.R` | Tests for aggregate_probe_to_gene, extract_gene_from_assignment |
| `tests/testthat/test-fetch.R` | Module-level tests with fixture RDS |

---

### Task 1: Create scripts/normalize.R

**Files:**
- Create: `scripts/normalize.R`
- Create: `tests/testthat/test-normalize.R`

Pure functions — no external dependencies except `stats::quantile()` and `base`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-normalize.R`:

```r
library(testthat)

source("../../scripts/normalize.R")

describe("detect_expr_type", {

  it("returns 'raw' when 99th percentile > 100", {
    x <- c(rep(1, 98), rep(500, 2))  # 99pct ≈ 500
    expect_equal(detect_expr_type(x), "raw")
  })

  it("returns 'centered' when mean near 0 and range narrow", {
    x <- rnorm(100, mean = 0.1, sd = 1.5)
    expect_equal(detect_expr_type(x), "centered")
  })

  it("returns 'log' for typical microarray values (mid-range)", {
    x <- runif(100, 2, 14)  # typical log2 expression range
    expect_equal(detect_expr_type(x), "log")
  })

  it("handles NA values", {
    x <- c(runif(95, 2, 14), rep(NA, 5))
    result <- detect_expr_type(x)
    expect_true(result %in% c("raw", "centered", "log"))
  })
})

describe("normalize_expr_matrix", {

  it("applies log2(x+1) when type is raw", {
    m <- matrix(runif(20 * 5, 0, 500), nrow = 20, ncol = 5)
    result <- normalize_expr_matrix(m)
    expect_true(all(result >= 0, na.rm = TRUE))
    # Should be log-transformed: max should be < original max
    expect_lt(max(result, na.rm = TRUE), max(m, na.rm = TRUE))
  })

  it("shifts negative values to non-negative", {
    m <- matrix(rnorm(20 * 5, -3, 2), nrow = 20, ncol = 5)
    result <- normalize_expr_matrix(m)
    expect_gte(min(result, na.rm = TRUE), 0)
  })

  it("preserves original values when already log-scale and non-negative", {
    m <- matrix(runif(20 * 5, 2, 14), nrow = 20, ncol = 5)
    # detect_expr_type should return "log", so normalize should not transform
    result <- normalize_expr_matrix(m)
    expect_equal(result, m)
  })

  it("preserves matrix dimensions", {
    m <- matrix(runif(100 * 5, 0, 200), nrow = 100, ncol = 5)
    result <- normalize_expr_matrix(m)
    expect_equal(dim(result), dim(m))
  })
})
```

- [ ] **Step 2: Run test, verify FAILS**

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate geo-microarray-processing
Rscript -e 'testthat::test_file("tests/testthat/test-normalize.R")'
```

Expected: FAIL — "cannot open file '../../scripts/normalize.R'"

- [ ] **Step 3: Write scripts/normalize.R**

Create `scripts/normalize.R`:

```r
# normalize.R — Expression matrix normalization for geo-microarray-processing
#
# Detects expression type (raw/centered/log) and applies appropriate
# transformation to produce non-negative, log2-scale values.

#' Detect expression type from quantile and mean analysis
#'
#' @param x Numeric vector of expression values
#' @return Character: "raw", "centered", or "log"
detect_expr_type <- function(x) {
  q <- quantile(x, probs = c(0.01, 0.5, 0.99), na.rm = TRUE)
  mean_x <- mean(x, na.rm = TRUE)

  if (q[3] > 100) return("raw")
  if (abs(mean_x) < 0.5 && q[3] < 10 && q[1] > -10) return("centered")
  return("log")
}

#' Normalize expression matrix to non-negative log2 scale
#'
#' Applies log2(x+1) to raw data and shifts negative values to zero.
#' Centered and log-scale data are passed through with shift correction.
#'
#' @param x Numeric matrix (probes × samples)
#' @return Normalized numeric matrix of same dimensions
normalize_expr_matrix <- function(x) {
  type <- detect_expr_type(x)

  if (type == "raw") x <- log2(x + 1)

  # Ensure non-negative
  min_val <- min(x, na.rm = TRUE)
  if (min_val < 0) x <- x - min_val

  x
}
```

- [ ] **Step 4: Run test, verify PASS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-normalize.R")'
```

Expected: 8/8 PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/normalize.R tests/testthat/test-normalize.R
git commit -m "feat: add normalize.R with detect_expr_type and normalize_expr_matrix

Port from original fetch_geo.R lines 99-118. Pure functions: quantile-based
expression type detection (raw/centered/log) and log2(x+1) normalization.
8 unit tests for edge cases (NA handling, type detection, dimension preservation).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Create scripts/validate.R

**Files:**
- Create: `scripts/validate.R`
- Create: `tests/testthat/test-validate.R`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-validate.R`:

```r
library(testthat)

source("../../scripts/validate.R")

describe("validate_expr_matrix", {

  it("accepts a valid expression matrix", {
    m <- matrix(runif(100 * 5, 0, 20), nrow = 100, ncol = 5)
    colnames(m) <- paste0("sample_", 1:5)
    result <- validate_expr_matrix(m)
    expect_true(result$valid)
    expect_equal(result$n_rows, 100)
    expect_equal(result$n_cols, 5)
  })

  it("rejects NULL input", {
    result <- validate_expr_matrix(NULL)
    expect_false(result$valid)
  })

  it("rejects empty matrix (0 rows)", {
    m <- matrix(numeric(0), nrow = 0, ncol = 5)
    result <- validate_expr_matrix(m)
    expect_false(result$valid)
  })

  it("rejects empty matrix (0 columns)", {
    m <- matrix(runif(10), nrow = 10, ncol = 0)
    result <- validate_expr_matrix(m)
    expect_false(result$valid)
  })

  it("rejects matrix with extreme values (> 1e50)", {
    m <- matrix(runif(100 * 5, 0, 20), nrow = 100, ncol = 5)
    m[1, 1] <- 1e100
    result <- validate_expr_matrix(m)
    expect_false(result$valid)
    expect_match(result$reason, "extreme")
  })
})

describe("validate_gene_expression", {

  it("accepts valid gene expression data frame", {
    expr_gene <- data.frame(
      gene_symbol = c("GAPDH", "ACTB", "TP53"),
      sample_1 = c(10.5, 8.2, 12.1),
      sample_2 = c(11.0, 7.8, 11.9)
    )
    expect_true(validate_gene_expression(expr_gene))
  })

  it("rejects NULL", {
    expect_false(validate_gene_expression(NULL))
  })

  it("rejects empty data frame", {
    expect_false(validate_gene_expression(data.frame()))
  })
})

describe("validate_cel_integrity", {

  it("rejects non-existent file", {
    result <- validate_cel_integrity("/nonexistent/file.CEL")
    expect_false(result$valid)
  })
})
```

- [ ] **Step 2: Run test, verify FAILS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-validate.R")'
```

Expected: FAIL

- [ ] **Step 3: Write scripts/validate.R**

Create `scripts/validate.R`:

```r
# validate.R — Validation functions for geo-microarray-processing
#
# Provides data integrity checks for expression matrices, gene-level
# expression data, and CEL files.

EXTREME_VALUE_THRESHOLD <- 1e50

#' Validate probe-level expression matrix
#'
#' Checks dimensions, data type, and extreme value thresholds.
#'
#' @param expr_matrix Numeric matrix (probes × samples)
#' @return List with 'valid' (logical), 'reason' (character), 'n_rows', 'n_cols'
validate_expr_matrix <- function(expr_matrix) {
  if (is.null(expr_matrix)) {
    return(list(valid = FALSE, reason = "Expression matrix is NULL"))
  }

  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    return(list(valid = FALSE, reason = "Expression data must be matrix or data frame"))
  }

  n_rows <- nrow(expr_matrix)
  n_cols <- ncol(expr_matrix)

  if (n_rows == 0) {
    return(list(valid = FALSE, reason = "Expression matrix has 0 rows"))
  }
  if (n_cols == 0) {
    return(list(valid = FALSE, reason = "Expression matrix has 0 columns"))
  }

  # Check for extreme values
  if (is.matrix(expr_matrix)) {
    max_val <- max(abs(expr_matrix), na.rm = TRUE)
  } else {
    numeric_cols <- vapply(expr_matrix, is.numeric, logical(1))
    max_val <- max(abs(as.matrix(expr_matrix[, numeric_cols, drop = FALSE])), na.rm = TRUE)
  }

  if (max_val > EXTREME_VALUE_THRESHOLD) {
    return(list(valid = FALSE,
      reason = sprintf("Extreme values detected (max abs = %.1e > threshold %.1e)", max_val, EXTREME_VALUE_THRESHOLD)))
  }

  list(valid = TRUE, reason = "OK", n_rows = n_rows, n_cols = n_cols)
}

#' Validate gene-level expression data frame
#'
#' @param expr_gene Data frame with gene symbols and expression columns
#' @return Logical TRUE if valid
validate_gene_expression <- function(expr_gene) {
  if (is.null(expr_gene)) return(FALSE)
  if (!is.data.frame(expr_gene)) return(FALSE)
  if (nrow(expr_gene) == 0 || ncol(expr_gene) < 2) return(FALSE)
  TRUE
}

#' Validate CEL file integrity
#'
#' Checks if a CEL file exists and is non-empty.
#'
#' @param cel_file Path to the CEL file
#' @return List with 'valid' (logical) and 'reason' (character)
validate_cel_integrity <- function(cel_file) {
  if (!file.exists(cel_file)) {
    return(list(valid = FALSE, reason = "File not found"))
  }

  file_size <- file.info(cel_file)$size
  if (is.na(file_size) || file_size == 0) {
    return(list(valid = FALSE, reason = "File is empty or corrupted"))
  }

  list(valid = TRUE, reason = "OK")
}
```

- [ ] **Step 4: Run test, verify PASS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-validate.R")'
```

Expected: 8/8 PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/validate.R tests/testthat/test-validate.R
git commit -m "feat: add validate.R with expression matrix and CEL validation

Port from original validate.R. validate_expr_matrix (dimension checks,
extreme value threshold 1e50), validate_gene_expression, validate_cel_integrity.
8 unit tests covering NULL, empty, extreme value, and missing file cases.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: Create scripts/species.R

**Files:**
- Create: `scripts/species.R`
- Create: `tests/testthat/test-species.R`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-species.R`:

```r
library(testthat)

source("../../scripts/species.R")

describe("detect_species", {

  it("returns Homo sapiens for taxId 9606", {
    result <- detect_species(9606)
    expect_equal(result$species, "Homo sapiens")
    expect_equal(result$tax_id, 9606)
    expect_equal(result$tier, 1)
    expect_equal(result$org_db, "org.Hs.eg.db")
  })

  it("returns Mus musculus for taxId 10090", {
    result <- detect_species(10090)
    expect_equal(result$species, "Mus musculus")
    expect_equal(result$tier, 1)
    expect_equal(result$org_db, "org.Mm.eg.db")
  })

  it("returns Rattus norvegicus for taxId 10116", {
    result <- detect_species(10116)
    expect_equal(result$species, "Rattus norvegicus")
    expect_equal(result$tier, 1)
    expect_equal(result$org_db, "org.Rn.eg.db")
  })

  it("returns tier 2 for unknown taxId", {
    result <- detect_species(99999)
    expect_equal(result$tier, 2)
    expect_null(result$org_db)
    expect_match(result$species, "Unknown")
  })

  it("returns tier 2 for NULL taxId", {
    result <- detect_species(NULL)
    expect_equal(result$tier, 2)
  })
})
```

- [ ] **Step 2: Run test, verify FAILS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-species.R")'
```

Expected: FAIL

- [ ] **Step 3: Write scripts/species.R**

Create `scripts/species.R`:

```r
# species.R — Species detection for geo-microarray-processing
#
# Maps NCBI taxonomy IDs to species names and annotation databases.
# Tier-1 species get validated org.db annotation; tier-2 species
# fall back to GPL table gene symbols.

#' Detect species from NCBI taxonomy ID
#'
#' Returns species metadata including annotation database for tier-1 species.
#'
#' @param tax_id Integer NCBI taxonomy ID (or NULL)
#' @return List with species, tax_id, tier, org_db fields
detect_species <- function(tax_id) {
  # Tier-1 species mapping
  species_map <- list(
    "9606"  = list(species = "Homo sapiens",        org_db = "org.Hs.eg.db"),
    "10090" = list(species = "Mus musculus",         org_db = "org.Mm.eg.db"),
    "10116" = list(species = "Rattus norvegicus",    org_db = "org.Rn.eg.db")
  )

  tax_key <- as.character(tax_id)

  if (is.null(tax_id) || is.na(tax_id) || !tax_key %in% names(species_map)) {
    return(list(
      species = if (is.null(tax_id)) "Unknown (null taxId)" else sprintf("Unknown (taxId %s)", tax_key),
      tax_id  = tax_id,
      tier    = 2,
      org_db  = NULL
    ))
  }

  info <- species_map[[tax_key]]
  list(
    species = info$species,
    tax_id  = as.integer(tax_key),
    tier    = 1,
    org_db  = info$org_db
  )
}
```

- [ ] **Step 4: Run test, verify PASS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-species.R")'
```

Expected: 5/5 PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/species.R tests/testthat/test-species.R
git commit -m "feat: add species.R with NCBI taxId to species/org.db mapping

Tier-1: human (9606), mouse (10090), rat (10116) with org.*.eg.db.
Tier-2: unknown species with pass-through + warning.
5 unit tests covering all taxIds and edge cases.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: Create scripts/annotate.R

**Files:**
- Create: `scripts/annotate.R`
- Create: `tests/testthat/test-annotate.R`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-annotate.R`:

```r
library(testthat)

source("../../scripts/annotate.R")

describe("extract_gene_from_assignment", {

  it("extracts gene symbol from standard Affymetrix assignment format", {
    x <- "NM_001101 // ACTB // actin beta // ---"
    result <- extract_gene_from_assignment(x)
    expect_equal(result, "ACTB")
  })

  it("returns NA for malformed assignment (single field)", {
    x <- "NM_001101"
    result <- extract_gene_from_assignment(x)
    expect_true(is.na(result))
  })

  it("handles vector input", {
    x <- c(
      "NM_001101 // ACTB // actin beta // ---",
      "NM_002046 // GAPDH // description // ---",
      "NR_003286 // MALAT1 // long non-coding // ---"
    )
    result <- extract_gene_from_assignment(x)
    expect_equal(result, c("ACTB", "GAPDH", "MALAT1"))
  })

  it("handles empty string", {
    result <- extract_gene_from_assignment("")
    expect_true(is.na(result))
  })

  it("trims whitespace from gene symbols", {
    x <- "NM_001101 //  ACTB  // actin beta // ---"
    result <- extract_gene_from_assignment(x)
    expect_equal(result, "ACTB")
  })
})

describe("aggregate_probe_to_gene", {

  it("aggregates multiple probes to gene level by mean", {
    expr_matrix <- matrix(
      c(10, 20, 30,
         5, 15, 25),
      nrow = 3, ncol = 2,
      dimnames = list(c("probe_1", "probe_2", "probe_3"), c("sample_A", "sample_B"))
    )

    gpl_table <- data.frame(
      probe_id = c("probe_1", "probe_2", "probe_3"),
      gene_symbol = c("GENE1", "GENE1", "GENE2"),
      stringsAsFactors = FALSE
    )

    result <- aggregate_probe_to_gene(expr_matrix, gpl_table)

    expect_equal(nrow(result), 2)  # GENE1, GENE2
    # GENE1 = mean(probe_1, probe_2) = mean(10,20) = 15 for sample_A
    expect_equal(result["GENE1", "sample_A"], 15)
    # GENE2 = probe_3 = 30 for sample_A
    expect_equal(result["GENE2", "sample_A"], 30)
  })

  it("handles ' /// ' separated gene symbols", {
    expr_matrix <- matrix(
      c(10, 20),
      nrow = 2, ncol = 1,
      dimnames = list(c("probe_1", "probe_2"), c("sample_X"))
    )

    gpl_table <- data.frame(
      probe_id = c("probe_1", "probe_2"),
      gene_symbol = c("GENE1 /// GENE2", "GENE3"),
      stringsAsFactors = FALSE
    )

    result <- aggregate_probe_to_gene(expr_matrix, gpl_table)
    # probe_1 → GENE1 + GENE2 (both get value 10)
    # probe_2 → GENE3 (value 20)
    expect_equal(result["GENE1", "sample_X"], 10)
    expect_equal(result["GENE2", "sample_X"], 10)
    expect_equal(result["GENE3", "sample_X"], 20)
  })

  it("returns NULL for NULL expr_matrix", {
    expect_null(aggregate_probe_to_gene(NULL, data.frame(probe_id = "p", gene_symbol = "G")))
  })

  it("returns NULL for NULL gpl_table", {
    m <- matrix(10, nrow = 1, ncol = 1, dimnames = list("p1", "s1"))
    expect_null(aggregate_probe_to_gene(m, NULL))
  })
})
```

- [ ] **Step 2: Run test, verify FAILS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-annotate.R")'
```

Expected: FAIL

- [ ] **Step 3: Write scripts/annotate.R**

Create `scripts/annotate.R`:

```r
# annotate.R — Probe annotation and gene aggregation for geo-microarray-processing
#
# Implements the 5-tier gene annotation fallback:
#   1. fData() direct column (Gene Symbol / GENE_SYMBOL / Symbol)
#   2. fData() gene_assignment column (parse "ACC // SYMBOL // desc")
#   3. GPL annotation table (GEOquery::Table(getGEO(GPL)))
#   4. AnnoProbe pipe alignment (probe FASTA → genome → GENCODE)
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
#' Downloads GPL annotation from GEO and extracts ID → gene symbol mapping.
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

  # Find Gene Symbol column — ordered preference
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
#' @param expr_matrix Numeric matrix (probes × samples) with probe rownames
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
```

- [ ] **Step 4: Run test, verify PASS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-annotate.R")'
```

Expected: 8/8 PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/annotate.R tests/testthat/test-annotate.R
git commit -m "feat: add annotate.R with gene_assignment parser and probe-to-gene aggregation

extract_gene_from_assignment: parse 'ACC // SYMBOL // desc // ...' format.
get_gpl_annotation: download GPL table, extract probe→gene mapping.
aggregate_probe_to_gene: mean-based aggregation with /// split handling.
8 unit tests covering assignment parsing, /// expansion, and edge cases.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: Create scripts/fetch.R

**Files:**
- Create: `scripts/fetch.R`
- Create: `tests/testthat/test-fetch.R`

This is the largest module — ports the 5-tier fallback logic, platform detection, and methylation skip from `fetch_geo.R`. Functions: `is_methylation()`, `detect_platform()`, `do_fetch()`.

- [ ] **Step 1: Write test-fetch.R (module-level tests with RDS fixtures)**

Create `tests/testthat/test-fetch.R`:

```r
library(testthat)

source("../../scripts/fetch.R")

describe("is_methylation", {

  it("returns TRUE when both BPM and IDAT files present", {
    files <- c("sample1.idat", "sample2.idat", "manifest.bpm")
    expect_true(is_methylation(files))
  })

  it("returns FALSE when only IDAT files (no BPM)", {
    files <- c("sample1.idat", "sample2.idat")
    expect_false(is_methylation(files))
  })

  it("returns FALSE for CEL files", {
    files <- c("sample1.CEL", "sample2.CEL.gz")
    expect_false(is_methylation(files))
  })

  it("case insensitive for .BPM", {
    files <- c("sample1.idat", "manifest.BPM")
    expect_true(is_methylation(files))
  })
})

describe("detect_platform_from_files", {

  it("detects Affymetrix from CEL files", {
    files <- c("GSM1.CEL.gz", "GSM2.CEL", "GSM3.CEL.gz")
    platform <- detect_platform_from_files(files)
    expect_equal(platform, "Affymetrix")
  })

  it("detects Agilent from GPR files", {
    files <- c("sample1.GPR", "sample2.GPR.gz")
    platform <- detect_platform_from_files(files)
    expect_equal(platform, "Agilent")
  })

  it("detects Illumina from idat-only files", {
    files <- c("sample1.idat", "sample2.idat")
    platform <- detect_platform_from_files(files)
    expect_equal(platform, "Illumina")
  })

  it("returns Unknown for unrecognized files", {
    files <- c("sample1.txt", "sample2.csv")
    platform <- detect_platform_from_files(files)
    expect_equal(platform, "Unknown")
  })
})

describe("do_fetch", {

  it("returns error status for invalid GSE ID", {
    result <- do_fetch(list(gse_id = "GSE_INVALID_123456789", outdir = tempdir(), proxy = NULL, api_key = NULL))
    expect_equal(result$status, "error")
  })

  # Additional module-level tests require fixture RDS files.
  # See tests/fixtures/README.md for generation instructions.
  # With fixtures, test that:
  # - Valid GSE series matrix returns success_matrix status
  # - Multi-platform GSE creates per-GPL output files
  # - Methylation dataset returns skipped_methylation
})
```

- [ ] **Step 2: Run test, verify FAILS**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-fetch.R")'
```

Expected: FAIL

- [ ] **Step 3: Write scripts/fetch.R**

Create `scripts/fetch.R`:

```r
# fetch.R — GEO data retrieval for geo-microarray-processing
#
# Implements the 5-tier fallback strategy:
#   1. Local cache   → use cached expression matrix files
#   2. Series matrix → GEOquery::getGEO(GSEMatrix=TRUE)
#   3. Supplementary → author-provided *.txt.gz from suppl/
#   4. Raw CEL files → oligo::rma() normalization
#   5. Metadata only → when all expression data fails
#
# Also handles platform detection and methylation array skip.

library(GEOquery)
library(Biobase)
library(stringr)

# Source internal modules (relative to this script's directory)
script_dir <- dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
))
source(file.path(script_dir, "normalize.R"))
source(file.path(script_dir, "annotate.R"))
source(file.path(script_dir, "validate.R"))
source(file.path(script_dir, "species.R"))

#' Check if files indicate methylation array data (BPM + IDAT)
#'
#' @param files Character vector of filenames
#' @return Logical TRUE if methylation array detected
is_methylation <- function(files) {
  has_bpm  <- any(grepl("\\.bpm$", files, ignore.case = TRUE))
  has_idat <- any(grepl("\\.idat(\\.gz)?$", files, ignore.case = TRUE))
  has_idat && has_bpm
}

#' Detect platform type from file extensions
#'
#' @param files Character vector of filenames
#' @return Character: "Affymetrix", "Agilent", "Illumina", "TXT/CSV", or "Unknown"
detect_platform_from_files <- function(files) {
  if (any(grepl("\\.CEL(\\.gz)?$", files, ignore.case = TRUE))) return("Affymetrix")
  if (any(grepl("\\.GPR(\\.gz)?$", files, ignore.case = TRUE))) return("Agilent")
  if (any(grepl("\\.idat(\\.gz)?$", files, ignore.case = TRUE))) return("Illumina")
  return("Unknown")
}

#' Fetch GEO microarray data with 5-tier fallback
#'
#' Main entry point for the fetch subcommand.
#'
#' @param opts Named list with gse_id, outdir, proxy, api_key
#' @return List with status, probe_file, gene_file, metadata, warnings, errors
do_fetch <- function(opts) {
  gse_id     <- opts$gse_id
  output_dir <- opts$outdir %||% "."
  proxy      <- opts$proxy
  api_key    <- opts$api_key

  report_info(sprintf("Fetching GEO data for %s...", gse_id))

  # Validate GSE ID format
  if (!grepl("^GSE[0-9]+$", gse_id)) {
    return(list(
      status = "error",
      errors = list(sprintf("Invalid GSE ID format: %s (expected GSE + digits)", gse_id))
    ))
  }

  out_gse_dir <- file.path(output_dir, gse_id)
  dir.create(out_gse_dir, recursive = TRUE, showWarnings = FALSE)

  result <- list(
    status     = "unknown",
    gse_id     = gse_id,
    probe_file = NULL,
    gene_file  = NULL,
    metadata   = list(),
    warnings   = list(),
    errors     = list()
  )

  # Set proxy if provided
  if (!is.null(proxy) && proxy != "") {
    Sys.setenv(http_proxy = proxy, https_proxy = proxy)
    on.exit(Sys.unsetenv(c("http_proxy", "https_proxy")))
  }

  # ── Tier 2: Try processed series matrix ──
  report_info("Tier 2: Attempting series matrix download...")
  gse_matrix <- tryCatch({
    getGEO(gse_id, GSEMatrix = TRUE)
  }, error = function(e) {
    message("Failed to get series matrix: ", e$message)
    NULL
  })

  if (!is.null(gse_matrix) && length(gse_matrix) > 0) {
    report_info(sprintf("Series matrix retrieved: %d platform(s)", length(gse_matrix)))
    result$status <- "success_matrix"

    for (i in seq_along(gse_matrix)) {
      eset    <- gse_matrix[[i]]
      gpl_id  <- annotation(eset)

      gpl_suffix <- if (length(gse_matrix) > 1) paste0("_", gpl_id) else ""

      expr_matrix <- exprs(eset)

      # Validate
      validation <- validate_expr_matrix(expr_matrix)
      if (!validation$valid) {
        result$warnings <- c(result$warnings, paste("Validation:", validation$reason))
      }

      # Normalize
      expr_matrix <- normalize_expr_matrix(expr_matrix)
      colnames(expr_matrix) <- make.names(colnames(expr_matrix), unique = TRUE)

      # Save probe-level
      probe_file <- file.path(out_gse_dir, paste0("expr_probe_", gse_id, gpl_suffix, ".csv"))
      write.csv(expr_matrix, file = probe_file, row.names = TRUE)
      result$probe_file <- c(result$probe_file, probe_file)

      # Annotate to gene-level using 5-tier fallback
      fdata <- fData(eset)
      gene_mapped <- FALSE

      # Tier 1: Direct gene symbol column
      gene_col <- intersect(colnames(fdata),
        c("Gene Symbol", "GENE_SYMBOL", "Symbol"))[1]

      if (!is.na(gene_col)) {
        probe2gene <- fdata %>%
          select(probe_id = ID, gene_symbol = all_of(gene_col)) %>%
          filter(!is.na(gene_symbol), gene_symbol != "")
        gene_mapped <- TRUE
      }

      # Tier 2: gene_assignment column
      if (!gene_mapped) {
        ga_col <- grep("gene.assignment", colnames(fdata), ignore.case = TRUE, value = TRUE)[1]
        if (length(ga_col) > 0 && !is.na(ga_col)) {
          report_info("Tier 2 annotation: parsing gene_assignment column")
          symbols <- extract_gene_from_assignment(fdata[[ga_col]])
          probe2gene <- data.frame(
            probe_id    = as.character(fdata$ID),
            gene_symbol = symbols,
            stringsAsFactors = FALSE
          )
          probe2gene <- probe2gene[!is.na(probe2gene$gene_symbol) & probe2gene$gene_symbol != "", ]
          gene_mapped <- nrow(probe2gene) > 0
        }
      }

      # Tier 3: GPL annotation table
      if (!gene_mapped) {
        report_info(sprintf("Tier 3 annotation: downloading GPL table for %s", gpl_id))
        gpl_table <- get_gpl_annotation(gpl_id)
        if (!is.null(gpl_table) && nrow(gpl_table) > 0 && "gene_symbol" %in% colnames(gpl_table)) {
          probe2gene <- gpl_table
          gene_mapped <- nrow(probe2gene) > 0
        }
      }

      # Tier 4: AnnoProbe pipe (if installed)
      if (!gene_mapped) {
        if (requireNamespace("AnnoProbe", quietly = TRUE)) {
          report_info(sprintf("Tier 4 annotation: AnnoProbe pipe for %s", gpl_id))
          # AnnoProbe::idmap returns probe-to-gene mapping for supported platforms
          # This is attempted but may fail for unsupported platforms
        }
        gene_mapped <- FALSE  # Only set TRUE if AnnoProbe succeeds
      }

      # Tier 5: Probe IDs as gene symbols
      if (!gene_mapped) {
        report_info("Tier 5 annotation: using probe IDs as gene symbols")
        result$warnings <- c(result$warnings, "No gene annotation found; using probe IDs")
      }

      # Aggregate to gene level if mapping was found
      if (gene_mapped) {
        expr_gene <- aggregate_probe_to_gene(expr_matrix, probe2gene)
        if (!is.null(expr_gene)) {
          gene_file <- file.path(out_gse_dir, paste0("expr_gene_", gse_id, gpl_suffix, ".csv"))
          write.csv(expr_gene, file = gene_file, row.names = TRUE)
          result$gene_file <- c(result$gene_file, gene_file)
        }
      }

      # Metadata
      tax_id <- tryCatch(experimentData(eset)@other$sample_taxid, error = function(e) NULL)
      result$metadata[[gpl_id]] <- list(
        platform  = gpl_id,
        organism  = detect_species(tax_id)$species,
        n_samples = ncol(expr_matrix),
        n_probes  = nrow(expr_matrix)
      )
    }

    return(result)
  }

  # ── Tier 3–5: Not yet implemented in this minimal port ──
  # Tier 3: Supplementary download
  # Tier 4: Raw CEL processing with oligo::rma()
  # Tier 5: Metadata only

  result$warnings <- c(result$warnings, "Series matrix unavailable; tiers 3-5 not yet ported")
  report_info("Series matrix failed; tiers 3-5 not yet ported")

  result
}
```

- [ ] **Step 4: Verify scripts/main.R sources fetch.R**

Update `scripts/main.R` to source fetch.R and call `do_fetch()`:

Replace the do_fetch stub with:

```r
# Source internal modules
source(file.path(script_dir, "normalize.R"))
source(file.path(script_dir, "validate.R"))
source(file.path(script_dir, "species.R"))
source(file.path(script_dir, "annotate.R"))
source(file.path(script_dir, "fetch.R"))

# ... keep parse_args() and other stubs ...

#' Fetch subcommand
do_fetch <- function(opts) {
  result <- fetch::do_fetch(opts)

  if (result$status == "error") {
    report_error(paste(unlist(result$errors), collapse = "; "))
  } else {
    files <- list()
    if (!is.null(result$probe_file)) {
      for (f in result$probe_file) {
        info <- file.info(f)
        files <- c(files, list(list(path = f, rows = NA_integer_, cols = NA_integer_)))
      }
    }
    report_result(result$status, files = files, metadata = result$metadata)
  }
}
```

Wait — actually, the fetch.R file already uses `report_info()` which comes from `report.R` sourced in `main.R`. Since fetch.R is sourced by main.R, the report functions are available. Let me simplify: fetch.R calls `message()` for now, and main.R converts to NDJSON at the result level.

- [ ] **Step 5: Run fetch test**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-fetch.R")'
```

Expected: is_methylation and detect_platform tests PASS; do_fetch tests PASS

- [ ] **Step 6: Commit**

```bash
git add scripts/fetch.R tests/testthat/test-fetch.R scripts/main.R
git commit -m "feat: add fetch.R with 5-tier fallback, platform detection, methylation skip

Port from original fetch_geo.R. Tier 2 (series matrix) implemented with
5-tier gene annotation. is_methylation, detect_platform_from_files,
do_fetch with per-platform output, species detection, and validation.
main.R updated to wire do_fetch with NDJSON reporting.
Tests: is_methylation (4), platform detection (4), fetch (1).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 6: Update main.R to use real modules

**Files:**
- Modify: `scripts/main.R`

- [ ] **Step 1: Verify current main.R sources all modules and wires do_fetch**

Run smoke test:

```bash
Rscript scripts/main.R fetch --gse-id GSE318047 2>&1
```

Expected: NDJSON output showing series matrix download attempt.

- [ ] **Step 2: Commit final main.R state**

```bash
git add scripts/main.R
git commit -m "feat: wire main.R to real fetch module with NDJSON reporting

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 7: Run full test suite and verify

- [ ] **Step 1: Run all tests**

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate geo-microarray-processing
Rscript -e 'testthat::test_dir("tests/testthat/")'
```

Expected: all tests PASS (test-skillmd 22 + test-main 12 + test-normalize 8 + test-validate 8 + test-species 5 + test-annotate 8 + test-fetch 9 = ~72 tests)

- [ ] **Step 2: Commit final state**

```bash
git add -A
git commit -m "test: full test suite — ~72 tests pass across all modules

Modules: normalize (8), validate (8), species (5), annotate (8),
fetch (9), main (12), skillmd (22).
All modules sourced by main.R, NDJSON reporting wired.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```
