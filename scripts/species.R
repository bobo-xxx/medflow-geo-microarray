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
