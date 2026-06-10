# species.R — Species detection for geo-microarray-processing
#
# Maps NCBI taxonomy IDs to human-readable species names.
# Currently recognizes the three most common species in GEO expression data.
# All species go through the same 5-tier gene annotation fallback chain
# regardless of taxId — the species name is for metadata reporting only.

#' Detect species from NCBI taxonomy ID
#'
#' @param tax_id Integer NCBI taxonomy ID (or NULL)
#' @return List with species, tax_id, tier fields
detect_species <- function(tax_id) {
  species_map <- list(
    "9606"  = "Homo sapiens",
    "10090" = "Mus musculus",
    "10116" = "Rattus norvegicus"
  )

  tax_key <- as.character(tax_id)

  if (is.null(tax_id) || is.na(tax_id) || !tax_key %in% names(species_map)) {
    return(list(
      species = if (is.null(tax_id)) "Unknown" else paste0("taxId:", tax_key),
      tax_id  = tax_id,
      tier    = 2
    ))
  }

  list(
    species = species_map[[tax_key]],
    tax_id  = as.integer(tax_key),
    tier    = 1
  )
}
