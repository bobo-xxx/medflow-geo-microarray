library(testthat)

source("../../scripts/species.R")

describe("detect_species", {

  it("returns Homo sapiens for taxId 9606", {
    result <- detect_species(9606)
    expect_equal(result$species, "Homo sapiens")
    expect_equal(result$tax_id, 9606)
    expect_equal(result$tier, 1)
  })

  it("returns Mus musculus for taxId 10090", {
    result <- detect_species(10090)
    expect_equal(result$species, "Mus musculus")
    expect_equal(result$tier, 1)
  })

  it("returns Rattus norvegicus for taxId 10116", {
    result <- detect_species(10116)
    expect_equal(result$species, "Rattus norvegicus")
    expect_equal(result$tier, 1)
  })

  it("returns tier 2 for unknown taxId", {
    result <- detect_species(99999)
    expect_equal(result$tier, 2)
    expect_match(result$species, "taxId")
  })

  it("returns tier 2 for NULL taxId", {
    result <- detect_species(NULL)
    expect_equal(result$tier, 2)
  })
})
