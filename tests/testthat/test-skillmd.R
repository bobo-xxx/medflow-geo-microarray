library(testthat)

# helpers in same directory as this test file
source("helpers.R")

# Path to SKILL.md relative to tests/testthat/
skillmd_path <- "../../SKILL.md"

describe("SKILL.md", {

  it("exists at the package root", {
    expect_true(file.exists(skillmd_path))
  })

  it("has valid YAML frontmatter", {
    content <- readLines(skillmd_path, warn = FALSE)

    # Find YAML frontmatter delimiters
    dashes <- which(content == "---")
    expect_gte(length(dashes), 2,
      "SKILL.md must have opening and closing '---' YAML delimiters")

    yaml_lines <- content[(dashes[1] + 1):(dashes[2] - 1)]
    expect_gt(length(yaml_lines), 0,
      "YAML frontmatter must not be empty")

    # Write YAML to temp file and parse with yaml package
    tmp <- tempfile(fileext = ".yaml")
    writeLines(yaml_lines, tmp)
    frontmatter <- yaml::read_yaml(tmp)

    expect_type(frontmatter, "list")
  })

  it("declares required frontmatter fields", {
    content <- readLines(skillmd_path, warn = FALSE)
    dashes <- which(content == "---")
    yaml_lines <- content[(dashes[1] + 1):(dashes[2] - 1)]
    tmp <- tempfile(fileext = ".yaml")
    writeLines(yaml_lines, tmp)
    fm <- yaml::read_yaml(tmp)

    required <- c("name", "description", "type", "inputs", "outputs",
                  "entry", "parameters", "exceptions", "hardware")
    for (field in required) {
      expect_true(field %in% names(fm),
        sprintf("SKILL.md frontmatter must have '%s' field", field))
    }
  })

  it("has type 'standard'", {
    fm <- parse_skillmd_frontmatter(skillmd_path)
    expect_equal(fm$type, "standard")
  })

  it("has entry pointing to scripts/main.R", {
    fm <- parse_skillmd_frontmatter(skillmd_path)
    expect_equal(fm$entry, "scripts/main.R")
  })

  it("has parameters for subcommand, gse-id, outdir, proxy, api-key", {
    fm <- parse_skillmd_frontmatter(skillmd_path)
    param_names <- vapply(fm$parameters, `[[`, "", "name")
    expected <- c("subcommand", "--gse-id", "--outdir", "--input",
                  "--proxy", "--api-key")
    for (p in expected) {
      expect_true(p %in% param_names,
        sprintf("SKILL.md must declare parameter '%s'", p))
    }
  })

  it("has 6 exception patterns", {
    fm <- parse_skillmd_frontmatter(skillmd_path)
    expect_gte(length(fm$exceptions), 6)
  })
})
