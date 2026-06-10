#' Parse SKILL.md frontmatter YAML
#' @param path Path to SKILL.md
#' @return List of frontmatter fields
parse_skillmd_frontmatter <- function(path = "../../node/SKILL.md") {
  content <- readLines(path, warn = FALSE)
  dashes <- which(content == "---")
  yaml_lines <- content[(dashes[1] + 1):(dashes[2] - 1)]
  tmp <- tempfile(fileext = ".yaml")
  writeLines(yaml_lines, tmp)
  yaml::read_yaml(tmp)
}
