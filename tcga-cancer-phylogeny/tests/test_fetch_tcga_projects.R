#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

projects_path <- file.path(root, "data", "interim", "tcga_projects.tsv")
config_path <- file.path(root, "config", "tcga_projects.yaml")

if (!file.exists(projects_path)) {
  fail(sprintf("Missing fetched project metadata: %s", projects_path))
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  fail("Package 'yaml' is required for fetch tests")
}

configured <- yaml::read_yaml(config_path)$projects
expected_project_ids <- vapply(configured, `[[`, character(1), "project_id")
actual <- read.delim(projects_path, stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c(
  "project_id",
  "project_code",
  "disease_type",
  "primary_site",
  "name",
  "program",
  "case_count",
  "file_count"
)

if (!all(required_columns %in% names(actual))) {
  fail("Fetched project metadata is missing one or more required columns")
}

if (nrow(actual) != length(expected_project_ids)) {
  fail(sprintf(
    "Expected %d TCGA projects but found %d rows",
    length(expected_project_ids),
    nrow(actual)
  ))
}

if (!identical(actual$project_id, expected_project_ids)) {
  fail("Fetched project IDs do not match the configured project list and order")
}

if (!all(actual$program == "TCGA")) {
  fail("Fetched project metadata includes non-TCGA projects")
}

if (any(is.na(actual$disease_type) | actual$disease_type == "")) {
  fail("One or more projects are missing disease_type")
}

if (any(is.na(actual$primary_site) | actual$primary_site == "")) {
  fail("One or more projects are missing primary_site")
}

cat("Fetched TCGA project metadata checks passed\n")
