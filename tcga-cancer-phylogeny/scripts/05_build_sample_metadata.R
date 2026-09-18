#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

build_sample_metadata <- function(clinical_path, mutation_path, copy_number_path, expression_path) {
  not_implemented(
    "Sample metadata construction",
    "join clinical, mutation, copy-number, and expression availability into one sample-level manifest"
  )
}

main <- function() {
  root <- find_project_root()
  log_info("Checking expected upstream sample metadata inputs")
  expected <- c(
    project_path("data", "interim", "tcga_projects.tsv", root = root),
    project_path("data", "processed", "mutations", "mc3_somatic_mutations.parquet", root = root)
  )
  missing <- expected[!file.exists(expected)]
  if (length(missing) > 0) {
    abort(sprintf("Missing upstream inputs: %s", paste(missing, collapse = ", ")))
  }
  build_sample_metadata(expected[[1]], expected[[2]], NA_character_, NA_character_)
}

if (identical(environment(), globalenv())) {
  main()
}
