#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

query_clinical_metadata <- function(project_ids) {
  not_implemented(
    "Clinical metadata download",
    "use TCGAbiolinks::GDCquery_clinic or the GDC API with caching per TCGA project"
  )
}

harmonize_clinical_metadata <- function(clinical_tables) {
  not_implemented(
    "Clinical metadata harmonization",
    "standardize case, sample, project, disease, primary site, and subtype fields"
  )
}

main <- function() {
  root <- find_project_root()
  projects <- load_tcga_project_config(root)
  log_info(sprintf("Clinical metadata planned for %d projects", nrow(projects)))
  query_clinical_metadata(projects$project_id)
}

if (identical(environment(), globalenv())) {
  main()
}
