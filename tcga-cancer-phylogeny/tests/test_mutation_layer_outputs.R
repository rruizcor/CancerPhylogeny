#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

required_paths <- c(
  file.path(root, "data", "processed", "mutations", "mc3_somatic_mutations.parquet"),
  file.path(root, "data", "processed", "features", "tmb_by_sample.tsv"),
  file.path(root, "data", "processed", "features", "mutation_prevalence_by_project.tsv"),
  file.path(root, "data", "processed", "features", "mutation_prevalence_by_project_driver_only.tsv"),
  file.path(root, "results", "tables", "mutation_layer_qc_summary.tsv")
)

missing_or_empty <- required_paths[!file.exists(required_paths) | file.info(required_paths)$size == 0]
if (length(missing_or_empty) > 0) {
  fail(sprintf("Missing or empty mutation-layer outputs: %s", paste(missing_or_empty, collapse = ", ")))
}

tmb <- read.delim(file.path(root, "data", "processed", "features", "tmb_by_sample.tsv"), stringsAsFactors = FALSE)
prevalence <- read.delim(file.path(root, "data", "processed", "features", "mutation_prevalence_by_project.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
driver_prevalence <- read.delim(file.path(root, "data", "processed", "features", "mutation_prevalence_by_project_driver_only.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
qc <- read.delim(file.path(root, "results", "tables", "mutation_layer_qc_summary.tsv"), stringsAsFactors = FALSE)

required_tmb_columns <- c(
  "sample_barcode",
  "patient_barcode",
  "project_id",
  "project_code",
  "nonsynonymous_count",
  "total_mutation_count",
  "metric_label"
)
if (!all(required_tmb_columns %in% names(tmb))) {
  fail("TMB table is missing required columns")
}

if (!all(tmb$metric_label == "mutation_count_proxy_no_callable_territory")) {
  fail("TMB table metric_label should indicate mutation-count proxy status")
}

if (!all(c("project_id", "project_code") %in% names(prevalence))) {
  fail("Project-level prevalence table is missing project identifiers")
}

if (!all(c("project_id", "project_code") %in% names(driver_prevalence))) {
  fail("Driver-only prevalence table is missing project identifiers")
}

required_qc_metrics <- c(
  "total_mutation_records",
  "unique_tumor_samples",
  "unique_patients",
  "tcga_projects_detected",
  "expected_tcga_projects_missing",
  "ref_alt_counts_present",
  "hgvsp_short_present"
)
if (!all(required_qc_metrics %in% qc$metric)) {
  fail("QC summary is missing required global metrics")
}

cat("Mutation-layer output checks passed\n")
