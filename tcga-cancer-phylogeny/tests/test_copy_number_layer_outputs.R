#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

paths <- c(
  purity = file.path(root, "data", "processed", "features", "purity_ploidy_by_sample.tsv"),
  aneuploidy = file.path(root, "data", "processed", "features", "aneuploidy_by_sample.tsv"),
  project = file.path(root, "data", "processed", "features", "aneuploidy_by_project.tsv"),
  qc = file.path(root, "results", "tables", "copy_number_layer_qc_summary.tsv")
)

missing_or_empty <- paths[!file.exists(paths) | file.info(paths)$size == 0]
if (length(missing_or_empty) > 0) {
  fail(sprintf("Missing or empty copy-number-layer outputs: %s", paste(missing_or_empty, collapse = ", ")))
}

purity <- read.delim(paths[["purity"]], stringsAsFactors = FALSE)
aneuploidy <- read.delim(paths[["aneuploidy"]], stringsAsFactors = FALSE)
project <- read.delim(paths[["project"]], stringsAsFactors = FALSE)
qc <- read.delim(paths[["qc"]], stringsAsFactors = FALSE)
tmb <- read.delim(file.path(root, "data", "processed", "features", "tmb_by_sample.tsv"), stringsAsFactors = FALSE)

required_purity_cols <- c(
  "sample_barcode",
  "patient_barcode",
  "project_id",
  "project_code",
  "purity",
  "ploidy",
  "whole_genome_doubling_status",
  "source",
  "notes"
)
required_aneuploidy_cols <- c(
  "sample_barcode",
  "patient_barcode",
  "project_id",
  "project_code",
  "aneuploidy_score",
  "fraction_genome_altered",
  "arm_gain_count",
  "arm_loss_count",
  "total_arm_alteration_count",
  "copy_number_complexity_score",
  "source",
  "notes"
)
required_project_cols <- c(
  "project_id",
  "project_code",
  "n_samples",
  "median_purity",
  "median_ploidy",
  "median_aneuploidy_score",
  "median_fraction_genome_altered",
  "median_arm_gain_count",
  "median_arm_loss_count",
  "median_total_arm_alteration_count",
  "median_copy_number_complexity_score",
  "missingness_purity",
  "missingness_ploidy",
  "missingness_aneuploidy_score",
  "missingness_fraction_genome_altered"
)

if (!all(required_purity_cols %in% names(purity))) fail("Purity/ploidy table is missing required columns")
if (!all(required_aneuploidy_cols %in% names(aneuploidy))) fail("Aneuploidy table is missing required columns")
if (!all(required_project_cols %in% names(project))) fail("Project-level aneuploidy table is missing required columns")

if (nrow(purity) == 0 || nrow(aneuploidy) == 0 || nrow(project) == 0 || nrow(qc) == 0) {
  fail("One or more copy-number output tables are completely empty")
}

if (any(!grepl("^TCGA-[A-Z0-9]+$", purity$project_id))) fail("Purity/ploidy table has invalid TCGA project IDs")
if (any(!grepl("^TCGA-[A-Z0-9]+$", aneuploidy$project_id))) fail("Aneuploidy table has invalid TCGA project IDs")
if (any(!grepl("^TCGA-[A-Z0-9]+$", project$project_id))) fail("Project table has invalid TCGA project IDs")

if (!all(purity$sample_barcode %in% tmb$sample_barcode)) {
  fail("Purity/ploidy sample barcodes do not align with the mutation-layer sample table")
}
if (!all(aneuploidy$sample_barcode %in% tmb$sample_barcode)) {
  fail("Aneuploidy sample barcodes do not align with the mutation-layer sample table")
}

numeric_columns <- c(
  "purity",
  "ploidy",
  "aneuploidy_score",
  "fraction_genome_altered",
  "arm_gain_count",
  "arm_loss_count",
  "total_arm_alteration_count",
  "copy_number_complexity_score"
)
for (column in intersect(numeric_columns, names(purity))) {
  if (any(!is.na(purity[[column]]) & is.na(suppressWarnings(as.numeric(purity[[column]]))))) {
    fail(sprintf("Column is not numeric/coercible: %s", column))
  }
}
for (column in intersect(numeric_columns, names(aneuploidy))) {
  if (any(!is.na(aneuploidy[[column]]) & is.na(suppressWarnings(as.numeric(aneuploidy[[column]]))))) {
    fail(sprintf("Column is not numeric/coercible: %s", column))
  }
}
for (column in grep("^missingness_", names(project), value = TRUE)) {
  if (any(!is.na(project[[column]]) & is.na(suppressWarnings(as.numeric(project[[column]]))))) {
    fail(sprintf("Missingness column is not numeric/coercible: %s", column))
  }
}

required_qc_metrics <- c(
  "absolute_purity_ploidy_file_status",
  "arm_calls_aneuploidy_file_status",
  "samples_with_purity",
  "samples_with_ploidy",
  "samples_with_aneuploidy_score",
  "samples_with_arm_level_calls",
  "samples_with_purity_ploidy_data",
  "samples_with_aneuploidy_or_copy_number_data",
  "mutation_layer_samples",
  "percent_mutation_layer_samples_with_purity_ploidy",
  "percent_mutation_layer_samples_with_copy_number_features"
)
if (!all(required_qc_metrics %in% qc$metric)) {
  fail("QC summary is missing required metrics")
}

cat("Copy-number-layer output checks passed\n")
