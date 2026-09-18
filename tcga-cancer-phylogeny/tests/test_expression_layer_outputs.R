#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

paths <- c(
  matrix = file.path(root, "data", "processed", "expression", "expression_matrix_by_sample.tsv.gz"),
  sample = file.path(root, "data", "processed", "features", "expression_pathway_scores_by_sample.tsv"),
  project = file.path(root, "data", "processed", "features", "expression_pathway_scores_by_project.tsv"),
  pca = file.path(root, "data", "processed", "features", "expression_pca_by_sample.tsv"),
  variance = file.path(root, "data", "processed", "features", "expression_pca_variance.tsv"),
  qc = file.path(root, "results", "tables", "expression_layer_qc_summary.tsv")
)

missing_or_empty <- paths[!file.exists(paths) | file.info(paths)$size == 0]
if (length(missing_or_empty) > 0) {
  fail(sprintf("Missing or empty expression-layer outputs: %s", paste(missing_or_empty, collapse = ", ")))
}

sample_scores <- read.delim(paths[["sample"]], stringsAsFactors = FALSE)
project_scores <- read.delim(paths[["project"]], stringsAsFactors = FALSE)
pca_scores <- read.delim(paths[["pca"]], stringsAsFactors = FALSE)
variance <- read.delim(paths[["variance"]], stringsAsFactors = FALSE)
qc <- read.delim(paths[["qc"]], stringsAsFactors = FALSE)
tmb <- read.delim(file.path(root, "data", "processed", "features", "tmb_by_sample.tsv"), stringsAsFactors = FALSE)

required_score_cols <- c(
  "sample_barcode",
  "patient_barcode",
  "project_id",
  "project_code",
  "proliferation_score",
  "immune_inflammatory_score",
  "interferon_gamma_score",
  "cytotoxic_t_cell_score",
  "stromal_score",
  "epithelial_score",
  "EMT_score",
  "hypoxia_score",
  "cell_cycle_score",
  "DNA_repair_score",
  "angiogenesis_score",
  "lineage_or_tissue_PC1",
  "lineage_or_tissue_PC2",
  "source",
  "notes"
)
required_project_cols <- c(
  "project_id",
  "project_code",
  "n_samples",
  "median_proliferation_score",
  "median_immune_inflammatory_score",
  "median_interferon_gamma_score",
  "median_cytotoxic_t_cell_score",
  "median_stromal_score",
  "median_epithelial_score",
  "median_EMT_score",
  "median_hypoxia_score",
  "median_cell_cycle_score",
  "median_DNA_repair_score",
  "median_angiogenesis_score",
  "median_lineage_or_tissue_PC1",
  "median_lineage_or_tissue_PC2"
)

if (!all(required_score_cols %in% names(sample_scores))) fail("Sample-level expression score table is missing required columns")
if (!all(required_project_cols %in% names(project_scores))) fail("Project-level expression score table is missing required columns")
if (!all(c("sample_barcode", "patient_barcode", "project_id", "project_code", paste0("PC", 1:20)) %in% names(pca_scores))) {
  fail("Expression PCA sample table is missing required columns")
}
if (!all(c("PC", "percent_variance_explained", "pca_basis") %in% names(variance))) {
  fail("Expression PCA variance table is missing required columns")
}

if (nrow(sample_scores) == 0 || nrow(project_scores) == 0 || nrow(pca_scores) == 0 || nrow(variance) == 0 || nrow(qc) == 0) {
  fail("One or more expression output tables are completely empty")
}

if (any(!grepl("^TCGA-[A-Z0-9]+$", sample_scores$project_id))) fail("Sample-level expression table has invalid TCGA project IDs")
if (any(!grepl("^TCGA-[A-Z0-9]+$", project_scores$project_id))) fail("Project-level expression table has invalid TCGA project IDs")
if (!all(sample_scores$sample_barcode %in% tmb$sample_barcode)) {
  fail("Expression sample barcodes do not align with the mutation-layer sample table")
}

numeric_columns <- c(
  grep("_score$|^lineage_or_tissue_PC", names(sample_scores), value = TRUE),
  grep("^median_|^missingness_", names(project_scores), value = TRUE),
  paste0("PC", 1:20),
  "percent_variance_explained"
)
for (column in intersect(numeric_columns, names(sample_scores))) {
  if (any(!is.na(sample_scores[[column]]) & is.na(suppressWarnings(as.numeric(sample_scores[[column]]))))) {
    fail(sprintf("Sample score column is not numeric/coercible: %s", column))
  }
}
for (column in intersect(numeric_columns, names(project_scores))) {
  if (any(!is.na(project_scores[[column]]) & is.na(suppressWarnings(as.numeric(project_scores[[column]]))))) {
    fail(sprintf("Project score column is not numeric/coercible: %s", column))
  }
}
for (column in intersect(numeric_columns, names(pca_scores))) {
  if (any(!is.na(pca_scores[[column]]) & is.na(suppressWarnings(as.numeric(pca_scores[[column]]))))) {
    fail(sprintf("PCA column is not numeric/coercible: %s", column))
  }
}
if (any(!is.na(variance$percent_variance_explained) & is.na(suppressWarnings(as.numeric(variance$percent_variance_explained))))) {
  fail("PCA variance values are not numeric/coercible")
}

required_qc_metrics <- c(
  "expression_source_file_status",
  "expression_data_type",
  "normalization_method",
  "genes_available_in_scoring_matrix",
  "expression_samples_processed",
  "tcga_projects_with_expression_data",
  "mutation_layer_samples",
  "samples_overlapping_purity_ploidy_layer",
  "samples_overlapping_aneuploidy_layer",
  "pca_status"
)
if (!all(required_qc_metrics %in% qc$metric)) {
  fail("Expression QC summary is missing required metrics")
}

cat("Expression-layer output checks passed\n")
