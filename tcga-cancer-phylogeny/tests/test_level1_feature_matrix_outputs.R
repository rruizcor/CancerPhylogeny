#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

paths <- c(
  matrix = file.path(root, "results", "tables", "level1_feature_matrix.tsv"),
  support = file.path(root, "results", "tables", "level1_feature_support.tsv"),
  qc = file.path(root, "results", "tables", "level1_feature_matrix_qc_summary.tsv")
)
missing_or_empty <- paths[!file.exists(paths) | file.info(paths)$size == 0]
if (length(missing_or_empty) > 0) {
  fail(sprintf("Missing or empty Level 1 feature-matrix outputs: %s", paste(missing_or_empty, collapse = ", ")))
}

matrix <- read.delim(paths[["matrix"]], stringsAsFactors = FALSE, check.names = FALSE)
support <- read.delim(paths[["support"]], stringsAsFactors = FALSE, check.names = FALSE)
qc <- read.delim(paths[["qc"]], stringsAsFactors = FALSE, check.names = FALSE)
projects <- read.delim(file.path(root, "data", "interim", "tcga_projects.tsv"), stringsAsFactors = FALSE)

required_annotation_cols <- c("project_id", "project_code", "cancer_name", "broad_group", "level2_group", "lineage_notes")
required_feature_cols <- c(
  "tmb_proxy_median_nonsynonymous_count",
  "tmb_proxy_median_total_mutation_count",
  "mutation_prevalence_TP53",
  "mutation_prevalence_KRAS",
  "median_purity",
  "median_ploidy",
  "median_aneuploidy_score",
  "median_total_arm_alteration_count",
  "median_proliferation_score",
  "median_immune_inflammatory_score",
  "median_EMT_score",
  "median_lineage_or_tissue_PC1",
  "median_lineage_or_tissue_PC2"
)
required_support_cols <- c("project_id", "project_code", "mutation_n_samples", "copy_number_n_samples", "expression_n_samples")

if (!all(required_annotation_cols %in% names(matrix))) fail("Level 1 feature matrix is missing required annotation columns")
if (!all(required_feature_cols %in% names(matrix))) fail("Level 1 feature matrix is missing required feature columns")
if (!all(required_support_cols %in% names(support))) fail("Level 1 support table is missing required support columns")
if (nrow(matrix) != nrow(projects)) fail("Level 1 feature matrix does not have one row per expected project")
if (anyDuplicated(matrix$project_id)) fail("Level 1 feature matrix has duplicate project IDs")
if (!setequal(matrix$project_id, projects$project_id)) fail("Level 1 feature matrix project IDs do not match expected TCGA projects")

feature_cols <- setdiff(names(matrix), required_annotation_cols)
for (column in feature_cols) {
  if (any(!is.na(matrix[[column]]) & is.na(suppressWarnings(as.numeric(matrix[[column]]))))) {
    fail(sprintf("Feature column is not numeric/coercible: %s", column))
  }
}
feature_missingness <- vapply(feature_cols, function(column) mean(is.na(matrix[[column]])), numeric(1))
if (any(feature_missingness > 0.30)) {
  fail(sprintf("Retained feature columns exceed missingness threshold: %s", paste(names(feature_missingness)[feature_missingness > 0.30], collapse = ", ")))
}

if (!all(c("project_rows", "retained_feature_columns", "dropped_feature_columns", "max_feature_missing_fraction") %in% qc$metric)) {
  fail("Level 1 feature-matrix QC summary is missing required global metrics")
}
if (!any(qc$qc_section == "per_feature" & qc$metric == "retained")) {
  fail("Level 1 feature-matrix QC summary is missing per-feature retained flags")
}

cat("Level 1 feature-matrix output checks passed\n")
