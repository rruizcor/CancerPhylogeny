#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

paths <- c(
  scaled = file.path(root, "results", "tables", "level1_feature_matrix_scaled.tsv"),
  gower_distance = file.path(root, "results", "tables", "level1_gower_distance_matrix.tsv"),
  correlation_distance = file.path(root, "results", "tables", "level1_correlation_distance_matrix.tsv"),
  bootstrap = file.path(root, "results", "tables", "level1_bootstrap_cluster_stability.tsv"),
  qc = file.path(root, "results", "tables", "level1_tree_qc_summary.tsv"),
  gower_tree = file.path(root, "results", "trees", "level1_pan_cancer_gower_hclust_tree.nwk"),
  correlation_tree = file.path(root, "results", "trees", "level1_pan_cancer_correlation_hclust_tree.nwk"),
  gower_figure = file.path(root, "results", "figures", "level1_pan_cancer_gower_tree.pdf"),
  correlation_figure = file.path(root, "results", "figures", "level1_pan_cancer_correlation_tree.pdf"),
  heatmap_figure = file.path(root, "results", "figures", "level1_pan_cancer_tree_with_trait_heatmap.pdf"),
  pca_figure = file.path(root, "results", "figures", "level1_pan_cancer_feature_pca.pdf")
)

missing_or_empty <- paths[!file.exists(paths) | file.info(paths)$size == 0]
if (length(missing_or_empty) > 0) {
  fail(sprintf("Missing or empty Level 1 tree outputs: %s", paste(missing_or_empty, collapse = ", ")))
}

scaled <- read.delim(paths[["scaled"]], stringsAsFactors = FALSE, check.names = FALSE)
projects <- read.delim(file.path(root, "data", "interim", "tcga_projects.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
qc <- read.delim(paths[["qc"]], stringsAsFactors = FALSE, check.names = FALSE)
bootstrap <- read.delim(paths[["bootstrap"]], stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(scaled) != 33) fail("Scaled Level 1 tree matrix does not include 33 TCGA projects")
if (!setequal(scaled$project_id, projects$project_id)) fail("Scaled Level 1 tree project IDs do not match expected TCGA projects")
if (!"major_cancer_group" %in% names(scaled)) fail("Scaled Level 1 tree matrix lacks major cancer group labels")
if (any(is.na(scaled$major_cancer_group) | scaled$major_cancer_group == "")) fail("Major cancer group labels contain missing values")

annotation_cols <- c("project_id", "project_code", "cancer_name", "disease_type", "primary_site", "broad_group", "major_cancer_group", "level2_group", "lineage_notes", "notes", "note", "source", "data_source")
feature_cols <- setdiff(names(scaled), annotation_cols)
if (length(feature_cols) == 0) fail("Scaled Level 1 tree matrix has no feature columns")
for (column in feature_cols) {
  values <- suppressWarnings(as.numeric(scaled[[column]]))
  if (any(is.na(values))) fail(sprintf("Scaled feature column is not numeric: %s", column))
}

read_distance <- function(path) {
  table <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"project_code" %in% names(table)) fail(sprintf("Distance matrix lacks project_code column: %s", path))
  matrix <- as.matrix(table[, setdiff(names(table), "project_code"), drop = FALSE])
  storage.mode(matrix) <- "numeric"
  rownames(matrix) <- table$project_code
  matrix
}

for (name in c("gower_distance", "correlation_distance")) {
  matrix <- read_distance(paths[[name]])
  if (nrow(matrix) != ncol(matrix)) fail(sprintf("%s is not square", name))
  if (nrow(matrix) != 33) fail(sprintf("%s does not include 33 projects", name))
  if (!identical(rownames(matrix), colnames(matrix))) fail(sprintf("%s row/column labels do not match", name))
  if (max(abs(matrix - t(matrix))) > 1e-8) fail(sprintf("%s is not symmetric", name))
  if (any(matrix < -1e-8)) fail(sprintf("%s contains negative distances", name))
}

required_qc_metrics <- c(
  "expected_tcga_projects",
  "included_tcga_projects",
  "raw_biologic_features",
  "features_retained_for_tree",
  "features_dropped",
  "imputed_values",
  "distance_methods_generated",
  "bootstrap_iterations_completed",
  "neighbor_joining_succeeded",
  "warnings_or_limitations"
)
if (!all(required_qc_metrics %in% qc$metric)) fail("Level 1 tree QC summary lacks required metrics")
if (nrow(bootstrap) == 0) fail("Bootstrap stability table is empty")
if (!all(c("project_code_a", "project_code_b", "co_clustering_frequency") %in% names(bootstrap))) {
  fail("Bootstrap stability table lacks required columns")
}

for (tree_name in c("gower_tree", "correlation_tree")) {
  text <- paste(readLines(paths[[tree_name]], warn = FALSE), collapse = "")
  if (nchar(text) == 0 || !grepl(";", text, fixed = TRUE)) fail(sprintf("Newick file is empty or malformed: %s", tree_name))
}

cat("Level 1 tree output checks passed\n")
