#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

paths <- c(
  gower_assignments = file.path(root, "results", "tables", "level1_gower_cluster_assignments.tsv"),
  correlation_assignments = file.path(root, "results", "tables", "level1_correlation_cluster_assignments.tsv"),
  gower_summaries = file.path(root, "results", "tables", "level1_gower_cluster_summaries.tsv"),
  correlation_summaries = file.path(root, "results", "tables", "level1_correlation_cluster_summaries.tsv"),
  pca_loadings = file.path(root, "results", "tables", "level1_feature_pca_loadings.tsv"),
  variance = file.path(root, "results", "tables", "level1_feature_variance_summary.tsv"),
  association = file.path(root, "results", "tables", "level1_feature_cluster_association.tsv"),
  class_summary = file.path(root, "results", "tables", "level1_feature_class_sensitivity_summary.tsv"),
  class_assignments = file.path(root, "results", "tables", "level1_feature_class_cluster_assignments.tsv"),
  weighting_summary = file.path(root, "results", "tables", "level1_feature_weighting_sensitivity_summary.tsv"),
  weighting_assignments = file.path(root, "results", "tables", "level1_weighted_cluster_assignments.tsv"),
  stable_pairs = file.path(root, "results", "tables", "level1_stable_project_pairs.tsv"),
  bootstrap_project_stability = file.path(root, "results", "tables", "level1_bootstrap_project_stability.tsv"),
  report = file.path(root, "results", "reports", "level1_interpretation_summary.md"),
  trait_heatmap = file.path(root, "results", "figures", "level1_cluster_trait_summary_heatmap.pdf"),
  driver_barplot = file.path(root, "results", "figures", "level1_feature_driver_barplot.pdf"),
  bootstrap_heatmap = file.path(root, "results", "figures", "level1_bootstrap_coclustering_heatmap.pdf"),
  class_heatmap = file.path(root, "results", "figures", "level1_feature_class_sensitivity_heatmap.pdf"),
  weighting_mds = file.path(root, "results", "figures", "level1_weighting_sensitivity_pca_or_mds.pdf")
)

missing_or_empty <- paths[!file.exists(paths) | file.info(paths)$size == 0]
if (length(missing_or_empty) > 0) {
  fail(sprintf("Missing or empty Level 1 interpretation outputs: %s", paste(missing_or_empty, collapse = ", ")))
}

projects <- read.delim(file.path(root, "data", "interim", "tcga_projects.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
gower <- read.delim(paths[["gower_assignments"]], stringsAsFactors = FALSE, check.names = FALSE)
correlation <- read.delim(paths[["correlation_assignments"]], stringsAsFactors = FALSE, check.names = FALSE)
required_assignment_cols <- c("project_id", "project_code", "major_cancer_group", "disease_type", "primary_site", paste0("cluster_k", 3:8))

for (name in c("gower", "correlation")) {
  table <- get(name)
  if (!all(required_assignment_cols %in% names(table))) fail(sprintf("%s assignment table lacks required columns", name))
  if (nrow(table) != 33) fail(sprintf("%s assignment table does not include 33 projects", name))
  if (!setequal(table$project_id, projects$project_id)) fail(sprintf("%s assignment project IDs do not match expected TCGA projects", name))
  for (column in paste0("cluster_k", 3:8)) {
    if (any(is.na(table[[column]]))) fail(sprintf("%s contains missing %s values", name, column))
  }
}

for (path_name in c("gower_summaries", "correlation_summaries", "pca_loadings", "variance", "association", "class_summary", "class_assignments", "weighting_summary", "weighting_assignments", "stable_pairs", "bootstrap_project_stability")) {
  table <- read.delim(paths[[path_name]], stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(table) == 0) fail(sprintf("%s is empty", path_name))
}

association <- read.delim(paths[["association"]], stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("feature_name", "kruskal_wallis_statistic", "p_value", "fdr", "epsilon_squared") %in% names(association))) {
  fail("Feature association table lacks required columns")
}

class_summary <- read.delim(paths[["class_summary"]], stringsAsFactors = FALSE, check.names = FALSE)
if (!"adjusted_rand_index_vs_all_feature_gower_k6" %in% names(class_summary)) fail("Feature-class summary lacks ARI column")
if (any(is.na(class_summary$adjusted_rand_index_vs_all_feature_gower_k6))) fail("Feature-class ARI contains NA values")

weighting_summary <- read.delim(paths[["weighting_summary"]], stringsAsFactors = FALSE, check.names = FALSE)
if (!"adjusted_rand_index_vs_unweighted_gower_k6" %in% names(weighting_summary)) fail("Weighting summary lacks ARI column")
if (any(is.na(weighting_summary$adjusted_rand_index_vs_unweighted_gower_k6))) fail("Weighting ARI contains NA values")

report <- paste(readLines(paths[["report"]], warn = FALSE), collapse = "\n")
if (!grepl("molecular similarity", report, fixed = TRUE)) fail("Interpretation report lacks molecular similarity language")
if (grepl("descended from|ancestral to|evolved from", report)) fail("Interpretation report contains overstrong evolutionary language")

cat("Level 1 interpretation output checks passed\n")
