#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

plot_level1_tree <- function(tree_path, feature_matrix_path, group_map_path, output_path) {
  not_implemented(
    "Level 1 tree plotting",
    "visualize tree tips by broad cancer group and add trait heatmap or sidebars for mapped features"
  )
}

plot_level2_trees <- function(tree_dir, output_dir) {
  not_implemented(
    "Level 2 tree plotting",
    "generate one figure per group with cluster stability annotations where available"
  )
}

plot_clonal_pilot_summary <- function(summary_path, output_path) {
  not_implemented(
    "Clonal pilot plotting",
    "summarize sample eligibility, clone counts, prevalence estimates, and insufficient-data warnings"
  )
}

main <- function() {
  root <- find_project_root()
  log_info("Planned figure generation for Level 1, Level 2, and clonal pilot summaries")
  plot_level1_tree(
    project_path("results", "trees", "level1_pan_cancer_gower_tree.nwk", root = root),
    project_path("results", "tables", "level1_feature_matrix.tsv", root = root),
    project_path("config", "cancer_group_map.csv", root = root),
    project_path("results", "figures", "level1_pan_cancer_tree.pdf", root = root)
  )
}

if (identical(environment(), globalenv())) {
  main()
}
