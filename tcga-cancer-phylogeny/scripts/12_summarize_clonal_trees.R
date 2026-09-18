#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

summarize_sample_clusters <- function(cluster_dir) {
  not_implemented(
    "Clonal cluster summarization",
    "read PyClone-VI cluster outputs and summarize clone counts, cellular prevalence ranges, and QC warnings"
  )
}

summarize_candidate_trees <- function(tree_dir) {
  not_implemented(
    "Candidate clonal tree summarization",
    "summarize PhyloWGS or other candidate tree structures when data support them"
  )
}

main <- function() {
  root <- find_project_root()
  output_path <- project_path("results", "tables", "clonal_summary_by_sample.tsv", root = root)
  log_warn("Single-bulk TCGA samples often support clonal clustering better than confident branching phylogenies")
  log_info(sprintf("Planned clonal summary output: %s", output_path))
  summarize_sample_clusters(project_path("results", "clonal", "sample_cluster_assignments", root = root))
}

if (identical(environment(), globalenv())) {
  main()
}
