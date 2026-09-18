#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

required <- c(
  group_definitions = file.path(root, "results", "tables", "level2_group_definitions.tsv"),
  qc = file.path(root, "results", "tables", "level2_tree_qc_summary.tsv"),
  report = file.path(root, "results", "reports", "level2_interpretation_summary.md"),
  carcinoma_raw = file.path(root, "results", "tables", "level2_carcinoma_all_feature_matrix.tsv"),
  carcinoma_scaled = file.path(root, "results", "tables", "level2_carcinoma_all_feature_matrix_scaled.tsv"),
  carcinoma_gower_distance = file.path(root, "results", "tables", "level2_carcinoma_all_gower_distance_matrix.tsv"),
  carcinoma_correlation_distance = file.path(root, "results", "tables", "level2_carcinoma_all_correlation_distance_matrix.tsv"),
  carcinoma_gower_tree = file.path(root, "results", "trees", "level2_carcinoma_all_gower_hclust_tree.nwk"),
  carcinoma_correlation_tree = file.path(root, "results", "trees", "level2_carcinoma_all_correlation_hclust_tree.nwk"),
  carcinoma_gower_figure = file.path(root, "results", "figures", "level2_carcinoma_all_gower_tree.pdf"),
  carcinoma_correlation_figure = file.path(root, "results", "figures", "level2_carcinoma_all_correlation_tree.pdf"),
  carcinoma_heatmap = file.path(root, "results", "figures", "level2_carcinoma_all_trait_heatmap.pdf"),
  carcinoma_pca = file.path(root, "results", "figures", "level2_carcinoma_all_feature_pca.pdf")
)

missing_or_empty <- required[!file.exists(required) | file.info(required)$size == 0]
if (length(missing_or_empty) > 0) {
  fail(sprintf("Missing or empty Level 2 required outputs: %s", paste(missing_or_empty, collapse = ", ")))
}

defs <- read.delim(required[["group_definitions"]], stringsAsFactors = FALSE, check.names = FALSE)
qc <- read.delim(required[["qc"]], stringsAsFactors = FALSE, check.names = FALSE)
projects <- read.delim(file.path(root, "data", "interim", "tcga_projects.tsv"), stringsAsFactors = FALSE, check.names = FALSE)

required_def_cols <- c("group_name", "group_type", "included_project_codes", "included_project_ids", "n_projects", "n_features_used", "tree_attempted", "reason_if_not_attempted")
if (!all(required_def_cols %in% names(defs))) fail("Level 2 group definitions lack required columns")
if (!all(c("carcinoma_all", "pan_squamous", "gi_pancancreatobiliary", "kidney", "gynecologic_breast", "cns_glial", "melanocytic", "hematolymphoid", "sarcoma") %in% defs$group_name)) {
  fail("Level 2 group definitions lack expected groups")
}
valid_codes <- projects$project_code
included_codes <- unique(unlist(strsplit(paste(defs$included_project_codes, collapse = ";"), ";", fixed = TRUE)))
included_codes <- included_codes[nzchar(included_codes)]
if (!all(included_codes %in% valid_codes)) fail("Level 2 group definitions contain invalid TCGA project codes")

eligible <- defs[defs$tree_attempted == "TRUE" | defs$tree_attempted == TRUE, ]
if (nrow(eligible) == 0) fail("No eligible Level 2 groups generated trees")
for (group_name in eligible$group_name) {
  for (method in c("gower", "correlation")) {
    tree_path <- file.path(root, "results", "trees", sprintf("level2_%s_%s_hclust_tree.nwk", group_name, method))
    if (!file.exists(tree_path) || file.info(tree_path)$size == 0) fail(sprintf("Missing tree for eligible group %s (%s)", group_name, method))
    text <- paste(readLines(tree_path, warn = FALSE), collapse = "")
    if (!grepl(";", text, fixed = TRUE)) fail(sprintf("Newick tree appears malformed for %s (%s)", group_name, method))
  }
}

underpowered <- defs[defs$tree_attempted == "FALSE" | defs$tree_attempted == FALSE, ]
if (nrow(underpowered) == 0) fail("No underpowered Level 2 groups were recorded")
if (any(is.na(underpowered$reason_if_not_attempted) | underpowered$reason_if_not_attempted == "")) fail("Underpowered groups lack reasons")
if (!any(grepl("underpowered_project_count", underpowered$reason_if_not_attempted))) fail("Underpowered project-count reason was not recorded")

if (nrow(qc) == 0) fail("Level 2 QC summary is empty")
if (!any(qc$qc_section == "per_feature" & qc$metric == "drop_reason")) fail("Level 2 QC summary lacks per-feature drop reasons")

report <- paste(readLines(required[["report"]], warn = FALSE), collapse = "\n")
if (!grepl("molecular similarity", report, fixed = TRUE)) fail("Level 2 report lacks molecular similarity framing")
if (grepl("evolved from|ancestral to|descended from", report)) fail("Level 2 report contains overstrong evolutionary language")

cat("Level 2 tree output checks passed\n")
