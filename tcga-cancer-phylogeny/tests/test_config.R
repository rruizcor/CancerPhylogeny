#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)

fail <- function(message) {
  stop(message, call. = FALSE)
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  fail("Package 'yaml' is required for config tests")
}

projects_cfg <- yaml::read_yaml(file.path(root, "config", "tcga_projects.yaml"))
group_map <- read.csv(file.path(root, "config", "cancer_group_map.csv"), stringsAsFactors = FALSE)
driver_warnings <- character()
drivers <- withCallingHandlers(
  read.csv(
    file.path(root, "config", "driver_genes.csv"),
    stringsAsFactors = FALSE,
    quote = "\"",
    fill = FALSE
  ),
  warning = function(warn) {
    driver_warnings <<- c(driver_warnings, conditionMessage(warn))
    invokeRestart("muffleWarning")
  }
)
features_cfg <- yaml::read_yaml(file.path(root, "config", "feature_sets.yaml"))

project_codes <- vapply(projects_cfg$projects, `[[`, character(1), "code")
project_ids <- vapply(projects_cfg$projects, `[[`, character(1), "project_id")

expected_codes <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
  "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
  "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

if (!setequal(project_codes, expected_codes)) {
  fail("Configured TCGA project codes do not match the expected project list")
}

if (anyDuplicated(project_codes) > 0 || anyDuplicated(project_ids) > 0) {
  fail("Duplicate project codes or project IDs found")
}

if (!all(project_codes %in% group_map$project_code)) {
  fail("Every configured project must appear in config/cancer_group_map.csv")
}

required_groups <- c(
  "Carcinoma",
  "Sarcoma",
  "Hematolymphoid",
  "Melanocytic",
  "CNS/glial",
  "Mesothelial",
  "Germ cell",
  "Neuroendocrine/adrenal/paraganglioma",
  "Other or mixed"
)

if (!all(required_groups %in% group_map$broad_group)) {
  fail("Cancer group map is missing one or more required broad groups")
}

if (length(driver_warnings) > 0) {
  fail(sprintf("Driver gene table emitted CSV parsing warnings: %s", paste(unique(driver_warnings), collapse = " | ")))
}

expected_driver_columns <- c("gene", "feature_group", "rationale")
if (!identical(names(drivers), expected_driver_columns)) {
  fail(sprintf(
    "Driver gene table must have exactly these columns: %s",
    paste(expected_driver_columns, collapse = ", ")
  ))
}

if (ncol(drivers) != 3) {
  fail("Driver gene table must have exactly three columns")
}

apc <- drivers[drivers$gene == "APC", , drop = FALSE]
if (nrow(apc) != 1 || !identical(apc$rationale, "WNT pathway driver, especially colorectal lineage.")) {
  fail("APC rationale was not preserved as one CSV field")
}

if (any(is.na(drivers$gene) | trimws(drivers$gene) == "")) {
  fail("Driver gene table contains missing or blank gene symbols")
}

if (anyDuplicated(drivers$gene) > 0) {
  fail("Driver gene table contains duplicate gene symbols")
}

expected_driver_genes <- c("TP53", "KRAS", "BRAF", "IDH1", "IDH2", "PIK3CA", "PTEN", "APC", "NF1", "BAP1", "RB1")
missing_expected_driver_genes <- setdiff(expected_driver_genes, drivers$gene)
if (length(missing_expected_driver_genes) > 0) {
  fail(sprintf("Driver gene table is missing expected genes: %s", paste(missing_expected_driver_genes, collapse = ", ")))
}

required_feature_categories <- c("mutation", "copy_number", "expression", "histogenetic_group")
if (!all(required_feature_categories %in% names(features_cfg$feature_categories))) {
  fail("Feature configuration is missing required feature categories")
}

cat("Configuration sanity checks passed\n")
