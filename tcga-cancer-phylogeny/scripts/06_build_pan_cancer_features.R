#!/usr/bin/env Rscript

source_common <- function() {
  frame_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE), error = function(err) NA_character_)
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  arg_file <- if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE) else NA_character_
  candidates <- unique(c(
    file.path(dirname(frame_file), "lib", "common.R"),
    file.path(dirname(arg_file), "lib", "common.R"),
    file.path(getwd(), "scripts", "lib", "common.R")
  ))
  candidates <- candidates[!is.na(candidates)]
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      source(candidate)
      return(invisible(candidate))
    }
  }
  stop("Could not locate scripts/lib/common.R", call. = FALSE)
}
source_common()

annotation_columns <- c("project_id", "project_code", "cancer_name", "broad_group", "level2_group", "lineage_notes")

require_feature_packages <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    abort("Package 'data.table' is required for Level 1 feature assembly")
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    abort("Package 'yaml' is required for Level 1 feature assembly")
  }
  suppressPackageStartupMessages(library(data.table))
}

read_tsv_required <- function(path, description) {
  require_file(path, description)
  data.table::fread(path, sep = "\t", quote = "\"", data.table = TRUE, showProgress = interactive())
}

write_tsv_table <- function(data, path) {
  ensure_dir(dirname(path))
  data.table::fwrite(data, path, sep = "\t", quote = FALSE, na = "NA")
  log_info(sprintf("Wrote %s", path))
  invisible(path)
}

read_feature_config <- function(path) {
  require_file(path, "feature-set configuration")
  cfg <- yaml::read_yaml(path)
  max_missing <- cfg$missingness$max_feature_missing_fraction
  if (is.null(max_missing) || !is.numeric(max_missing) || length(max_missing) != 1L) {
    max_missing <- 0.30
  }
  list(max_feature_missing_fraction = max_missing)
}

read_projects <- function(path) {
  projects <- read_tsv_required(path, "TCGA project metadata")
  required <- c("project_id", "project_code")
  missing <- setdiff(required, names(projects))
  if (length(missing) > 0) {
    abort(sprintf("Project metadata is missing required columns: %s", paste(missing, collapse = ", ")))
  }
  projects <- unique(projects[, ..required])
  data.table::setorder(projects, project_id)
  projects
}

read_cancer_group_map <- function(path) {
  require_file(path, "cancer group map")
  group_map <- data.table::fread(path, sep = ",", quote = "\"", data.table = TRUE, showProgress = FALSE)
  required <- annotation_columns
  missing <- setdiff(required, names(group_map))
  if (length(missing) > 0) {
    abort(sprintf("Cancer group map is missing required columns: %s", paste(missing, collapse = ", ")))
  }
  group_map <- unique(group_map[, ..required])
  duplicated_projects <- group_map$project_id[duplicated(group_map$project_id)]
  if (length(duplicated_projects) > 0) {
    abort(sprintf("Cancer group map has duplicate project IDs: %s", paste(unique(duplicated_projects), collapse = ", ")))
  }
  data.table::setorder(group_map, project_id)
  group_map
}

numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

summarize_tmb_by_project <- function(tmb_by_sample) {
  required <- c("sample_barcode", "project_id", "project_code", "nonsynonymous_count", "total_mutation_count")
  missing <- setdiff(required, names(tmb_by_sample))
  if (length(missing) > 0) {
    abort(sprintf("TMB sample table is missing required columns: %s", paste(missing, collapse = ", ")))
  }

  tmb_by_sample[, nonsynonymous_count := numeric_or_na(nonsynonymous_count)]
  tmb_by_sample[, total_mutation_count := numeric_or_na(total_mutation_count)]
  tmb_by_sample[, .(
    mutation_n_samples = data.table::uniqueN(sample_barcode),
    tmb_proxy_median_nonsynonymous_count = stats::median(nonsynonymous_count, na.rm = TRUE),
    tmb_proxy_mean_nonsynonymous_count = mean(nonsynonymous_count, na.rm = TRUE),
    tmb_proxy_median_total_mutation_count = stats::median(total_mutation_count, na.rm = TRUE),
    tmb_proxy_mean_total_mutation_count = mean(total_mutation_count, na.rm = TRUE)
  ), by = .(project_id, project_code)]
}

prepare_driver_prevalence <- function(driver_prevalence) {
  required <- c("project_id", "project_code")
  missing <- setdiff(required, names(driver_prevalence))
  if (length(missing) > 0) {
    abort(sprintf("Driver prevalence table is missing required columns: %s", paste(missing, collapse = ", ")))
  }

  gene_cols <- setdiff(names(driver_prevalence), required)
  for (column in gene_cols) {
    driver_prevalence[, (column) := numeric_or_na(get(column))]
  }
  data.table::setnames(driver_prevalence, gene_cols, paste0("mutation_prevalence_", gene_cols))
  driver_prevalence
}

prepare_copy_number_project_features <- function(copy_number_project) {
  required <- c("project_id", "project_code", "n_samples")
  missing <- setdiff(required, names(copy_number_project))
  if (length(missing) > 0) {
    abort(sprintf("Copy-number project table is missing required columns: %s", paste(missing, collapse = ", ")))
  }

  support <- copy_number_project[, .(
    project_id,
    project_code,
    copy_number_n_samples = numeric_or_na(n_samples)
  )]
  feature_cols <- setdiff(names(copy_number_project), c("project_id", "project_code", "n_samples", grep("^missingness_", names(copy_number_project), value = TRUE)))
  features <- copy_number_project[, c("project_id", "project_code", feature_cols), with = FALSE]
  for (column in feature_cols) {
    features[, (column) := numeric_or_na(get(column))]
  }
  list(features = features, support = support)
}

prepare_expression_project_features <- function(expression_project) {
  required <- c("project_id", "project_code", "n_samples")
  missing <- setdiff(required, names(expression_project))
  if (length(missing) > 0) {
    abort(sprintf("Expression project table is missing required columns: %s", paste(missing, collapse = ", ")))
  }

  support <- expression_project[, .(
    project_id,
    project_code,
    expression_n_samples = numeric_or_na(n_samples)
  )]
  feature_cols <- setdiff(names(expression_project), c("project_id", "project_code", "n_samples", grep("^missingness_", names(expression_project), value = TRUE)))
  features <- expression_project[, c("project_id", "project_code", feature_cols), with = FALSE]
  for (column in feature_cols) {
    features[, (column) := numeric_or_na(get(column))]
  }
  list(features = features, support = support)
}

left_join_project <- function(base, update) {
  merge(base, update, by = c("project_id", "project_code"), all.x = TRUE, sort = FALSE)
}

feature_category <- function(feature_name) {
  if (grepl("^tmb_proxy_", feature_name)) return("mutation_count_proxy")
  if (grepl("^mutation_prevalence_", feature_name)) return("driver_mutation_prevalence")
  if (feature_name %in% c(
    "median_purity",
    "median_ploidy",
    "median_aneuploidy_score",
    "median_fraction_genome_altered",
    "median_arm_gain_count",
    "median_arm_loss_count",
    "median_total_arm_alteration_count",
    "median_copy_number_complexity_score"
  )) return("copy_number_purity_ploidy")
  if (grepl("^median_.*score$|^median_lineage_or_tissue_PC", feature_name)) return("expression")
  "other"
}

filter_features_by_missingness <- function(matrix, max_missing_fraction) {
  feature_cols <- setdiff(names(matrix), annotation_columns)
  missingness <- vapply(feature_cols, function(column) mean(is.na(matrix[[column]])), numeric(1))
  retained <- names(missingness)[missingness <= max_missing_fraction]
  dropped <- names(missingness)[missingness > max_missing_fraction]

  filtered <- matrix[, c(annotation_columns, retained), with = FALSE]
  qc <- data.table::data.table(
    feature_name = names(missingness),
    feature_category = vapply(names(missingness), feature_category, character(1)),
    n_missing_projects = vapply(names(missingness), function(column) sum(is.na(matrix[[column]])), integer(1)),
    missing_fraction = as.numeric(missingness),
    retained = names(missingness) %in% retained,
    drop_reason = ifelse(names(missingness) %in% dropped, sprintf("missing_fraction_gt_%.2f", max_missing_fraction), "")
  )
  data.table::setorder(qc, retained, feature_category, feature_name)
  list(matrix = filtered, feature_qc = qc)
}

build_level1_feature_matrix <- function(projects, group_map, tmb_by_sample, driver_prevalence, copy_number_project, expression_project, max_missing_fraction = 0.30) {
  missing_group_projects <- setdiff(projects$project_id, group_map$project_id)
  if (length(missing_group_projects) > 0) {
    abort(sprintf("Cancer group map is missing project IDs: %s", paste(missing_group_projects, collapse = ", ")))
  }

  base <- left_join_project(projects, group_map)
  tmb_project <- summarize_tmb_by_project(copy(tmb_by_sample))
  driver_features <- prepare_driver_prevalence(copy(driver_prevalence))
  copy_number <- prepare_copy_number_project_features(copy(copy_number_project))
  expression <- prepare_expression_project_features(copy(expression_project))

  matrix <- base
  matrix <- left_join_project(matrix, tmb_project[, setdiff(names(tmb_project), "mutation_n_samples"), with = FALSE])
  matrix <- left_join_project(matrix, driver_features)
  matrix <- left_join_project(matrix, copy_number$features)
  matrix <- left_join_project(matrix, expression$features)

  for (column in setdiff(names(matrix), annotation_columns)) {
    matrix[, (column) := numeric_or_na(get(column))]
  }

  filtered <- filter_features_by_missingness(matrix, max_missing_fraction)

  support <- base[, .(project_id, project_code, cancer_name, broad_group, level2_group)]
  support <- left_join_project(support, tmb_project[, .(project_id, project_code, mutation_n_samples)])
  support <- left_join_project(support, copy_number$support)
  support <- left_join_project(support, expression$support)

  list(
    feature_matrix = filtered$matrix,
    feature_qc = filtered$feature_qc,
    sample_support = support,
    unfiltered_feature_matrix = matrix
  )
}

write_level1_qc_summary <- function(features, projects, max_missing_fraction, output_path) {
  retained_count <- features$feature_qc[retained == TRUE, .N]
  dropped_count <- features$feature_qc[retained == FALSE, .N]
  matrix <- features$feature_matrix
  missing_projects <- setdiff(projects$project_id, matrix$project_id)

  global <- data.table::data.table(
    qc_section = "global",
    metric = c(
      "project_rows",
      "expected_tcga_projects",
      "expected_tcga_projects_missing",
      "missing_tcga_projects",
      "retained_feature_columns",
      "dropped_feature_columns",
      "max_feature_missing_fraction",
      "scaling_status",
      "support_columns_location"
    ),
    value = c(
      as.character(nrow(matrix)),
      as.character(nrow(projects)),
      as.character(length(missing_projects)),
      ifelse(length(missing_projects) == 0, "none", paste(missing_projects, collapse = ";")),
      as.character(retained_count),
      as.character(dropped_count),
      as.character(max_missing_fraction),
      "raw_unscaled_features_scale_before_distance_calculation",
      "results/tables/level1_feature_support.tsv"
    )
  )

  feature_qc <- copy(features$feature_qc)
  feature_qc[, qc_section := "per_feature"]
  feature_measure_cols <- setdiff(names(feature_qc), c("qc_section", "feature_name", "feature_category"))
  for (column in feature_measure_cols) {
    feature_qc[, (column) := as.character(get(column))]
  }
  feature_qc_long <- data.table::melt(
    feature_qc,
    id.vars = c("qc_section", "feature_name", "feature_category"),
    variable.name = "metric",
    value.name = "value"
  )
  feature_qc_long[, value := as.character(value)]

  qc <- data.table::rbindlist(
    list(
      global[, .(qc_section, feature_name = "ALL", feature_category = "ALL", metric, value)],
      feature_qc_long[, .(qc_section, feature_name, feature_category, metric, value)]
    ),
    use.names = TRUE,
    fill = TRUE
  )
  data.table::setorder(qc, qc_section, feature_category, feature_name, metric)
  write_tsv_table(qc, output_path)
  invisible(qc)
}

validate_level1_outputs <- function(feature_matrix, projects) {
  if (nrow(feature_matrix) == 0) {
    abort("Level 1 feature matrix has zero rows")
  }
  if (anyDuplicated(feature_matrix$project_id)) {
    abort("Level 1 feature matrix has duplicate project IDs")
  }
  missing_projects <- setdiff(projects$project_id, feature_matrix$project_id)
  if (length(missing_projects) > 0) {
    abort(sprintf("Level 1 feature matrix is missing expected projects: %s", paste(missing_projects, collapse = ", ")))
  }
  feature_cols <- setdiff(names(feature_matrix), annotation_columns)
  if (length(feature_cols) == 0) {
    abort("Level 1 feature matrix has no retained feature columns")
  }
  invisible(TRUE)
}

main <- function() {
  require_feature_packages()
  root <- find_project_root()

  config <- read_feature_config(project_path("config", "feature_sets.yaml", root = root))
  projects <- read_projects(project_path("data", "interim", "tcga_projects.tsv", root = root))
  group_map <- read_cancer_group_map(project_path("config", "cancer_group_map.csv", root = root))
  tmb_by_sample <- read_tsv_required(project_path("data", "processed", "features", "tmb_by_sample.tsv", root = root), "sample-level mutation-count proxy table")
  driver_prevalence <- read_tsv_required(project_path("data", "processed", "features", "mutation_prevalence_by_project_driver_only.tsv", root = root), "driver mutation prevalence table")
  copy_number_project <- read_tsv_required(project_path("data", "processed", "features", "aneuploidy_by_project.tsv", root = root), "project-level copy-number table")
  expression_project <- read_tsv_required(project_path("data", "processed", "features", "expression_pathway_scores_by_project.tsv", root = root), "project-level expression score table")

  features <- build_level1_feature_matrix(
    projects,
    group_map,
    tmb_by_sample,
    driver_prevalence,
    copy_number_project,
    expression_project,
    max_missing_fraction = config$max_feature_missing_fraction
  )
  validate_level1_outputs(features$feature_matrix, projects)

  matrix_output <- project_path("results", "tables", "level1_feature_matrix.tsv", root = root)
  support_output <- project_path("results", "tables", "level1_feature_support.tsv", root = root)
  qc_output <- project_path("results", "tables", "level1_feature_matrix_qc_summary.tsv", root = root)

  write_tsv_table(features$feature_matrix, matrix_output)
  write_tsv_table(features$sample_support, support_output)
  qc <- write_level1_qc_summary(features, projects, config$max_feature_missing_fraction, qc_output)

  log_info(sprintf("Level 1 feature matrix complete: %d projects, %d retained features", nrow(features$feature_matrix), length(setdiff(names(features$feature_matrix), annotation_columns))))
  log_info(sprintf("Dropped %d features for missingness > %.2f", features$feature_qc[retained == FALSE, .N], config$max_feature_missing_fraction))
  invisible(qc)
}

if (identical(environment(), globalenv())) {
  main()
}
