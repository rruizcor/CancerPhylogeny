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

tree_annotation_columns <- c(
  "project_id",
  "project_code",
  "cancer_name",
  "disease_type",
  "primary_site",
  "broad_group",
  "major_cancer_group",
  "level2_group",
  "lineage_notes",
  "notes",
  "note",
  "source",
  "data_source"
)

require_tree_packages <- function() {
  required <- c("data.table", "cluster", "ape", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    abort(sprintf("Missing R package(s) required for Level 1 tree construction: %s", paste(missing, collapse = ", ")))
  }
  suppressPackageStartupMessages({
    library(data.table)
  })
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

parse_tree_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  iterations <- 100L
  seed <- 1L
  max_missing_fraction <- 0.30
  for (arg in args) {
    if (grepl("^--bootstrap-iterations=", arg)) {
      iterations <- as.integer(sub("^--bootstrap-iterations=", "", arg))
    }
    if (grepl("^--seed=", arg)) {
      seed <- as.integer(sub("^--seed=", "", arg))
    }
    if (grepl("^--max-missing-fraction=", arg)) {
      max_missing_fraction <- as.numeric(sub("^--max-missing-fraction=", "", arg))
    }
  }
  if (is.na(iterations) || iterations < 0L) abort("--bootstrap-iterations must be a non-negative integer")
  if (is.na(seed)) abort("--seed must be an integer")
  if (is.na(max_missing_fraction) || max_missing_fraction < 0 || max_missing_fraction > 1) {
    abort("--max-missing-fraction must be between 0 and 1")
  }
  list(bootstrap_iterations = iterations, seed = seed, max_missing_fraction = max_missing_fraction)
}

standardize_tree_annotations <- function(feature_matrix) {
  if (!"major_cancer_group" %in% names(feature_matrix) && "broad_group" %in% names(feature_matrix)) {
    feature_matrix[, major_cancer_group := broad_group]
  }
  if (!"major_cancer_group" %in% names(feature_matrix)) {
    feature_matrix[, major_cancer_group := "Unannotated"]
  }
  feature_matrix[is.na(major_cancer_group) | major_cancer_group == "", major_cancer_group := "Unannotated"]
  feature_matrix
}

numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

select_biologic_feature_columns <- function(feature_matrix, annotation_columns = tree_annotation_columns) {
  candidates <- setdiff(names(feature_matrix), annotation_columns)
  candidates <- candidates[!grepl("(^|_)n_samples$|^missingness_|_support$|^support_", candidates)]
  numeric_candidates <- candidates[vapply(candidates, function(column) {
    values <- feature_matrix[[column]]
    converted <- numeric_or_na(values)
    all(is.na(values) | !is.na(converted))
  }, logical(1))]
  numeric_candidates
}

prepare_tree_features <- function(feature_matrix, max_missing_fraction = 0.30, annotation_columns = tree_annotation_columns) {
  matrix <- data.table::copy(feature_matrix)
  matrix <- standardize_tree_annotations(matrix)
  feature_cols <- select_biologic_feature_columns(matrix, annotation_columns)
  if (length(feature_cols) == 0) {
    abort("No numeric biologic features were available for Level 1 tree construction")
  }

  for (column in feature_cols) {
    matrix[, (column) := numeric_or_na(get(column))]
  }

  missing_fraction <- vapply(feature_cols, function(column) mean(is.na(matrix[[column]])), numeric(1))
  dropped_missing <- names(missing_fraction)[missing_fraction > max_missing_fraction]
  retained <- setdiff(feature_cols, dropped_missing)
  if (length(retained) == 0) {
    abort("All biologic features were dropped by the tree missingness filter")
  }

  medians <- vapply(retained, function(column) {
    stats::median(matrix[[column]], na.rm = TRUE)
  }, numeric(1))
  all_missing_medians <- names(medians)[is.na(medians)]
  if (length(all_missing_medians) > 0) {
    retained <- setdiff(retained, all_missing_medians)
    dropped_missing <- unique(c(dropped_missing, all_missing_medians))
    medians <- medians[retained]
  }
  if (length(retained) == 0) {
    abort("No biologic features remained after removing all-missing features")
  }

  imputed_counts <- vapply(retained, function(column) sum(is.na(matrix[[column]])), integer(1))
  for (column in retained) {
    if (imputed_counts[[column]] > 0) {
      matrix[is.na(get(column)), (column) := medians[[column]]]
    }
  }

  zero_variance <- retained[vapply(retained, function(column) {
    stats::sd(matrix[[column]], na.rm = TRUE) == 0
  }, logical(1))]
  retained <- setdiff(retained, zero_variance)
  if (length(retained) == 0) {
    abort("No biologic features remained after removing zero-variance features")
  }

  annotation_present <- intersect(annotation_columns, names(matrix))
  annotations <- matrix[, ..annotation_present]
  feature_matrix_retained <- as.matrix(matrix[, ..retained])
  rownames(feature_matrix_retained) <- matrix$project_code

  feature_qc <- data.table::data.table(
    feature_name = feature_cols,
    missing_fraction = as.numeric(missing_fraction[feature_cols]),
    n_missing_projects = vapply(feature_cols, function(column) sum(is.na(feature_matrix[[column]])), integer(1)),
    imputed_value = NA_real_,
    n_imputed_values = 0L,
    retained_for_tree = feature_cols %in% retained,
    drop_reason = ""
  )
  feature_qc[feature_name %in% names(medians), imputed_value := as.numeric(medians[feature_name])]
  feature_qc[feature_name %in% names(imputed_counts), n_imputed_values := as.integer(imputed_counts[feature_name])]
  feature_qc[feature_name %in% dropped_missing, drop_reason := sprintf("missing_fraction_gt_%.2f_or_all_missing", max_missing_fraction)]
  feature_qc[feature_name %in% zero_variance, drop_reason := "zero_variance_after_imputation"]

  list(
    annotations = annotations,
    raw_feature_matrix = feature_matrix_retained,
    feature_columns_raw = feature_cols,
    feature_columns_retained = retained,
    feature_qc = feature_qc
  )
}

scale_tree_feature_matrix <- function(prepared) {
  scaled <- scale(prepared$raw_feature_matrix)
  scaled[is.na(scaled)] <- 0
  scaled <- as.matrix(scaled)
  rownames(scaled) <- rownames(prepared$raw_feature_matrix)
  list(
    annotations = prepared$annotations,
    scaled_matrix = scaled
  )
}

write_scaled_feature_matrix <- function(scaled, path) {
  annotations <- data.table::copy(scaled$annotations)
  features <- data.table::as.data.table(scaled$scaled_matrix)
  output <- cbind(annotations, features)
  write_tsv_table(output, path)
}

validate_distance_matrix <- function(distance_matrix, label) {
  if (!is.matrix(distance_matrix)) abort(sprintf("%s distance output is not a matrix", label))
  if (nrow(distance_matrix) != ncol(distance_matrix)) abort(sprintf("%s distance matrix is not square", label))
  if (any(is.na(distance_matrix))) abort(sprintf("%s distance matrix contains NA values", label))
  if (max(abs(distance_matrix - t(distance_matrix))) > 1e-8) abort(sprintf("%s distance matrix is not symmetric", label))
  if (any(distance_matrix < -1e-8)) abort(sprintf("%s distance matrix contains negative distances", label))
  diag(distance_matrix) <- 0
  distance_matrix[distance_matrix < 0] <- 0
  distance_matrix
}

compute_gower_distance_matrix <- function(scaled_matrix) {
  distance <- cluster::daisy(as.data.frame(scaled_matrix), metric = "gower")
  matrix <- as.matrix(distance)
  rownames(matrix) <- rownames(scaled_matrix)
  colnames(matrix) <- rownames(scaled_matrix)
  validate_distance_matrix(matrix, "Gower")
}

compute_correlation_distance_matrix <- function(scaled_matrix) {
  correlation <- suppressWarnings(stats::cor(t(scaled_matrix), use = "pairwise.complete.obs", method = "pearson"))
  correlation[is.na(correlation)] <- 0
  correlation[correlation > 1] <- 1
  correlation[correlation < -1] <- -1
  distance <- 1 - correlation
  diag(distance) <- 0
  distance <- (distance + t(distance)) / 2
  distance[distance < 0 & distance > -1e-8] <- 0
  rownames(distance) <- rownames(scaled_matrix)
  colnames(distance) <- rownames(scaled_matrix)
  validate_distance_matrix(distance, "Correlation")
}

write_distance_matrix <- function(distance_matrix, path) {
  output <- data.table::as.data.table(distance_matrix, keep.rownames = "project_code")
  write_tsv_table(output, path)
}

build_hclust_tree <- function(distance_matrix, method = "average") {
  stats::hclust(stats::as.dist(distance_matrix), method = method)
}

write_hclust_newick <- function(hclust_object, path) {
  ensure_dir(dirname(path))
  phylo <- ape::as.phylo(hclust_object)
  ape::write.tree(phylo, file = path)
  log_info(sprintf("Wrote %s", path))
  invisible(phylo)
}

attempt_neighbor_joining <- function(distance_matrix, path) {
  ensure_dir(dirname(path))
  result <- tryCatch({
    phylo <- ape::nj(stats::as.dist(distance_matrix))
    ape::write.tree(phylo, file = path)
    list(success = TRUE, message = "neighbor_joining_succeeded", phylo = phylo)
  }, error = function(err) {
    list(success = FALSE, message = conditionMessage(err), phylo = NULL)
  })
  if (isTRUE(result$success)) {
    log_info(sprintf("Wrote %s", path))
  } else {
    log_warn(sprintf("Neighbor-joining tree skipped: %s", result$message))
  }
  result
}

group_palette <- function(groups) {
  groups <- sort(unique(groups))
  n <- length(groups)
  if (n <= 8) {
    colors <- RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)]
  } else if (n <= 12) {
    colors <- RColorBrewer::brewer.pal(n, "Paired")
  } else {
    colors <- grDevices::rainbow(n)
  }
  stats::setNames(colors, groups)
}

plot_tree_pdf <- function(phylo, annotations, path, title) {
  ensure_dir(dirname(path))
  groups <- annotations$major_cancer_group
  names(groups) <- annotations$project_code
  palette <- group_palette(groups)
  tip_colors <- palette[groups[phylo$tip.label]]
  grDevices::pdf(path, width = 11, height = 8.5)
  op <- par(mar = c(2, 2, 4, 1))
  on.exit({
    par(op)
    grDevices::dev.off()
  }, add = TRUE)
  ape::plot.phylo(
    phylo,
    type = "phylogram",
    cex = 0.72,
    tip.color = tip_colors,
    main = title,
    no.margin = FALSE
  )
  legend(
    "topleft",
    legend = names(palette),
    col = palette,
    pch = 19,
    cex = 0.58,
    bty = "n",
    title = "Major cancer group"
  )
  invisible(path)
}

selected_heatmap_traits <- function(scaled_matrix) {
  preferred <- c(
    "tmb_proxy_median_nonsynonymous_count",
    "tmb_proxy_median_total_mutation_count",
    "mutation_prevalence_TP53",
    "mutation_prevalence_KRAS",
    "mutation_prevalence_BRAF",
    "mutation_prevalence_IDH1",
    "mutation_prevalence_IDH2",
    "median_aneuploidy_score",
    "median_ploidy",
    "median_purity",
    "median_immune_inflammatory_score",
    "median_proliferation_score",
    "median_EMT_score"
  )
  intersect(preferred, colnames(scaled_matrix))
}

plot_trait_heatmap_pdf <- function(scaled_matrix, hclust_object, annotations, path) {
  ensure_dir(dirname(path))
  traits <- selected_heatmap_traits(scaled_matrix)
  if (length(traits) == 0) {
    log_warn("Trait heatmap skipped because none of the requested interpretable traits were present")
    return(invisible(FALSE))
  }
  heatmap_matrix <- scaled_matrix[, traits, drop = FALSE]
  row_annotation <- data.frame(
    major_cancer_group = annotations$major_cancer_group,
    row.names = annotations$project_code,
    check.names = FALSE
  )
  palette <- group_palette(annotations$major_cancer_group)
  annotation_colors <- list(major_cancer_group = palette)
  pheatmap::pheatmap(
    heatmap_matrix,
    filename = path,
    width = 11,
    height = 8.5,
    cluster_rows = hclust_object,
    cluster_cols = TRUE,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    main = "Level 1 pan-cancer molecular trait tree with selected traits",
    fontsize_row = 7,
    fontsize_col = 8,
    border_color = NA
  )
  invisible(TRUE)
}

plot_pca_pdf <- function(scaled_matrix, annotations, path) {
  ensure_dir(dirname(path))
  pca <- stats::prcomp(scaled_matrix, center = FALSE, scale. = FALSE)
  variance <- (pca$sdev^2) / sum(pca$sdev^2)
  scores <- data.table::data.table(
    project_code = rownames(scaled_matrix),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  plot_data <- merge(scores, annotations, by = "project_code", all.x = TRUE, sort = FALSE)
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC1, y = PC2, color = major_cancer_group, label = project_code)) +
    ggplot2::geom_point(size = 2.8) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::labs(
      title = "Level 1 pan-cancer molecular feature PCA",
      subtitle = "Scaled project-level molecular traits; colors show major cancer groups",
      x = sprintf("PC1 (%.1f%%)", 100 * variance[[1]]),
      y = sprintf("PC2 (%.1f%%)", 100 * variance[[2]]),
      color = "Major cancer group"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )
  ggplot2::ggsave(path, plot, width = 11, height = 8.5, device = "pdf")
  invisible(path)
}

bootstrap_cluster_stability <- function(scaled_matrix, iterations = 100L, seed = 1L, k = NULL) {
  projects <- rownames(scaled_matrix)
  n_projects <- length(projects)
  if (is.null(k)) {
    k <- max(2L, min(n_projects - 1L, round(sqrt(n_projects))))
  }
  pair_grid <- utils::combn(projects, 2)
  counts <- integer(ncol(pair_grid))
  if (iterations == 0L) {
    return(data.table::data.table(
      project_code_a = pair_grid[1, ],
      project_code_b = pair_grid[2, ],
      bootstrap_iterations = 0L,
      cutree_k = k,
      n_same_cluster = 0L,
      co_clustering_frequency = NA_real_
    ))
  }
  set.seed(seed)
  for (iteration in seq_len(iterations)) {
    columns <- sample(seq_len(ncol(scaled_matrix)), size = ncol(scaled_matrix), replace = TRUE)
    boot_matrix <- scaled_matrix[, columns, drop = FALSE]
    colnames(boot_matrix) <- make.unique(colnames(boot_matrix))
    distance <- compute_gower_distance_matrix(boot_matrix)
    clusters <- stats::cutree(build_hclust_tree(distance), k = k)
    same <- clusters[pair_grid[1, ]] == clusters[pair_grid[2, ]]
    counts <- counts + as.integer(same)
  }
  stability <- data.table::data.table(
    project_code_a = pair_grid[1, ],
    project_code_b = pair_grid[2, ],
    bootstrap_iterations = iterations,
    cutree_k = k,
    n_same_cluster = counts,
    co_clustering_frequency = counts / iterations
  )
  data.table::setorder(stability, -co_clustering_frequency, project_code_a, project_code_b)
  stability
}

write_tree_qc_summary <- function(qc, output_path) {
  global_qc <- qc[setdiff(names(qc), "feature_qc")]
  global <- data.table::data.table(
    qc_section = "global",
    feature_name = "ALL",
    metric = names(global_qc),
    value = vapply(global_qc, function(x) paste(as.character(x), collapse = ";"), character(1))
  )
  feature_qc <- data.table::copy(qc$feature_qc)
  feature_qc[, qc_section := "per_feature"]
  measure_cols <- setdiff(names(feature_qc), c("qc_section", "feature_name"))
  for (column in measure_cols) {
    feature_qc[, (column) := as.character(get(column))]
  }
  feature_qc_long <- data.table::melt(
    feature_qc,
    id.vars = c("qc_section", "feature_name"),
    variable.name = "metric",
    value.name = "value"
  )
  summary <- data.table::rbindlist(
    list(global, feature_qc_long[, .(qc_section, feature_name, metric, value = as.character(value))]),
    use.names = TRUE,
    fill = TRUE
  )
  data.table::setorder(summary, qc_section, feature_name, metric)
  write_tsv_table(summary, output_path)
  invisible(summary)
}

build_level1_tree_outputs <- function(feature_matrix, projects, max_missing_fraction = 0.30, bootstrap_iterations = 100L, seed = 1L, root = find_project_root()) {
  feature_matrix <- standardize_tree_annotations(data.table::copy(feature_matrix))
  missing_projects <- setdiff(projects$project_id, feature_matrix$project_id)
  if (anyDuplicated(feature_matrix$project_id)) abort("Level 1 feature matrix has duplicate project IDs")
  if (length(missing_projects) > 0) {
    abort(sprintf("Level 1 tree input is missing expected TCGA projects: %s", paste(missing_projects, collapse = ", ")))
  }

  prepared <- prepare_tree_features(feature_matrix, max_missing_fraction = max_missing_fraction)
  scaled <- scale_tree_feature_matrix(prepared)

  table_dir <- project_path("results", "tables", root = root)
  tree_dir <- project_path("results", "trees", root = root)
  figure_dir <- project_path("results", "figures", root = root)

  scaled_path <- file.path(table_dir, "level1_feature_matrix_scaled.tsv")
  gower_distance_path <- file.path(table_dir, "level1_gower_distance_matrix.tsv")
  correlation_distance_path <- file.path(table_dir, "level1_correlation_distance_matrix.tsv")
  bootstrap_path <- file.path(table_dir, "level1_bootstrap_cluster_stability.tsv")
  qc_path <- file.path(table_dir, "level1_tree_qc_summary.tsv")
  gower_tree_path <- file.path(tree_dir, "level1_pan_cancer_gower_hclust_tree.nwk")
  correlation_tree_path <- file.path(tree_dir, "level1_pan_cancer_correlation_hclust_tree.nwk")
  nj_tree_path <- file.path(tree_dir, "level1_pan_cancer_neighbor_joining_tree.nwk")
  gower_figure_path <- file.path(figure_dir, "level1_pan_cancer_gower_tree.pdf")
  correlation_figure_path <- file.path(figure_dir, "level1_pan_cancer_correlation_tree.pdf")
  heatmap_figure_path <- file.path(figure_dir, "level1_pan_cancer_tree_with_trait_heatmap.pdf")
  pca_figure_path <- file.path(figure_dir, "level1_pan_cancer_feature_pca.pdf")

  write_scaled_feature_matrix(scaled, scaled_path)

  gower_distance <- compute_gower_distance_matrix(scaled$scaled_matrix)
  correlation_distance <- compute_correlation_distance_matrix(scaled$scaled_matrix)
  write_distance_matrix(gower_distance, gower_distance_path)
  write_distance_matrix(correlation_distance, correlation_distance_path)

  gower_hclust <- build_hclust_tree(gower_distance, method = "average")
  correlation_hclust <- build_hclust_tree(correlation_distance, method = "average")
  gower_phylo <- write_hclust_newick(gower_hclust, gower_tree_path)
  correlation_phylo <- write_hclust_newick(correlation_hclust, correlation_tree_path)
  nj_result <- attempt_neighbor_joining(gower_distance, nj_tree_path)

  plot_tree_pdf(
    gower_phylo,
    scaled$annotations,
    gower_figure_path,
    "Level 1 pan-cancer molecular similarity dendrogram (Gower distance)"
  )
  plot_tree_pdf(
    correlation_phylo,
    scaled$annotations,
    correlation_figure_path,
    "Level 1 pan-cancer molecular similarity dendrogram (correlation distance)"
  )
  plot_trait_heatmap_pdf(scaled$scaled_matrix, gower_hclust, scaled$annotations, heatmap_figure_path)
  plot_pca_pdf(scaled$scaled_matrix, scaled$annotations, pca_figure_path)

  stability <- bootstrap_cluster_stability(scaled$scaled_matrix, iterations = bootstrap_iterations, seed = seed)
  write_tsv_table(stability, bootstrap_path)

  dropped_features <- prepared$feature_qc[retained_for_tree == FALSE, feature_name]
  warnings <- character()
  if (length(dropped_features) > 0) {
    warnings <- c(warnings, sprintf("Dropped features: %s", paste(dropped_features, collapse = ",")))
  }
  imputed_total <- sum(prepared$feature_qc$n_imputed_values)
  if (imputed_total > 0) {
    warnings <- c(warnings, sprintf("Median-imputed %d feature values", imputed_total))
  }
  if (!isTRUE(nj_result$success)) {
    warnings <- c(warnings, sprintf("Neighbor joining skipped: %s", nj_result$message))
  }
  warnings <- c(warnings, "This is a pan-cancer molecular similarity dendrogram / molecular trait tree, not a literal species-like cancer phylogeny")

  tree_files <- c(gower_tree_path, correlation_tree_path)
  if (isTRUE(nj_result$success)) tree_files <- c(tree_files, nj_tree_path)
  qc_values <- list(
    expected_tcga_projects = nrow(projects),
    included_tcga_projects = nrow(feature_matrix),
    missing_expected_projects = ifelse(length(missing_projects) == 0, "none", paste(missing_projects, collapse = ";")),
    raw_biologic_features = length(prepared$feature_columns_raw),
    features_retained_for_tree = length(prepared$feature_columns_retained),
    features_dropped = length(dropped_features),
    dropped_features = ifelse(length(dropped_features) == 0, "none", paste(dropped_features, collapse = ";")),
    imputed_values = imputed_total,
    distance_methods_generated = "gower;correlation",
    tree_files_generated = paste(basename(tree_files), collapse = ";"),
    figure_files_generated = paste(basename(c(gower_figure_path, correlation_figure_path, heatmap_figure_path, pca_figure_path)), collapse = ";"),
    bootstrap_iterations_completed = bootstrap_iterations,
    neighbor_joining_succeeded = isTRUE(nj_result$success),
    neighbor_joining_message = nj_result$message,
    warnings_or_limitations = paste(warnings, collapse = " | "),
    feature_qc = prepared$feature_qc
  )
  write_tree_qc_summary(qc_values, qc_path)

  list(
    prepared = prepared,
    scaled = scaled,
    gower_distance = gower_distance,
    correlation_distance = correlation_distance,
    gower_hclust = gower_hclust,
    correlation_hclust = correlation_hclust,
    neighbor_joining = nj_result,
    stability = stability,
    qc = qc_values,
    paths = list(
      scaled = scaled_path,
      gower_distance = gower_distance_path,
      correlation_distance = correlation_distance_path,
      bootstrap = bootstrap_path,
      qc = qc_path,
      gower_tree = gower_tree_path,
      correlation_tree = correlation_tree_path,
      neighbor_joining_tree = if (isTRUE(nj_result$success)) nj_tree_path else NA_character_,
      gower_figure = gower_figure_path,
      correlation_figure = correlation_figure_path,
      heatmap_figure = heatmap_figure_path,
      pca_figure = pca_figure_path
    )
  )
}

main <- function() {
  require_tree_packages()
  args <- parse_tree_args()
  root <- find_project_root()

  feature_matrix_path <- project_path("results", "tables", "level1_feature_matrix.tsv", root = root)
  support_path <- project_path("results", "tables", "level1_feature_support.tsv", root = root)
  feature_qc_path <- project_path("results", "tables", "level1_feature_matrix_qc_summary.tsv", root = root)
  group_map_path <- project_path("config", "cancer_group_map.csv", root = root)
  projects_path <- project_path("data", "interim", "tcga_projects.tsv", root = root)

  feature_matrix <- read_tsv_required(feature_matrix_path, "Level 1 feature matrix")
  read_tsv_required(support_path, "Level 1 feature support table")
  read_tsv_required(feature_qc_path, "Level 1 feature matrix QC summary")
  read_tsv_required(group_map_path, "cancer group map")
  projects <- read_tsv_required(projects_path, "TCGA project metadata")

  required_project_columns <- c("project_id", "project_code")
  missing_project_columns <- setdiff(required_project_columns, names(projects))
  if (length(missing_project_columns) > 0) {
    abort(sprintf("TCGA project metadata is missing required columns: %s", paste(missing_project_columns, collapse = ", ")))
  }
  projects <- unique(projects[, ..required_project_columns])

  log_info("Building Level 1 pan-cancer molecular similarity dendrograms; these are molecular trait trees, not literal species-like phylogenies")
  outputs <- build_level1_tree_outputs(
    feature_matrix = feature_matrix,
    projects = projects,
    max_missing_fraction = args$max_missing_fraction,
    bootstrap_iterations = args$bootstrap_iterations,
    seed = args$seed,
    root = root
  )
  log_info(sprintf(
    "Level 1 tree stage complete: %d projects, %d features retained, %d imputed values, %d bootstrap iterations",
    outputs$qc$included_tcga_projects,
    outputs$qc$features_retained_for_tree,
    outputs$qc$imputed_values,
    outputs$qc$bootstrap_iterations_completed
  ))
  invisible(outputs)
}

if (identical(environment(), globalenv())) {
  main()
}
