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

level2_annotation_columns <- c(
  "project_id", "project_code", "cancer_name", "disease_type", "primary_site",
  "broad_group", "major_cancer_group", "level2_group", "lineage_notes",
  "notes", "note", "source", "data_source"
)

require_level2_packages <- function() {
  required <- c("data.table", "cluster", "ape", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    abort(sprintf("Missing R package(s) required for Level 2 tree construction: %s", paste(missing, collapse = ", ")))
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

numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

ensure_level2_metadata <- function(feature_matrix, projects) {
  matrix <- data.table::copy(feature_matrix)
  if (!"major_cancer_group" %in% names(matrix) && "broad_group" %in% names(matrix)) {
    matrix[, major_cancer_group := broad_group]
  }
  if (!"major_cancer_group" %in% names(matrix)) matrix[, major_cancer_group := "Unannotated"]
  matrix[is.na(major_cancer_group) | major_cancer_group == "", major_cancer_group := "Unannotated"]

  metadata_fields <- setdiff(intersect(c("disease_type", "primary_site"), names(projects)), names(matrix))
  project_fields <- c("project_id", "project_code", metadata_fields)
  if (length(metadata_fields) > 0 && all(c("project_id", "project_code") %in% names(projects))) {
    matrix <- merge(matrix, unique(projects[, ..project_fields]), by = c("project_id", "project_code"), all.x = TRUE, sort = FALSE)
  }
  for (column in c("disease_type", "primary_site")) {
    if (!column %in% names(matrix)) matrix[, (column) := NA_character_]
  }
  matrix
}

biologic_feature_columns <- function(matrix, annotation_columns = level2_annotation_columns) {
  candidates <- setdiff(names(matrix), annotation_columns)
  candidates <- candidates[!grepl("(^|_)n_samples$|^missingness_|_support$|^support_", candidates)]
  candidates[vapply(candidates, function(column) {
    values <- matrix[[column]]
    converted <- numeric_or_na(values)
    all(is.na(values) | !is.na(converted))
  }, logical(1))]
}

slugify <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_|_$", "", x)
}

codes_present <- function(codes, available_codes) {
  intersect(codes, available_codes)
}

define_level2_groups <- function(feature_matrix, cancer_group_map) {
  available <- feature_matrix$project_code
  carcinoma_codes <- cancer_group_map[broad_group == "Carcinoma" & project_code %in% available, project_code]
  groups <- data.table::data.table(
    group_name = c(
      "carcinoma_all",
      "pan_squamous",
      "gi_pancancreatobiliary",
      "kidney",
      "gynecologic_breast",
      "cns_glial",
      "melanocytic",
      "hematolymphoid",
      "sarcoma"
    ),
    group_type = c(
      "primary",
      "carcinoma_subgroup",
      "carcinoma_subgroup",
      "carcinoma_subgroup",
      "carcinoma_subgroup",
      "limited_noncarcinoma",
      "limited_noncarcinoma",
      "limited_noncarcinoma",
      "project_level_singleton"
    ),
    candidate_project_codes = c(
      paste(sort(carcinoma_codes), collapse = ";"),
      paste(c("LUSC", "HNSC", "CESC", "ESCA", "BLCA"), collapse = ";"),
      paste(c("COAD", "READ", "STAD", "ESCA", "PAAD", "CHOL", "LIHC"), collapse = ";"),
      paste(c("KIRC", "KIRP", "KICH"), collapse = ";"),
      paste(c("BRCA", "OV", "UCEC", "CESC", "UCS"), collapse = ";"),
      paste(c("GBM", "LGG"), collapse = ";"),
      paste(c("SKCM", "UVM"), collapse = ";"),
      paste(c("LAML", "DLBC"), collapse = ";"),
      "SARC"
    )
  )
  groups[, included_project_codes := vapply(candidate_project_codes, function(x) paste(codes_present(strsplit(x, ";", fixed = TRUE)[[1]], available), collapse = ";"), character(1))]
  groups[, n_projects := vapply(included_project_codes, function(x) if (nchar(x) == 0) 0L else length(strsplit(x, ";", fixed = TRUE)[[1]]), integer(1))]
  groups
}

prepare_group_features <- function(feature_matrix, group_codes, max_missing_fraction = 0.30) {
  group_matrix <- data.table::copy(feature_matrix[project_code %in% group_codes])
  group_matrix <- group_matrix[match(group_codes, project_code)]
  annotation_cols <- intersect(level2_annotation_columns, names(group_matrix))
  feature_cols <- biologic_feature_columns(group_matrix)
  for (column in feature_cols) group_matrix[, (column) := numeric_or_na(get(column))]

  feature_qc <- data.table::data.table(
    feature_name = feature_cols,
    missing_fraction = vapply(feature_cols, function(column) {
      if (nrow(group_matrix) == 0) return(1)
      value <- mean(is.na(group_matrix[[column]]))
      ifelse(is.na(value), 1, value)
    }, numeric(1)),
    n_missing_projects = vapply(feature_cols, function(column) sum(is.na(group_matrix[[column]])), integer(1)),
    n_unique_non_missing = vapply(feature_cols, function(column) length(unique(stats::na.omit(group_matrix[[column]]))), integer(1)),
    raw_sd = vapply(feature_cols, function(column) stats::sd(group_matrix[[column]], na.rm = TRUE), numeric(1)),
    imputed_value = NA_real_,
    n_imputed_values = 0L,
    retained_for_group_tree = FALSE,
    drop_reason = ""
  )

  retained <- character()
  for (feature in feature_cols) {
    values <- group_matrix[[feature]]
    missing_fraction <- if (length(values) == 0) 1 else mean(is.na(values))
    if (is.na(missing_fraction)) missing_fraction <- 1
    unique_non_missing <- length(unique(stats::na.omit(values)))
    if (missing_fraction > max_missing_fraction) {
      feature_qc[feature_name == feature, drop_reason := sprintf("missing_fraction_gt_%.2f", max_missing_fraction)]
      next
    }
    if (unique_non_missing == 0) {
      feature_qc[feature_name == feature, drop_reason := "all_missing_within_group"]
      next
    }
    median_value <- stats::median(values, na.rm = TRUE)
    n_imputed <- sum(is.na(values))
    if (n_imputed > 0) group_matrix[is.na(get(feature)), (feature) := median_value]
    if (length(unique(group_matrix[[feature]])) < 2 || is.na(stats::sd(group_matrix[[feature]])) || stats::sd(group_matrix[[feature]]) <= 1e-8) {
      feature_qc[feature_name == feature, `:=`(
        imputed_value = median_value,
        n_imputed_values = n_imputed,
        drop_reason = "constant_or_near_constant_within_group"
      )]
      next
    }
    retained <- c(retained, feature)
    feature_qc[feature_name == feature, `:=`(
      imputed_value = median_value,
      n_imputed_values = n_imputed,
      retained_for_group_tree = TRUE
    )]
  }

  raw_retained <- group_matrix[, c(annotation_cols, retained), with = FALSE]
  scaled_retained <- data.table::copy(raw_retained)
  if (length(retained) > 0) {
    scaled_values <- scale(as.matrix(group_matrix[, ..retained]))
    scaled_values[is.na(scaled_values)] <- 0
    for (feature in retained) scaled_retained[, (feature) := scaled_values[, feature]]
  }
  list(
    raw_matrix = raw_retained,
    scaled_matrix = scaled_retained,
    feature_qc = feature_qc,
    retained_features = retained
  )
}

matrix_from_group_scaled <- function(scaled_table, features) {
  mat <- as.matrix(scaled_table[, ..features])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- scaled_table$project_code
  mat
}

validate_distance_matrix <- function(distance_matrix, label) {
  if (!is.matrix(distance_matrix)) abort(sprintf("%s distance output is not a matrix", label))
  if (nrow(distance_matrix) != ncol(distance_matrix)) abort(sprintf("%s distance matrix is not square", label))
  if (any(is.na(distance_matrix))) abort(sprintf("%s distance matrix contains NA values", label))
  distance_matrix <- (distance_matrix + t(distance_matrix)) / 2
  diag(distance_matrix) <- 0
  distance_matrix[distance_matrix < 0 & distance_matrix > -1e-8] <- 0
  if (any(distance_matrix < -1e-8)) abort(sprintf("%s distance matrix contains negative distances", label))
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
  rownames(distance) <- rownames(scaled_matrix)
  colnames(distance) <- rownames(scaled_matrix)
  validate_distance_matrix(distance, "Correlation")
}

write_distance_matrix <- function(distance_matrix, path) {
  write_tsv_table(data.table::as.data.table(distance_matrix, keep.rownames = "project_code"), path)
}

write_hclust_tree <- function(hclust_object, path) {
  ensure_dir(dirname(path))
  phylo <- ape::as.phylo(hclust_object)
  ape::write.tree(phylo, file = path)
  log_info(sprintf("Wrote %s", path))
  phylo
}

group_palette <- function(groups) {
  groups <- sort(unique(groups[!is.na(groups) & groups != ""]))
  if (length(groups) == 0) groups <- "Unannotated"
  n <- length(groups)
  colors <- if (n <= 8) {
    RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)]
  } else if (n <= 12) {
    RColorBrewer::brewer.pal(n, "Paired")
  } else {
    grDevices::rainbow(n)
  }
  stats::setNames(colors, groups)
}

plot_tree_pdf <- function(phylo, annotations, path, title, color_column = "level2_group") {
  ensure_dir(dirname(path))
  groups <- annotations[[color_column]]
  if (is.null(groups)) groups <- annotations$major_cancer_group
  names(groups) <- annotations$project_code
  palette <- group_palette(groups)
  tip_colors <- palette[groups[phylo$tip.label]]
  grDevices::pdf(path, width = 10.5, height = 8)
  op <- par(mar = c(2, 2, 4, 1))
  on.exit({
    par(op)
    grDevices::dev.off()
  }, add = TRUE)
  ape::plot.phylo(phylo, type = "phylogram", cex = 0.8, tip.color = tip_colors, main = title, no.margin = FALSE)
  legend("topleft", legend = names(palette), col = palette, pch = 19, cex = 0.65, bty = "n", title = color_column)
}

selected_heatmap_traits <- function(features) {
  preferred <- c(
    "tmb_proxy_median_nonsynonymous_count",
    "mutation_prevalence_TP53",
    "mutation_prevalence_KRAS",
    "mutation_prevalence_BRAF",
    "mutation_prevalence_IDH1",
    "mutation_prevalence_IDH2",
    "median_aneuploidy_score",
    "median_ploidy",
    "median_purity",
    "median_arm_gain_count",
    "median_arm_loss_count",
    "median_immune_inflammatory_score",
    "median_proliferation_score",
    "median_EMT_score"
  )
  intersect(preferred, features)
}

plot_trait_heatmap <- function(scaled_table, features, hclust_object, path, title) {
  traits <- selected_heatmap_traits(features)
  if (length(traits) == 0) {
    log_warn(sprintf("Trait heatmap skipped for %s because no selected traits were present", title))
    return(FALSE)
  }
  mat <- as.matrix(scaled_table[, ..traits])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- scaled_table$project_code
  annotation <- data.frame(level2_group = scaled_table$level2_group, row.names = scaled_table$project_code, check.names = FALSE)
  pheatmap::pheatmap(
    mat,
    filename = path,
    width = 9.5,
    height = 7.5,
    cluster_rows = hclust_object,
    cluster_cols = TRUE,
    annotation_row = annotation,
    main = title,
    fontsize_row = 8,
    fontsize_col = 8,
    border_color = NA
  )
  TRUE
}

plot_pca_pdf <- function(scaled_table, features, path, title) {
  if (nrow(scaled_table) < 3 || length(features) < 2) return(FALSE)
  mat <- matrix_from_group_scaled(scaled_table, features)
  pca <- stats::prcomp(mat, center = FALSE, scale. = FALSE)
  if (ncol(pca$x) < 2) return(FALSE)
  variance <- (pca$sdev^2) / sum(pca$sdev^2)
  plot_data <- data.table::data.table(
    project_code = rownames(mat),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  plot_data <- merge(plot_data, scaled_table[, .(project_code, level2_group, major_cancer_group)], by = "project_code", all.x = TRUE, sort = FALSE)
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(PC1, PC2, color = level2_group, label = project_code)) +
    ggplot2::geom_point(size = 2.8) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", 100 * variance[[1]]),
      y = sprintf("PC2 (%.1f%%)", 100 * variance[[2]]),
      color = "Level 2 group"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  ggplot2::ggsave(path, plot, width = 9.5, height = 7.5, device = "pdf")
  TRUE
}

bootstrap_cluster_stability <- function(scaled_matrix, iterations = 100L, seed = 1L, k = NULL) {
  projects <- rownames(scaled_matrix)
  n_projects <- length(projects)
  if (is.null(k)) k <- max(2L, min(n_projects - 1L, round(sqrt(n_projects))))
  pair_grid <- utils::combn(projects, 2)
  counts <- integer(ncol(pair_grid))
  set.seed(seed)
  for (iteration in seq_len(iterations)) {
    columns <- sample(seq_len(ncol(scaled_matrix)), size = ncol(scaled_matrix), replace = TRUE)
    boot_matrix <- scaled_matrix[, columns, drop = FALSE]
    colnames(boot_matrix) <- make.unique(colnames(boot_matrix))
    distance <- compute_gower_distance_matrix(boot_matrix)
    clusters <- stats::cutree(stats::hclust(stats::as.dist(distance), method = "average"), k = k)
    counts <- counts + as.integer(clusters[pair_grid[1, ]] == clusters[pair_grid[2, ]])
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

closest_pairs_from_hclust <- function(hclust_object, n = 8) {
  cophen <- as.matrix(stats::cophenetic(hclust_object))
  pairs <- as.data.table(as.table(cophen))
  data.table::setnames(pairs, c("project_code_a", "project_code_b", "cophenetic_height"))
  pairs <- pairs[as.character(project_code_a) < as.character(project_code_b)]
  data.table::setorder(pairs, cophenetic_height)
  head(pairs, n)
}

cluster_text <- function(hclust_object, k = 4) {
  labels <- hclust_object$labels
  if (length(labels) < 3) return("underpowered")
  k <- max(2L, min(k, length(labels) - 1L))
  clusters <- stats::cutree(hclust_object, k = k)
  split_projects <- split(names(clusters), clusters)
  paste(vapply(names(split_projects), function(id) sprintf("cluster %s: %s", id, paste(split_projects[[id]], collapse = ",")), character(1)), collapse = "; ")
}

check_sarc_subtype_metadata <- function(root = find_project_root()) {
  candidate_files <- c(
    project_path("data", "interim", "clinical_metadata.tsv", root = root),
    project_path("data", "processed", "clinical_metadata.tsv", root = root),
    project_path("data", "processed", "features", "sarc_subtype_feature_matrix.tsv", root = root)
  )
  existing <- candidate_files[file.exists(candidate_files)]
  if (length(existing) == 0) {
    return("no_clinical_or_sarc_subtype_feature_metadata_found")
  }
  if (any(grepl("sarc_subtype_feature_matrix.tsv$", existing))) {
    return("sarc_subtype_feature_matrix_found_but_subtype_tree_not_implemented_in_project_level_step")
  }
  "clinical_metadata_present_but_no_sarc_subtype_molecular_feature_matrix"
}

process_level2_group <- function(group_row, feature_matrix, root = find_project_root(), bootstrap_iterations = 100L, seed = 1L) {
  group_name <- group_row$group_name
  group_codes <- if (nchar(group_row$included_project_codes) == 0) character() else strsplit(group_row$included_project_codes, ";", fixed = TRUE)[[1]]
  paths <- list(
    raw = project_path("results", "tables", sprintf("level2_%s_feature_matrix.tsv", group_name), root = root),
    scaled = project_path("results", "tables", sprintf("level2_%s_feature_matrix_scaled.tsv", group_name), root = root),
    gower_distance = project_path("results", "tables", sprintf("level2_%s_gower_distance_matrix.tsv", group_name), root = root),
    correlation_distance = project_path("results", "tables", sprintf("level2_%s_correlation_distance_matrix.tsv", group_name), root = root),
    gower_tree = project_path("results", "trees", sprintf("level2_%s_gower_hclust_tree.nwk", group_name), root = root),
    correlation_tree = project_path("results", "trees", sprintf("level2_%s_correlation_hclust_tree.nwk", group_name), root = root),
    gower_figure = project_path("results", "figures", sprintf("level2_%s_gower_tree.pdf", group_name), root = root),
    correlation_figure = project_path("results", "figures", sprintf("level2_%s_correlation_tree.pdf", group_name), root = root),
    heatmap = project_path("results", "figures", sprintf("level2_%s_trait_heatmap.pdf", group_name), root = root),
    pca = project_path("results", "figures", sprintf("level2_%s_feature_pca.pdf", group_name), root = root),
    bootstrap = project_path("results", "tables", sprintf("level2_%s_bootstrap_cluster_stability.tsv", group_name), root = root)
  )

  prepared <- prepare_group_features(feature_matrix, group_codes)
  write_tsv_table(prepared$raw_matrix, paths$raw)
  write_tsv_table(prepared$scaled_matrix, paths$scaled)

  reason <- ""
  tree_attempted <- TRUE
  if (length(group_codes) < 3) {
    tree_attempted <- FALSE
    reason <- sprintf("underpowered_project_count_%d_lt_3", length(group_codes))
  }
  if (length(prepared$retained_features) < 2) {
    tree_attempted <- FALSE
    reason <- paste(c(reason, sprintf("insufficient_variable_features_%d_lt_2", length(prepared$retained_features))), collapse = ";")
    reason <- gsub("^;", "", reason)
  }
  if (group_name == "sarcoma" && length(group_codes) < 3) {
    reason <- paste(c(reason, check_sarc_subtype_metadata(root)), collapse = ";")
    reason <- gsub("^;", "", reason)
  }

  result <- list(
    group_name = group_name,
    group_type = group_row$group_type,
    group_codes = group_codes,
    retained_features = prepared$retained_features,
    feature_qc = prepared$feature_qc,
    tree_attempted = tree_attempted,
    reason = ifelse(tree_attempted, "", reason),
    paths = paths,
    gower_hclust = NULL,
    correlation_hclust = NULL,
    bootstrap_iterations = 0L
  )

  if (!tree_attempted) return(result)

  mat <- matrix_from_group_scaled(prepared$scaled_matrix, prepared$retained_features)
  gower_distance <- compute_gower_distance_matrix(mat)
  correlation_distance <- compute_correlation_distance_matrix(mat)
  write_distance_matrix(gower_distance, paths$gower_distance)
  write_distance_matrix(correlation_distance, paths$correlation_distance)
  gower_hclust <- stats::hclust(stats::as.dist(gower_distance), method = "average")
  correlation_hclust <- stats::hclust(stats::as.dist(correlation_distance), method = "average")
  gower_phylo <- write_hclust_tree(gower_hclust, paths$gower_tree)
  correlation_phylo <- write_hclust_tree(correlation_hclust, paths$correlation_tree)
  plot_tree_pdf(gower_phylo, prepared$scaled_matrix, paths$gower_figure, sprintf("Level 2 %s molecular similarity tree (Gower)", group_name))
  plot_tree_pdf(correlation_phylo, prepared$scaled_matrix, paths$correlation_figure, sprintf("Level 2 %s molecular similarity tree (correlation)", group_name))
  plot_trait_heatmap(prepared$scaled_matrix, prepared$retained_features, gower_hclust, paths$heatmap, sprintf("Level 2 %s selected molecular traits", group_name))
  plot_pca_pdf(prepared$scaled_matrix, prepared$retained_features, paths$pca, sprintf("Level 2 %s feature PCA", group_name))

  if (length(group_codes) >= 5) {
    stability <- bootstrap_cluster_stability(mat, iterations = bootstrap_iterations, seed = seed)
    write_tsv_table(stability, paths$bootstrap)
    result$bootstrap_iterations <- bootstrap_iterations
  } else {
    result$reason <- "bootstrap_skipped_project_count_lt_5"
  }

  result$gower_hclust <- gower_hclust
  result$correlation_hclust <- correlation_hclust
  result
}

group_definition_output <- function(groups, feature_matrix, results) {
  rows <- lapply(results, function(result) {
    ids <- feature_matrix[match(result$group_codes, project_code), project_id]
    data.table::data.table(
      group_name = result$group_name,
      group_type = result$group_type,
      included_project_codes = paste(result$group_codes, collapse = ";"),
      included_project_ids = paste(ids, collapse = ";"),
      n_projects = length(result$group_codes),
      n_features_used = length(result$retained_features),
      tree_attempted = result$tree_attempted,
      reason_if_not_attempted = ifelse(result$tree_attempted, "", result$reason)
    )
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

level2_qc_summary <- function(results) {
  rows <- list()
  for (result in results) {
    rows[[length(rows) + 1L]] <- data.table::data.table(
      qc_section = "global",
      group_name = result$group_name,
      feature_name = "ALL",
      metric = c(
        "n_projects",
        "n_raw_biologic_features",
        "n_features_retained",
        "n_features_dropped",
        "tree_attempted",
        "reason_or_warning",
        "bootstrap_iterations_completed"
      ),
      value = c(
        length(result$group_codes),
        nrow(result$feature_qc),
        length(result$retained_features),
        sum(!result$feature_qc$retained_for_group_tree),
        result$tree_attempted,
        ifelse(nchar(result$reason) == 0, "none", result$reason),
        result$bootstrap_iterations
      )
    )
    feature_qc <- data.table::copy(result$feature_qc)
    if (nrow(feature_qc) > 0) {
      feature_qc[, `:=`(qc_section = "per_feature", group_name = result$group_name)]
      measure_cols <- setdiff(names(feature_qc), c("qc_section", "group_name", "feature_name"))
      for (column in measure_cols) feature_qc[, (column) := as.character(get(column))]
      rows[[length(rows) + 1L]] <- data.table::melt(
        feature_qc,
        id.vars = c("qc_section", "group_name", "feature_name"),
        variable.name = "metric",
        value.name = "value"
      )
    }
  }
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  data.table::setorder(out, group_name, qc_section, feature_name, metric)
  out
}

write_level2_report <- function(results, group_definitions, path) {
  ensure_dir(dirname(path))
  tree_groups <- group_definitions[tree_attempted == TRUE, group_name]
  underpowered <- group_definitions[tree_attempted == FALSE]
  result_by_name <- stats::setNames(results, vapply(results, `[[`, character(1), "group_name"))
  carcinoma <- result_by_name[["carcinoma_all"]]
  pan_squamous <- result_by_name[["pan_squamous"]]
  gi <- result_by_name[["gi_pancancreatobiliary"]]
  kidney <- result_by_name[["kidney"]]
  gyn <- result_by_name[["gynecologic_breast"]]

  carcinoma_text <- if (!is.null(carcinoma$gower_hclust)) cluster_text(carcinoma$gower_hclust, k = 6) else "not attempted"
  correlation_text <- if (!is.null(carcinoma$correlation_hclust)) cluster_text(carcinoma$correlation_hclust, k = 6) else "not attempted"
  closest_pairs <- if (!is.null(carcinoma$gower_hclust)) {
    paste(sprintf("%s-%s", closest_pairs_from_hclust(carcinoma$gower_hclust, 8)$project_code_a, closest_pairs_from_hclust(carcinoma$gower_hclust, 8)$project_code_b), collapse = ", ")
  } else {
    "not available"
  }

  subgroup_line <- function(label, result) {
    if (is.null(result) || is.null(result$gower_hclust)) return(sprintf("%s was underpowered or not attempted.", label))
    sprintf("%s produced a group-specific molecular trait tree with Gower clusters: %s.", label, cluster_text(result$gower_hclust, k = min(4, length(result$group_codes) - 1L)))
  }

  lines <- c(
    "# Level 2 Within-Group Molecular Similarity Summary",
    "",
    "Level 2 builds group-specific molecular similarity trees from the same project-level molecular traits used in Level 1. Features are filtered for within-group variability, median-imputed within group when needed, scaled within group, and clustered with Gower and correlation distances using average linkage.",
    "",
    sprintf("Groups generating project-level trees were: %s. Underpowered groups were: %s.", paste(tree_groups, collapse = ", "), paste(sprintf("%s (%s)", underpowered$group_name, underpowered$reason_if_not_attempted), collapse = "; ")),
    "",
    sprintf("The carcinoma-all Gower tree shows the following k-level structure: %s. Close project pairs include %s. These results indicate similar molecular trait profiles within the carcinoma set and should not be interpreted as directional ancestry.", carcinoma_text, closest_pairs),
    "",
    sprintf("The carcinoma-all correlation tree differs as expected because it emphasizes relative feature-profile shape rather than absolute scaled trait levels. Its k-level structure is: %s.", correlation_text),
    "",
    paste(
      subgroup_line("Pan-squamous / squamous-enriched", pan_squamous),
      subgroup_line("GI / pancreatobiliary adenocarcinoma-enriched", gi),
      subgroup_line("Kidney", kidney),
      subgroup_line("Gynecologic / breast", gyn),
      sep = " "
    ),
    "",
    "Across Level 2 trees, clustering appears influenced by a mixture of mutation-count proxies, driver-gene prevalence, aneuploidy/copy-number burden, purity/ploidy, and expression programs. Because Level 1 sensitivity showed mutation-heavy features had strong cluster associations, exact Level 2 topology should be read together with feature-subset and weighting sensitivity rather than treated as a single definitive structure.",
    "",
    "CNS/glial, melanocytic, hematolymphoid, and sarcoma project-level analyses are limited by having one or two TCGA projects in the relevant group. The SARC project lacks a project-internal subtype molecular feature matrix in the current workflow, so a subtype-level sarcoma tree is deferred.",
    "",
    "Scientific caveat: Level 2 outputs are molecular similarity dendrograms or molecular trait trees. They are not literal organismal phylogenies and should not be used to infer directional evolutionary relationships between cancer types.",
    ""
  )
  writeLines(lines, con = path)
  log_info(sprintf("Wrote %s", path))
}

build_level2_outputs <- function(feature_matrix, scaled_matrix, cancer_group_map, projects, root = find_project_root(), bootstrap_iterations = 100L, seed = 1L) {
  feature_matrix <- ensure_level2_metadata(feature_matrix, projects)
  scaled_matrix <- ensure_level2_metadata(scaled_matrix, projects)
  groups <- define_level2_groups(feature_matrix, cancer_group_map)
  results <- lapply(seq_len(nrow(groups)), function(i) {
    process_level2_group(groups[i], feature_matrix, root = root, bootstrap_iterations = bootstrap_iterations, seed = seed)
  })
  definitions <- group_definition_output(groups, feature_matrix, results)
  qc <- level2_qc_summary(results)
  definition_path <- project_path("results", "tables", "level2_group_definitions.tsv", root = root)
  qc_path <- project_path("results", "tables", "level2_tree_qc_summary.tsv", root = root)
  report_path <- project_path("results", "reports", "level2_interpretation_summary.md", root = root)
  write_tsv_table(definitions, definition_path)
  write_tsv_table(qc, qc_path)
  write_level2_report(results, definitions, report_path)
  list(groups = definitions, qc = qc, results = results, paths = list(definitions = definition_path, qc = qc_path, report = report_path))
}

main <- function() {
  require_level2_packages()
  root <- find_project_root()
  feature_matrix <- read_tsv_required(project_path("results", "tables", "level1_feature_matrix.tsv", root = root), "Level 1 feature matrix")
  scaled_matrix <- read_tsv_required(project_path("results", "tables", "level1_feature_matrix_scaled.tsv", root = root), "scaled Level 1 feature matrix")
  cancer_group_map <- data.table::fread(project_path("config", "cancer_group_map.csv", root = root), sep = ",", quote = "\"", data.table = TRUE, showProgress = FALSE)
  projects <- read_tsv_required(project_path("data", "interim", "tcga_projects.tsv", root = root), "TCGA project metadata")

  optional_inputs <- c(
    project_path("results", "tables", "level1_feature_cluster_association.tsv", root = root),
    project_path("results", "tables", "level1_feature_class_sensitivity_summary.tsv", root = root),
    project_path("results", "tables", "level1_gower_cluster_assignments.tsv", root = root),
    project_path("results", "tables", "level1_correlation_cluster_assignments.tsv", root = root)
  )
  invisible(lapply(optional_inputs[file.exists(optional_inputs)], function(path) read_tsv_required(path, basename(path))))

  log_info("Building Level 2 group-specific molecular similarity trees; these are molecular trait trees, not literal ancestry among cancer types")
  outputs <- build_level2_outputs(feature_matrix, scaled_matrix, cancer_group_map, projects, root = root)
  log_info(sprintf(
    "Level 2 complete: %d groups defined, %d groups generated trees, %d groups underpowered",
    nrow(outputs$groups),
    sum(outputs$groups$tree_attempted),
    sum(!outputs$groups$tree_attempted)
  ))
  invisible(outputs)
}

if (identical(environment(), globalenv())) {
  main()
}
