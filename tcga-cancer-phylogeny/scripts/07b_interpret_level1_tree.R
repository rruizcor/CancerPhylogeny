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

interpretation_annotation_columns <- c(
  "project_id", "project_code", "cancer_name", "disease_type", "primary_site",
  "broad_group", "major_cancer_group", "level2_group", "lineage_notes",
  "notes", "note", "source", "data_source"
)

require_interpretation_packages <- function() {
  required <- c("data.table", "cluster", "ggplot2", "pheatmap", "RColorBrewer")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    abort(sprintf("Missing R package(s) required for Level 1 interpretation: %s", paste(missing, collapse = ", ")))
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

ensure_metadata <- function(matrix, projects) {
  matrix <- data.table::copy(matrix)
  if (!"major_cancer_group" %in% names(matrix) && "broad_group" %in% names(matrix)) {
    matrix[, major_cancer_group := broad_group]
  }
  if (!"major_cancer_group" %in% names(matrix)) {
    matrix[, major_cancer_group := "Unannotated"]
  }
  matrix[is.na(major_cancer_group) | major_cancer_group == "", major_cancer_group := "Unannotated"]

  project_fields <- intersect(c("project_id", "project_code", "disease_type", "primary_site"), names(projects))
  matrix <- merge(matrix, unique(projects[, ..project_fields]), by = c("project_id", "project_code"), all.x = TRUE, sort = FALSE)
  for (column in c("disease_type", "primary_site")) {
    if (!column %in% names(matrix)) matrix[, (column) := NA_character_]
  }
  matrix
}

feature_columns <- function(matrix, annotation_columns = interpretation_annotation_columns) {
  candidates <- setdiff(names(matrix), annotation_columns)
  candidates[vapply(candidates, function(column) {
    values <- matrix[[column]]
    converted <- numeric_or_na(values)
    all(is.na(values) | !is.na(converted))
  }, logical(1))]
}

matrix_from_scaled <- function(scaled_matrix) {
  features <- feature_columns(scaled_matrix)
  mat <- as.matrix(scaled_matrix[, ..features])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- scaled_matrix$project_code
  mat
}

read_distance_matrix <- function(path) {
  table <- read_tsv_required(path, basename(path))
  if (!"project_code" %in% names(table)) abort(sprintf("Distance matrix lacks project_code column: %s", path))
  matrix <- as.matrix(table[, setdiff(names(table), "project_code"), with = FALSE])
  storage.mode(matrix) <- "numeric"
  rownames(matrix) <- table$project_code
  matrix <- (matrix + t(matrix)) / 2
  diag(matrix) <- 0
  matrix
}

cluster_assignments_from_distance <- function(distance_matrix, metadata, k_values = 3:8) {
  hc <- stats::hclust(stats::as.dist(distance_matrix), method = "average")
  assignments <- metadata[, .(project_id, project_code, major_cancer_group, disease_type, primary_site)]
  assignments <- assignments[match(rownames(distance_matrix), project_code)]
  for (k in k_values) {
    clusters <- stats::cutree(hc, k = k)
    assignments[, (sprintf("cluster_k%d", k)) := as.integer(clusters[project_code])]
  }
  list(assignments = assignments, hclust = hc)
}

median_feature <- function(data, feature) {
  if (!feature %in% names(data)) return(NA_real_)
  stats::median(numeric_or_na(data[[feature]]), na.rm = TRUE)
}

top_group_text <- function(values) {
  values <- unique(values[!is.na(values) & values != ""])
  if (length(values) == 0) "none" else paste(sort(values), collapse = ";")
}

interpret_cluster_label <- function(cluster_data) {
  projects <- cluster_data$project_code
  groups <- cluster_data$major_cancer_group
  tmb <- median_feature(cluster_data, "tmb_proxy_median_nonsynonymous_count")
  tp53 <- median_feature(cluster_data, "mutation_prevalence_TP53")
  kras <- median_feature(cluster_data, "mutation_prevalence_KRAS")
  aneuploidy <- median_feature(cluster_data, "median_aneuploidy_score")
  immune <- median_feature(cluster_data, "median_immune_inflammatory_score")

  if (length(projects) == 1 && projects %in% c("DLBC", "LAML")) return("hematolymphoid outlier")
  if (length(projects) == 1 && projects == "SKCM") return("high mutation burden / immune-active outlier")
  if (length(projects) == 1 && projects %in% c("UCEC", "UCS")) return("gynecologic outlier with distinctive molecular profile")
  if (all(c("COAD", "READ") %in% projects) || sum(projects %in% c("COAD", "READ", "STAD", "ESCA", "PAAD", "CHOL", "LIHC")) >= 3) {
    return("GI carcinoma enriched")
  }
  if (sum(projects %in% c("CESC", "HNSC", "LUSC", "ESCA", "BLCA")) >= 3) return("squamous-like carcinoma enriched")
  if (sum(projects %in% c("KICH", "KIRC", "KIRP", "GBM", "LGG", "ACC", "PCPG", "THCA", "PRAD", "UVM")) >= 5) {
    return("kidney/CNS/endocrine low-mutation cluster")
  }
  if (!is.na(tmb) && tmb >= 150 && !is.na(immune) && immune > 0) return("high mutation burden / immune-active")
  if (!is.na(aneuploidy) && aneuploidy >= 12 && !is.na(tp53) && tp53 >= 0.35) return("high aneuploidy carcinoma cluster")
  if (sum(groups == "Carcinoma", na.rm = TRUE) >= max(2, ceiling(0.6 * length(projects)))) return("carcinoma enriched")
  if (!is.na(kras) && kras >= 0.20) return("KRAS-altered carcinoma enriched")
  "mixed lineage cluster"
}

cluster_summaries <- function(assignments, raw_matrix, method_name, k_values = 3:8) {
  merged <- merge(assignments, raw_matrix, by = c("project_id", "project_code"), all.x = TRUE, sort = FALSE)
  summaries <- list()
  traits <- c(
    median_mutation_count_proxy = "tmb_proxy_median_nonsynonymous_count",
    median_TP53_prevalence = "mutation_prevalence_TP53",
    median_KRAS_prevalence = "mutation_prevalence_KRAS",
    median_BRAF_prevalence = "mutation_prevalence_BRAF",
    median_IDH1_prevalence = "mutation_prevalence_IDH1",
    median_IDH2_prevalence = "mutation_prevalence_IDH2",
    median_purity = "median_purity",
    median_ploidy = "median_ploidy",
    median_aneuploidy_score = "median_aneuploidy_score",
    median_arm_gain_count = "median_arm_gain_count",
    median_arm_loss_count = "median_arm_loss_count",
    median_immune_inflammatory_score = "median_immune_inflammatory_score",
    median_proliferation_score = "median_proliferation_score",
    median_EMT_score = "median_EMT_score"
  )
  for (k in k_values) {
    cluster_col <- sprintf("cluster_k%d", k)
    for (cluster_id in sort(unique(merged[[cluster_col]]))) {
      cluster_data <- merged[get(cluster_col) == cluster_id]
      row <- data.table::data.table(
        distance_method = method_name,
        k = k,
        cluster_id = cluster_id,
        n_projects = nrow(cluster_data),
        projects_included = paste(cluster_data$project_code, collapse = ";"),
        major_cancer_groups_represented = top_group_text(cluster_data$major_cancer_group)
      )
      for (trait_name in names(traits)) {
        row[, (trait_name) := median_feature(cluster_data, traits[[trait_name]])]
      }
      row[, interpretation_label := interpret_cluster_label(cluster_data)]
      summaries[[length(summaries) + 1L]] <- row
    }
  }
  data.table::rbindlist(summaries, use.names = TRUE, fill = TRUE)
}

feature_class <- function(feature) {
  if (grepl("^tmb_proxy_|^mutation_prevalence_", feature)) return("mutation")
  if (grepl("purity|ploidy|aneuploidy|arm_|copy_number|fraction_genome_altered", feature)) return("copy_number_aneuploidy")
  if (grepl("^median_.*score$|^median_lineage_or_tissue_PC", feature)) return("expression")
  "other"
}

feature_pca_loadings <- function(scaled_matrix) {
  pca <- stats::prcomp(scaled_matrix, center = FALSE, scale. = FALSE)
  pcs <- seq_len(min(6L, ncol(pca$rotation)))
  rows <- list()
  for (pc in pcs) {
    loadings <- pca$rotation[, pc]
    rows[[length(rows) + 1L]] <- data.table::data.table(
      feature_name = names(loadings),
      feature_class = vapply(names(loadings), feature_class, character(1)),
      pc = paste0("PC", pc),
      loading = as.numeric(loadings),
      abs_loading = abs(as.numeric(loadings))
    )
  }
  out <- data.table::rbindlist(rows)
  data.table::setorder(out, pc, -abs_loading)
  out
}

feature_variance_summary <- function(raw_matrix, scaled_matrix) {
  features <- colnames(scaled_matrix)
  out <- data.table::data.table(
    feature_name = features,
    feature_class = vapply(features, feature_class, character(1)),
    raw_variance = vapply(features, function(feature) stats::var(numeric_or_na(raw_matrix[[feature]]), na.rm = TRUE), numeric(1)),
    scaled_variance = vapply(features, function(feature) stats::var(scaled_matrix[, feature], na.rm = TRUE), numeric(1)),
    raw_median = vapply(features, function(feature) stats::median(numeric_or_na(raw_matrix[[feature]]), na.rm = TRUE), numeric(1)),
    raw_iqr = vapply(features, function(feature) stats::IQR(numeric_or_na(raw_matrix[[feature]]), na.rm = TRUE), numeric(1))
  )
  data.table::setorder(out, feature_class, -scaled_variance, feature_name)
  out
}

feature_cluster_association <- function(raw_matrix, scaled_matrix, assignments, cluster_col = "cluster_k6") {
  features <- colnames(scaled_matrix)
  clusters <- assignments[[cluster_col]][match(rownames(scaled_matrix), assignments$project_code)]
  rows <- lapply(features, function(feature) {
    values <- scaled_matrix[, feature]
    test <- tryCatch(stats::kruskal.test(values, as.factor(clusters)), error = function(err) NULL)
    h <- if (is.null(test)) NA_real_ else as.numeric(test$statistic)
    p <- if (is.null(test)) NA_real_ else test$p.value
    n <- length(values)
    k <- length(unique(clusters))
    epsilon_squared <- if (!is.na(h) && n > k) max(0, (h - k + 1) / (n - k)) else NA_real_
    data.table::data.table(
      feature_name = feature,
      feature_class = feature_class(feature),
      kruskal_wallis_statistic = h,
      p_value = p,
      epsilon_squared = epsilon_squared,
      raw_variance = stats::var(numeric_or_na(raw_matrix[[feature]]), na.rm = TRUE)
    )
  })
  out <- data.table::rbindlist(rows)
  out[, fdr := stats::p.adjust(p_value, method = "BH")]
  data.table::setorder(out, fdr, -epsilon_squared, feature_name)
  out
}

choose2 <- function(x) {
  x * (x - 1) / 2
}

adjusted_rand_index <- function(labels_a, labels_b) {
  labels_a <- as.factor(labels_a)
  labels_b <- as.factor(labels_b)
  tab <- table(labels_a, labels_b)
  sum_comb <- sum(choose2(tab))
  row_comb <- sum(choose2(rowSums(tab)))
  col_comb <- sum(choose2(colSums(tab)))
  total_comb <- choose2(sum(tab))
  expected <- row_comb * col_comb / total_comb
  max_index <- 0.5 * (row_comb + col_comb)
  denom <- max_index - expected
  if (denom == 0) return(1)
  (sum_comb - expected) / denom
}

compute_gower_distance <- function(matrix, weights = NULL) {
  if (is.null(weights)) {
    distance <- cluster::daisy(as.data.frame(matrix), metric = "gower")
  } else {
    distance <- cluster::daisy(as.data.frame(matrix), metric = "gower", weights = weights)
  }
  out <- as.matrix(distance)
  rownames(out) <- rownames(matrix)
  colnames(out) <- rownames(matrix)
  diag(out) <- 0
  (out + t(out)) / 2
}

cluster_k6_from_matrix <- function(matrix, weights = NULL) {
  distance <- compute_gower_distance(matrix, weights = weights)
  hc <- stats::hclust(stats::as.dist(distance), method = "average")
  stats::cutree(hc, k = 6)
}

changed_projects <- function(reference, comparison) {
  projects <- names(reference)
  changed <- projects[reference != comparison[projects]]
  if (length(changed) == 0) "none" else paste(changed, collapse = ";")
}

outlier_status <- function(clusters, project) {
  if (!project %in% names(clusters)) return(NA_character_)
  if (sum(clusters == clusters[[project]]) == 1) "singleton" else "grouped"
}

set_preserved <- function(clusters, project_set) {
  present <- intersect(project_set, names(clusters))
  if (length(present) < 2) return(NA)
  ref_cluster <- clusters[[present[[1]]]]
  all(clusters[present] == ref_cluster)
}

feature_class_sensitivity <- function(scaled_matrix, metadata, reference_clusters) {
  classes <- vapply(colnames(scaled_matrix), feature_class, character(1))
  subset_defs <- list(
    all_features = names(classes),
    mutation_only = names(classes)[classes == "mutation"],
    copy_number_aneuploidy_only = names(classes)[classes == "copy_number_aneuploidy"],
    expression_only = names(classes)[classes == "expression"],
    mutation_plus_copy_number = names(classes)[classes %in% c("mutation", "copy_number_aneuploidy")],
    mutation_plus_expression = names(classes)[classes %in% c("mutation", "expression")],
    copy_number_plus_expression = names(classes)[classes %in% c("copy_number_aneuploidy", "expression")]
  )
  subset_defs <- subset_defs[vapply(subset_defs, length, integer(1)) > 0]

  summaries <- list()
  assignment_rows <- list()
  for (name in names(subset_defs)) {
    features <- subset_defs[[name]]
    clusters <- cluster_k6_from_matrix(scaled_matrix[, features, drop = FALSE])
    ari <- adjusted_rand_index(reference_clusters[names(clusters)], clusters)
    changed <- changed_projects(reference_clusters, clusters)
    summaries[[length(summaries) + 1L]] <- data.table::data.table(
      feature_class_subset = name,
      n_features_used = length(features),
      distance_method_used = "gower_average_linkage_k6",
      adjusted_rand_index_vs_all_feature_gower_k6 = ari,
      major_clustering_differences = changed,
      dlbc_status = outlier_status(clusters, "DLBC"),
      skcm_status = outlier_status(clusters, "SKCM"),
      ucec_status = outlier_status(clusters, "UCEC"),
      ucs_status = outlier_status(clusters, "UCS"),
      carcinoma_grouping_preserved = set_preserved(clusters, c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "HNSC", "LUAD", "LUSC", "OV", "READ", "STAD")),
      kidney_cns_endocrine_low_mutation_cluster_preserved = set_preserved(clusters, c("KICH", "KIRC", "KIRP", "GBM", "LGG", "ACC", "PCPG", "THCA", "PRAD", "UVM"))
    )
    assignment_rows[[length(assignment_rows) + 1L]] <- data.table::data.table(
      project_code = names(clusters),
      feature_class_subset = name,
      n_features_used = length(features),
      cluster_k6 = as.integer(clusters)
    )
  }
  assignments <- data.table::rbindlist(assignment_rows)
  assignments <- merge(metadata[, .(project_id, project_code, major_cancer_group)], assignments, by = "project_code", all.y = TRUE, sort = FALSE)
  list(summary = data.table::rbindlist(summaries), assignments = assignments)
}

weighting_sensitivity <- function(scaled_matrix, metadata, reference_clusters) {
  classes <- vapply(colnames(scaled_matrix), feature_class, character(1))
  class_counts <- table(classes)
  expression_exists <- any(classes == "expression")
  scenario_weights <- list(unweighted_all_features = rep(1, length(classes)))
  names(scenario_weights$unweighted_all_features) <- names(classes)

  equal_weights <- rep(1, length(classes))
  for (class in names(class_counts)) {
    equal_weights[classes == class] <- 1 / as.numeric(class_counts[[class]])
  }
  scenario_weights$equal_weight_mutation_copy_number_expression <- equal_weights

  mutation_down <- rep(1, length(classes)); mutation_down[classes == "mutation"] <- 0.5
  copy_down <- rep(1, length(classes)); copy_down[classes == "copy_number_aneuploidy"] <- 0.5
  scenario_weights$mutation_down_weighted <- mutation_down
  scenario_weights$copy_number_down_weighted <- copy_down
  if (expression_exists) {
    expression_down <- rep(1, length(classes)); expression_down[classes == "expression"] <- 0.5
    scenario_weights$expression_down_weighted <- expression_down
  }

  summaries <- list()
  assignment_rows <- list()
  for (scenario in names(scenario_weights)) {
    weights <- scenario_weights[[scenario]]
    clusters <- cluster_k6_from_matrix(scaled_matrix, weights = weights)
    ari <- adjusted_rand_index(reference_clusters[names(clusters)], clusters)
    changed <- changed_projects(reference_clusters, clusters)
    summaries[[length(summaries) + 1L]] <- data.table::data.table(
      weighting_scenario = scenario,
      adjusted_rand_index_vs_unweighted_gower_k6 = ari,
      projects_changing_cluster = changed,
      n_projects_changing_cluster = ifelse(changed == "none", 0L, length(strsplit(changed, ";", fixed = TRUE)[[1]])),
      dlbc_status = outlier_status(clusters, "DLBC"),
      skcm_status = outlier_status(clusters, "SKCM"),
      ucec_status = outlier_status(clusters, "UCEC"),
      ucs_status = outlier_status(clusters, "UCS"),
      major_groups_stable = ari >= 0.75
    )
    assignment_rows[[length(assignment_rows) + 1L]] <- data.table::data.table(
      project_code = names(clusters),
      weighting_scenario = scenario,
      cluster_k6 = as.integer(clusters)
    )
  }
  assignments <- data.table::rbindlist(assignment_rows)
  assignments <- merge(metadata[, .(project_id, project_code, major_cancer_group)], assignments, by = "project_code", all.y = TRUE, sort = FALSE)
  list(summary = data.table::rbindlist(summaries), assignments = assignments)
}

bootstrap_interpretation <- function(bootstrap, metadata) {
  stable <- bootstrap[co_clustering_frequency >= 0.90]
  data.table::setorder(stable, -co_clustering_frequency, project_code_a, project_code_b)
  project_stats_a <- bootstrap[, .(
    mean_coclustering_frequency = mean(co_clustering_frequency, na.rm = TRUE),
    max_coclustering_frequency = max(co_clustering_frequency, na.rm = TRUE),
    n_high_stability_partners = sum(co_clustering_frequency >= 0.90, na.rm = TRUE)
  ), by = .(project_code = project_code_a)]
  project_stats_b <- bootstrap[, .(
    mean_coclustering_frequency = mean(co_clustering_frequency, na.rm = TRUE),
    max_coclustering_frequency = max(co_clustering_frequency, na.rm = TRUE),
    n_high_stability_partners = sum(co_clustering_frequency >= 0.90, na.rm = TRUE)
  ), by = .(project_code = project_code_b)]
  project_stats <- data.table::rbindlist(list(project_stats_a, project_stats_b))
  project_stats <- project_stats[, .(
    mean_coclustering_frequency = mean(mean_coclustering_frequency, na.rm = TRUE),
    max_coclustering_frequency = max(max_coclustering_frequency, na.rm = TRUE),
    n_high_stability_partners = sum(n_high_stability_partners, na.rm = TRUE)
  ), by = project_code]
  project_stats[, stability_label := fifelse(n_high_stability_partners == 0 | mean_coclustering_frequency < 0.20, "unstable_or_switching", fifelse(n_high_stability_partners >= 3, "stable_group_member", "limited_stable_partners"))]
  project_stats <- merge(metadata[, .(project_id, project_code, major_cancer_group)], project_stats, by = "project_code", all.y = TRUE, sort = FALSE)
  list(stable_pairs = stable, project_stats = project_stats)
}

plot_cluster_trait_summary <- function(cluster_summary, path) {
  ensure_dir(dirname(path))
  data <- cluster_summary[distance_method == "gower" & k == 6]
  traits <- c(
    "median_mutation_count_proxy", "median_TP53_prevalence", "median_KRAS_prevalence",
    "median_BRAF_prevalence", "median_purity", "median_ploidy", "median_aneuploidy_score",
    "median_immune_inflammatory_score", "median_proliferation_score", "median_EMT_score"
  )
  traits <- intersect(traits, names(data))
  mat <- as.matrix(data[, ..traits])
  rownames(mat) <- paste0("cluster_", data$cluster_id, "_", data$interpretation_label)
  mat <- scale(mat)
  mat[is.na(mat)] <- 0
  pheatmap::pheatmap(mat, filename = path, width = 11, height = 7, main = "Gower k=6 cluster trait summary", border_color = NA)
}

plot_feature_driver_barplot <- function(association, path) {
  ensure_dir(dirname(path))
  top <- head(association[order(fdr, -epsilon_squared)], 20)
  top[, feature_name := factor(feature_name, levels = rev(feature_name))]
  plot <- ggplot2::ggplot(top, ggplot2::aes(x = feature_name, y = epsilon_squared, fill = feature_class)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Top Level 1 molecular trait drivers of Gower k=6 clusters",
      x = NULL,
      y = "Kruskal-Wallis epsilon-squared",
      fill = "Feature class"
    ) +
    ggplot2::theme_bw(base_size = 11)
  ggplot2::ggsave(path, plot, width = 11, height = 8.5, device = "pdf")
}

plot_bootstrap_heatmap <- function(bootstrap, path) {
  ensure_dir(dirname(path))
  projects <- sort(unique(c(bootstrap$project_code_a, bootstrap$project_code_b)))
  mat <- matrix(0, nrow = length(projects), ncol = length(projects), dimnames = list(projects, projects))
  for (i in seq_len(nrow(bootstrap))) {
    a <- bootstrap$project_code_a[[i]]
    b <- bootstrap$project_code_b[[i]]
    mat[a, b] <- bootstrap$co_clustering_frequency[[i]]
    mat[b, a] <- bootstrap$co_clustering_frequency[[i]]
  }
  diag(mat) <- 1
  pheatmap::pheatmap(mat, filename = path, width = 10, height = 10, main = "Feature-bootstrap co-clustering frequency", border_color = NA)
}

plot_sensitivity_heatmap <- function(assignments, scenario_col, path, title) {
  ensure_dir(dirname(path))
  formula <- stats::as.formula(sprintf("project_code ~ %s", scenario_col))
  wide <- data.table::dcast(assignments, formula, value.var = "cluster_k6")
  project_codes <- wide$project_code
  mat <- as.matrix(wide[, -"project_code"])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- project_codes
  pheatmap::pheatmap(mat, filename = path, width = 9, height = 10, main = title, cluster_cols = FALSE, border_color = NA)
}

plot_weighting_mds <- function(assignments, path) {
  ensure_dir(dirname(path))
  scenarios <- unique(assignments$weighting_scenario)
  projects <- unique(assignments$project_code)
  mat <- matrix(0, nrow = length(projects), ncol = length(scenarios), dimnames = list(projects, scenarios))
  for (scenario in scenarios) {
    vals <- assignments[weighting_scenario == scenario]
    mat[vals$project_code, scenario] <- vals$cluster_k6
  }
  distance <- stats::dist(mat)
  coords <- as.data.frame(stats::cmdscale(distance, k = 2))
  names(coords) <- c("MDS1", "MDS2")
  coords$project_code <- rownames(coords)
  plot <- ggplot2::ggplot(coords, ggplot2::aes(MDS1, MDS2, label = project_code)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_text(vjust = -0.6, size = 3) +
    ggplot2::labs(title = "Project movement across feature-weighting sensitivity scenarios", x = "MDS1", y = "MDS2") +
    ggplot2::theme_bw(base_size = 11)
  ggplot2::ggsave(path, plot, width = 10, height = 8, device = "pdf")
}

top_n_text <- function(values, n = 8) {
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) return("none")
  paste(head(values, n), collapse = ", ")
}

write_interpretation_report <- function(path, gower_summary, correlation_summary, association, class_summary, weighting_summary, bootstrap_info) {
  ensure_dir(dirname(path))
  gower_k6 <- gower_summary[k == 6]
  correlation_k6 <- correlation_summary[k == 6]
  top_features <- head(association[order(fdr, -epsilon_squared), feature_name], 10)
  stable_pairs <- head(bootstrap_info$stable_pairs, 10)
  unstable <- bootstrap_info$project_stats[stability_label == "unstable_or_switching", project_code]
  class_best <- class_summary[order(-adjusted_rand_index_vs_all_feature_gower_k6)]
  weighting_changes <- weighting_summary[weighting_scenario != "unweighted_all_features"]

  lines <- c(
    "# Level 1 Interpretation and Sensitivity Summary",
    "",
    "This analysis interprets the Level 1 TCGA pan-cancer molecular similarity dendrograms by cutting the Gower and correlation average-linkage trees at k = 3 through k = 8, summarizing molecular traits within clusters, testing feature associations with the Gower k = 6 clustering, and comparing feature-class and weighting sensitivity runs. Cluster stability is summarized from feature-bootstrap project-pair co-clustering frequencies.",
    "",
    sprintf("The Gower k = 6 tree separates projects into clusters labeled as: %s. These groups should be read as projects that cluster with similar scaled molecular trait values, not as evolutionary lineages.", paste(sprintf("cluster %s (%s: %s)", gower_k6$cluster_id, gower_k6$interpretation_label, gower_k6$projects_included), collapse = "; ")),
    "",
    sprintf("The correlation k = 6 tree separates projects into clusters labeled as: %s. Correlation distance emphasizes similarity in each project's relative feature-profile shape across traits rather than absolute scaled levels.", paste(sprintf("cluster %s (%s: %s)", correlation_k6$cluster_id, correlation_k6$interpretation_label, correlation_k6$projects_included), collapse = "; ")),
    "",
    "The Gower and correlation trees differ because Gower distance is sensitive to absolute trait differences after feature-level normalization, whereas correlation distance can cluster projects that have parallel molecular profiles even when their overall mutation burden, copy-number burden, or expression-program magnitude differs.",
    "",
    sprintf("Feature-driver analysis suggests that the strongest Gower k = 6 associations include: %s. These associations are descriptive Kruskal-Wallis tests across tree-derived clusters, so they indicate traits that separate the observed clusters rather than causal drivers.", top_n_text(top_features, 10)),
    "",
    sprintf("Feature-class sensitivity shows the closest agreement with the all-feature reference for: %s. Weighting sensitivity changed clusters for %s across non-reference scenarios, indicating that the broad interpretation should emphasize robust groupings over exact branch placement.", paste(sprintf("%s (ARI %.2f)", class_best$feature_class_subset, class_best$adjusted_rand_index_vs_all_feature_gower_k6), collapse = "; "), paste(sprintf("%s: %s", weighting_changes$weighting_scenario, weighting_changes$projects_changing_cluster), collapse = "; ")),
    "",
    sprintf("Bootstrap co-clustering highlights stable pairs such as: %s. Projects with relatively low or diffuse co-clustering support include: %s.", paste(sprintf("%s-%s (%.2f)", stable_pairs$project_code_a, stable_pairs$project_code_b, stable_pairs$co_clustering_frequency), collapse = "; "), top_n_text(unstable, 12)),
    "",
    "Limitations: the current tree uses project-level aggregated traits and first-pass pathway scores, mutation burden is still a mutation-count proxy rather than true mutations per megabase, and the feature set mixes mutation, copy-number, purity/ploidy, and expression summaries with different biological meanings. These results should be treated as molecular similarity analysis and hypothesis generation, not as a literal species-like cancer phylogeny, and they should not be used to infer directional ancestry between cancer types.",
    ""
  )
  writeLines(lines, con = path)
  log_info(sprintf("Wrote %s", path))
}

build_level1_interpretation_outputs <- function(raw_matrix, scaled_table, gower_distance, correlation_distance, bootstrap, projects, root = find_project_root()) {
  raw_matrix <- ensure_metadata(raw_matrix, projects)
  scaled_table <- ensure_metadata(scaled_table, projects)
  scaled_matrix <- matrix_from_scaled(scaled_table)
  raw_features <- feature_columns(raw_matrix)
  for (feature in raw_features) raw_matrix[, (feature) := numeric_or_na(get(feature))]

  metadata <- scaled_table[, .(project_id, project_code, major_cancer_group, disease_type, primary_site)]
  gower_clusters <- cluster_assignments_from_distance(gower_distance, metadata)
  correlation_clusters <- cluster_assignments_from_distance(correlation_distance, metadata)
  gower_summary <- cluster_summaries(gower_clusters$assignments, raw_matrix, "gower")
  correlation_summary <- cluster_summaries(correlation_clusters$assignments, raw_matrix, "correlation")

  pca_loadings <- feature_pca_loadings(scaled_matrix)
  variance_summary <- feature_variance_summary(raw_matrix, scaled_matrix)
  association <- feature_cluster_association(raw_matrix, scaled_matrix, gower_clusters$assignments)
  reference_clusters <- gower_clusters$assignments$cluster_k6
  names(reference_clusters) <- gower_clusters$assignments$project_code
  class_sensitivity <- feature_class_sensitivity(scaled_matrix, metadata, reference_clusters)
  weighting <- weighting_sensitivity(scaled_matrix, metadata, reference_clusters)
  bootstrap_info <- bootstrap_interpretation(bootstrap, metadata)

  table_dir <- project_path("results", "tables", root = root)
  figure_dir <- project_path("results", "figures", root = root)
  report_dir <- project_path("results", "reports", root = root)

  paths <- list(
    gower_assignments = file.path(table_dir, "level1_gower_cluster_assignments.tsv"),
    correlation_assignments = file.path(table_dir, "level1_correlation_cluster_assignments.tsv"),
    gower_summaries = file.path(table_dir, "level1_gower_cluster_summaries.tsv"),
    correlation_summaries = file.path(table_dir, "level1_correlation_cluster_summaries.tsv"),
    pca_loadings = file.path(table_dir, "level1_feature_pca_loadings.tsv"),
    variance_summary = file.path(table_dir, "level1_feature_variance_summary.tsv"),
    association = file.path(table_dir, "level1_feature_cluster_association.tsv"),
    class_summary = file.path(table_dir, "level1_feature_class_sensitivity_summary.tsv"),
    class_assignments = file.path(table_dir, "level1_feature_class_cluster_assignments.tsv"),
    weighting_summary = file.path(table_dir, "level1_feature_weighting_sensitivity_summary.tsv"),
    weighting_assignments = file.path(table_dir, "level1_weighted_cluster_assignments.tsv"),
    stable_pairs = file.path(table_dir, "level1_stable_project_pairs.tsv"),
    bootstrap_project_stability = file.path(table_dir, "level1_bootstrap_project_stability.tsv"),
    report = file.path(report_dir, "level1_interpretation_summary.md"),
    cluster_trait_heatmap = file.path(figure_dir, "level1_cluster_trait_summary_heatmap.pdf"),
    feature_driver_barplot = file.path(figure_dir, "level1_feature_driver_barplot.pdf"),
    bootstrap_heatmap = file.path(figure_dir, "level1_bootstrap_coclustering_heatmap.pdf"),
    class_heatmap = file.path(figure_dir, "level1_feature_class_sensitivity_heatmap.pdf"),
    weighting_mds = file.path(figure_dir, "level1_weighting_sensitivity_pca_or_mds.pdf")
  )

  write_tsv_table(gower_clusters$assignments, paths$gower_assignments)
  write_tsv_table(correlation_clusters$assignments, paths$correlation_assignments)
  write_tsv_table(gower_summary, paths$gower_summaries)
  write_tsv_table(correlation_summary, paths$correlation_summaries)
  write_tsv_table(pca_loadings, paths$pca_loadings)
  write_tsv_table(variance_summary, paths$variance_summary)
  write_tsv_table(association, paths$association)
  write_tsv_table(class_sensitivity$summary, paths$class_summary)
  write_tsv_table(class_sensitivity$assignments, paths$class_assignments)
  write_tsv_table(weighting$summary, paths$weighting_summary)
  write_tsv_table(weighting$assignments, paths$weighting_assignments)
  write_tsv_table(bootstrap_info$stable_pairs, paths$stable_pairs)
  write_tsv_table(bootstrap_info$project_stats, paths$bootstrap_project_stability)

  all_summaries <- data.table::rbindlist(list(gower_summary, correlation_summary), use.names = TRUE, fill = TRUE)
  plot_cluster_trait_summary(all_summaries, paths$cluster_trait_heatmap)
  plot_feature_driver_barplot(association, paths$feature_driver_barplot)
  plot_bootstrap_heatmap(bootstrap, paths$bootstrap_heatmap)
  plot_sensitivity_heatmap(class_sensitivity$assignments, "feature_class_subset", paths$class_heatmap, "Feature-class sensitivity cluster assignments")
  plot_weighting_mds(weighting$assignments, paths$weighting_mds)
  write_interpretation_report(paths$report, gower_summary, correlation_summary, association, class_sensitivity$summary, weighting$summary, bootstrap_info)

  list(
    gower_assignments = gower_clusters$assignments,
    correlation_assignments = correlation_clusters$assignments,
    gower_summary = gower_summary,
    correlation_summary = correlation_summary,
    association = association,
    class_sensitivity = class_sensitivity,
    weighting = weighting,
    bootstrap_info = bootstrap_info,
    paths = paths
  )
}

main <- function() {
  require_interpretation_packages()
  root <- find_project_root()
  raw_matrix <- read_tsv_required(project_path("results", "tables", "level1_feature_matrix.tsv", root = root), "Level 1 feature matrix")
  scaled_table <- read_tsv_required(project_path("results", "tables", "level1_feature_matrix_scaled.tsv", root = root), "scaled Level 1 feature matrix")
  gower_distance <- read_distance_matrix(project_path("results", "tables", "level1_gower_distance_matrix.tsv", root = root))
  correlation_distance <- read_distance_matrix(project_path("results", "tables", "level1_correlation_distance_matrix.tsv", root = root))
  bootstrap <- read_tsv_required(project_path("results", "tables", "level1_bootstrap_cluster_stability.tsv", root = root), "Level 1 bootstrap cluster stability table")
  read_tsv_required(project_path("results", "tables", "level1_tree_qc_summary.tsv", root = root), "Level 1 tree QC summary")
  require_file(project_path("results", "trees", "level1_pan_cancer_gower_hclust_tree.nwk", root = root), "Gower Newick tree")
  require_file(project_path("results", "trees", "level1_pan_cancer_correlation_hclust_tree.nwk", root = root), "correlation Newick tree")
  require_file(project_path("config", "cancer_group_map.csv", root = root), "cancer group map")
  projects <- read_tsv_required(project_path("data", "interim", "tcga_projects.tsv", root = root), "TCGA project metadata")

  log_info("Interpreting Level 1 pan-cancer molecular similarity trees and sensitivity analyses")
  outputs <- build_level1_interpretation_outputs(raw_matrix, scaled_table, gower_distance, correlation_distance, bootstrap, projects, root = root)
  log_info(sprintf(
    "Level 1 interpretation complete: %d Gower k=6 clusters, %d top-level feature associations, %d stable bootstrap pairs",
    length(unique(outputs$gower_assignments$cluster_k6)),
    nrow(outputs$association),
    nrow(outputs$bootstrap_info$stable_pairs)
  ))
  invisible(outputs)
}

if (identical(environment(), globalenv())) {
  main()
}
