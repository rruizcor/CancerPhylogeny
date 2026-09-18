#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "07b_interpret_level1_tree.R"), envir = env)
env$require_interpretation_packages()

fail <- function(message) {
  stop(message, call. = FALSE)
}

projects <- data.table::data.table(
  project_id = paste0("TCGA-", LETTERS[1:8]),
  project_code = LETTERS[1:8],
  disease_type = rep(c("Epithelial", "Mesenchymal"), each = 4),
  primary_site = paste("Site", LETTERS[1:8])
)
raw <- data.table::data.table(
  project_id = projects$project_id,
  project_code = projects$project_code,
  cancer_name = paste("Cancer", LETTERS[1:8]),
  broad_group = rep(c("Carcinoma", "Sarcoma"), each = 4),
  level2_group = letters[1:8],
  lineage_notes = "synthetic",
  tmb_proxy_median_nonsynonymous_count = c(10, 11, 12, 13, 80, 82, 85, 90),
  mutation_prevalence_TP53 = c(0.1, 0.1, 0.2, 0.2, 0.7, 0.7, 0.8, 0.8),
  mutation_prevalence_KRAS = c(0.01, 0.02, 0.03, 0.04, 0.2, 0.22, 0.25, 0.28),
  median_purity = c(0.6, 0.61, 0.62, 0.63, 0.55, 0.54, 0.53, 0.52),
  median_ploidy = c(2.1, 2.1, 2.2, 2.2, 3.5, 3.6, 3.7, 3.8),
  median_aneuploidy_score = c(4, 4, 5, 5, 12, 13, 14, 15),
  median_proliferation_score = c(-1, -0.9, -0.8, -0.7, 0.7, 0.8, 0.9, 1.0),
  median_immune_inflammatory_score = c(-0.5, -0.4, -0.3, -0.2, 0.2, 0.3, 0.4, 0.5),
  median_EMT_score = c(-0.2, -0.1, 0, 0.1, 0.3, 0.4, 0.5, 0.6)
)

raw <- env$ensure_metadata(raw, projects)
features <- env$feature_columns(raw)
scaled <- data.table::copy(raw)
scaled_matrix <- scale(as.matrix(raw[, ..features]))
scaled_matrix[is.na(scaled_matrix)] <- 0
for (feature in features) scaled[, (feature) := scaled_matrix[, feature]]
scaled_matrix <- env$matrix_from_scaled(scaled)

distance <- env$compute_gower_distance(scaled_matrix)
clusters <- env$cluster_assignments_from_distance(distance, raw, k_values = 3:8)
required_cluster_cols <- paste0("cluster_k", 3:8)
if (!all(required_cluster_cols %in% names(clusters$assignments))) fail("Cluster assignments for k=3 through k=8 were not created")
if (nrow(clusters$assignments) != 8) fail("Synthetic assignments do not include all projects")

summary <- env$cluster_summaries(clusters$assignments, raw, "synthetic", k_values = 3:8)
if (nrow(summary) == 0) fail("Synthetic cluster summaries are empty")
if (!"interpretation_label" %in% names(summary)) fail("Synthetic cluster summaries lack interpretation labels")

association <- env$feature_cluster_association(raw, scaled_matrix, clusters$assignments, cluster_col = "cluster_k6")
if (nrow(association) == 0 || !"fdr" %in% names(association)) fail("Synthetic feature association table is empty or lacks FDR")

reference <- clusters$assignments$cluster_k6
names(reference) <- clusters$assignments$project_code
class_sensitivity <- env$feature_class_sensitivity(scaled_matrix, raw[, .(project_id, project_code, major_cancer_group)], reference)
if (nrow(class_sensitivity$summary) == 0) fail("Synthetic feature-class sensitivity summary is empty")
if (any(is.na(class_sensitivity$summary$adjusted_rand_index_vs_all_feature_gower_k6))) fail("Synthetic ARI calculation returned NA")

weighting <- env$weighting_sensitivity(scaled_matrix, raw[, .(project_id, project_code, major_cancer_group)], reference)
if (nrow(weighting$summary) == 0) fail("Synthetic weighting sensitivity summary is empty")
if (!isTRUE(all.equal(env$adjusted_rand_index(reference, reference), 1))) fail("ARI identity check failed")

bootstrap <- data.table::data.table(
  project_code_a = rep(LETTERS[1:4], each = 2),
  project_code_b = rep(LETTERS[5:6], times = 4),
  bootstrap_iterations = 10,
  cutree_k = 6,
  n_same_cluster = c(9, 1, 8, 2, 10, 0, 7, 3),
  co_clustering_frequency = c(0.9, 0.1, 0.8, 0.2, 1, 0, 0.7, 0.3)
)
boot <- env$bootstrap_interpretation(bootstrap, raw[, .(project_id, project_code, major_cancer_group)])
if (nrow(boot$stable_pairs) == 0) fail("Synthetic stable bootstrap pairs were not identified")

cat("Synthetic Level 1 interpretation checks passed\n")
