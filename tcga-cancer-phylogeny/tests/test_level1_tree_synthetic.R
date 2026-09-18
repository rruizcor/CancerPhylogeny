#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "07_level1_pan_cancer_tree.R"), envir = env)
env$require_tree_packages()

fail <- function(message) {
  stop(message, call. = FALSE)
}

feature_matrix <- data.table::data.table(
  project_id = c("TCGA-A", "TCGA-B", "TCGA-C", "TCGA-D"),
  project_code = c("A", "B", "C", "D"),
  cancer_name = c("Cancer A", "Cancer B", "Cancer C", "Cancer D"),
  broad_group = c("Carcinoma", "Carcinoma", "Sarcoma", "Sarcoma"),
  level2_group = c("a", "b", "c", "d"),
  lineage_notes = c("note", "note", "note", "note"),
  source = c("synthetic", "synthetic", "synthetic", "synthetic"),
  tmb_proxy_median_nonsynonymous_count = c(10, 12, 100, 105),
  mutation_prevalence_TP53 = c(0.2, NA_real_, 0.7, 0.75),
  median_purity = c(0.6, 0.65, 0.5, 0.52),
  median_aneuploidy_score = c(4, 5, 12, 14),
  mostly_missing_feature = c(NA_real_, NA_real_, NA_real_, 1),
  constant_feature = c(3, 3, 3, 3)
)

selected <- env$select_biologic_feature_columns(feature_matrix)
if (!"tmb_proxy_median_nonsynonymous_count" %in% selected) fail("Numeric biologic feature was not selected")
if ("project_code" %in% selected || "broad_group" %in% selected || "source" %in% selected) {
  fail("Annotation/source columns were selected as biologic features")
}

prepared <- env$prepare_tree_features(feature_matrix, max_missing_fraction = 0.30)
if (!"mutation_prevalence_TP53" %in% prepared$feature_columns_retained) fail("Feature with allowable missingness was not retained")
if ("mostly_missing_feature" %in% prepared$feature_columns_retained) fail("High-missingness feature was retained")
if ("constant_feature" %in% prepared$feature_columns_retained) fail("Zero-variance feature was retained")
if (prepared$feature_qc[feature_name == "mutation_prevalence_TP53", n_imputed_values] != 1L) {
  fail("Median imputation count was not recorded")
}
imputed_matrix <- prepared$raw_feature_matrix
expected_median <- stats::median(c(0.2, 0.7, 0.75))
if (!isTRUE(all.equal(imputed_matrix["B", "mutation_prevalence_TP53"], expected_median))) {
  fail("Median imputation used an unexpected value")
}

scaled <- env$scale_tree_feature_matrix(prepared)
if (!is.matrix(scaled$scaled_matrix) || !is.numeric(scaled$scaled_matrix)) fail("Scaled matrix is not numeric")
if (any(is.na(scaled$scaled_matrix))) fail("Scaled matrix contains NA values")

gower <- env$compute_gower_distance_matrix(scaled$scaled_matrix)
correlation <- env$compute_correlation_distance_matrix(scaled$scaled_matrix)
if (nrow(gower) != ncol(gower) || nrow(gower) != 4) fail("Gower distance matrix is not square")
if (nrow(correlation) != ncol(correlation) || nrow(correlation) != 4) fail("Correlation distance matrix is not square")
if (max(abs(gower - t(gower))) > 1e-8) fail("Gower distance matrix is not symmetric")
if (max(abs(correlation - t(correlation))) > 1e-8) fail("Correlation distance matrix is not symmetric")

tree <- env$build_hclust_tree(gower)
newick_path <- tempfile(fileext = ".nwk")
env$write_hclust_newick(tree, newick_path)
if (!file.exists(newick_path) || file.info(newick_path)$size == 0) fail("Synthetic Newick output was empty")

stability <- env$bootstrap_cluster_stability(scaled$scaled_matrix, iterations = 5L, seed = 42L, k = 2L)
if (nrow(stability) != choose(4, 2)) fail("Bootstrap pairwise stability table has unexpected row count")
if (!all(stability$bootstrap_iterations == 5L)) fail("Bootstrap iteration count was not recorded")

cat("Synthetic Level 1 tree checks passed\n")
