#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "08_level2_within_group_trees.R"), envir = env)
env$require_level2_packages()

fail <- function(message) {
  stop(message, call. = FALSE)
}

codes <- c(
  "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP",
  "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PRAD", "READ", "STAD", "THCA", "UCEC",
  "UCS", "GBM", "LGG", "SARC"
)
project_ids <- paste0("TCGA-", codes)
projects <- data.table::data.table(
  project_id = project_ids,
  project_code = codes,
  disease_type = paste("Disease", codes),
  primary_site = paste("Site", codes)
)
group_map <- data.table::data.table(
  project_code = codes,
  project_id = project_ids,
  cancer_name = paste("Cancer", codes),
  broad_group = c(
    rep("Carcinoma", 20), "Other or mixed", "CNS/glial", "CNS/glial", "Sarcoma"
  ),
  level2_group = c(
    "urothelial", "breast", "gynecologic", "gi_hepatobiliary", "gi_colorectal", "gi_esophageal",
    "pan_squamous", "kidney", "kidney", "kidney", "gi_hepatobiliary", "lung_adeno",
    "pan_squamous", "gynecologic", "gi_pancreatic", "prostate", "gi_colorectal", "gi_gastric",
    "thyroid", "gynecologic", "gynecologic_mixed", "adult_glioma", "adult_glioma", "sarcoma"
  ),
  lineage_notes = "synthetic"
)
idx <- seq_along(codes)
feature_matrix <- data.table::data.table(
  project_id = project_ids,
  project_code = codes,
  cancer_name = paste("Cancer", codes),
  broad_group = group_map$broad_group,
  major_cancer_group = group_map$broad_group,
  level2_group = group_map$level2_group,
  lineage_notes = "synthetic",
  tmb_proxy_median_nonsynonymous_count = idx * 5,
  mutation_prevalence_TP53 = idx / 30,
  mutation_prevalence_KRAS = rev(idx) / 40,
  median_purity = 0.45 + idx / 100,
  median_ploidy = 2 + idx / 20,
  median_aneuploidy_score = idx %% 12,
  median_arm_gain_count = idx %% 7,
  median_arm_loss_count = rev(idx) %% 7,
  median_immune_inflammatory_score = scale(idx)[, 1],
  median_proliferation_score = scale(rev(idx))[, 1],
  median_EMT_score = sin(idx),
  constant_feature = 1
)
scaled_matrix <- data.table::copy(feature_matrix)

groups <- env$define_level2_groups(feature_matrix, group_map)
if (!all(c("carcinoma_all", "pan_squamous", "kidney", "cns_glial", "sarcoma") %in% groups$group_name)) {
  fail("Synthetic group definitions are missing expected groups")
}

prepared <- env$prepare_group_features(feature_matrix, c("BLCA", "BRCA", "CESC", "CHOL"))
if ("constant_feature" %in% prepared$retained_features) fail("Constant feature was retained")
if (!any(prepared$feature_qc$drop_reason == "constant_or_near_constant_within_group")) fail("Constant feature drop reason was not recorded")

tmp_root <- tempfile("level2_synthetic_")
dir.create(tmp_root, recursive = TRUE)
outputs <- env$build_level2_outputs(feature_matrix, scaled_matrix, group_map, projects, root = tmp_root, bootstrap_iterations = 5L, seed = 42L)
defs <- outputs$groups
if (!"carcinoma_all" %in% defs$group_name) fail("Synthetic output lacks carcinoma_all definition")
if (!isTRUE(defs[group_name == "carcinoma_all", tree_attempted])) fail("Synthetic carcinoma_all tree was not attempted")
if (isTRUE(defs[group_name == "cns_glial", tree_attempted])) fail("Synthetic cns_glial should be underpowered")
if (!grepl("underpowered_project_count", defs[group_name == "cns_glial", reason_if_not_attempted])) fail("Underpowered reason was not recorded")

carcinoma_tree <- file.path(tmp_root, "results", "trees", "level2_carcinoma_all_gower_hclust_tree.nwk")
if (!file.exists(carcinoma_tree) || file.info(carcinoma_tree)$size == 0) fail("Synthetic eligible Newick tree was not created")
distance_path <- file.path(tmp_root, "results", "tables", "level2_carcinoma_all_gower_distance_matrix.tsv")
distance_table <- read.delim(distance_path, stringsAsFactors = FALSE, check.names = FALSE)
distance_matrix <- as.matrix(distance_table[, setdiff(names(distance_table), "project_code"), drop = FALSE])
storage.mode(distance_matrix) <- "numeric"
if (nrow(distance_matrix) != ncol(distance_matrix)) fail("Synthetic distance matrix is not square")
if (max(abs(distance_matrix - t(distance_matrix))) > 1e-8) fail("Synthetic distance matrix is not symmetric")

cat("Synthetic Level 2 tree checks passed\n")
