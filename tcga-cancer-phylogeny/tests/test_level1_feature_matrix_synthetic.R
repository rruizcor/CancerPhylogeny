#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "06_build_pan_cancer_features.R"), envir = env)
env$require_feature_packages()

projects <- data.table::data.table(
  project_id = c("TCGA-BRCA", "TCGA-LUAD"),
  project_code = c("BRCA", "LUAD")
)
group_map <- data.table::data.table(
  project_id = c("TCGA-BRCA", "TCGA-LUAD"),
  project_code = c("BRCA", "LUAD"),
  cancer_name = c("Breast Invasive Carcinoma", "Lung Adenocarcinoma"),
  broad_group = c("Carcinoma", "Carcinoma"),
  level2_group = c("breast", "lung_adeno"),
  lineage_notes = c("Synthetic breast", "Synthetic lung")
)
tmb <- data.table::data.table(
  sample_barcode = c("TCGA-AB-1234-01A", "TCGA-CD-5678-01A", "TCGA-EF-9999-01A"),
  project_id = c("TCGA-BRCA", "TCGA-BRCA", "TCGA-LUAD"),
  project_code = c("BRCA", "BRCA", "LUAD"),
  nonsynonymous_count = c(10, 20, 30),
  total_mutation_count = c(15, 25, 45)
)
drivers <- data.table::data.table(
  project_id = c("TCGA-BRCA", "TCGA-LUAD"),
  project_code = c("BRCA", "LUAD"),
  TP53 = c(0.30, 0.50),
  KRAS = c(0.02, 0.25)
)
copy_number <- data.table::data.table(
  project_id = c("TCGA-BRCA", "TCGA-LUAD"),
  project_code = c("BRCA", "LUAD"),
  n_samples = c(2, 1),
  median_purity = c(0.70, 0.55),
  median_ploidy = c(2.4, 3.1),
  median_aneuploidy_score = c(8, 12),
  median_fraction_genome_altered = c(NA_real_, NA_real_),
  median_arm_gain_count = c(3, 5),
  median_arm_loss_count = c(2, 7),
  median_total_arm_alteration_count = c(5, 12),
  median_copy_number_complexity_score = c(NA_real_, NA_real_),
  missingness_purity = c(0, 0),
  missingness_ploidy = c(0, 0),
  missingness_aneuploidy_score = c(0, 0),
  missingness_fraction_genome_altered = c(1, 1)
)
expression <- data.table::data.table(
  project_id = c("TCGA-BRCA", "TCGA-LUAD"),
  project_code = c("BRCA", "LUAD"),
  n_samples = c(2, 1),
  median_proliferation_score = c(0.1, 0.2),
  median_immune_inflammatory_score = c(-0.1, -0.3),
  median_interferon_gamma_score = c(-0.2, -0.1),
  median_cytotoxic_t_cell_score = c(-0.4, -0.2),
  median_stromal_score = c(0.3, 0.1),
  median_epithelial_score = c(0.5, 0.4),
  median_EMT_score = c(-0.2, 0.6),
  median_hypoxia_score = c(0.1, 0.7),
  median_cell_cycle_score = c(0.2, 0.4),
  median_DNA_repair_score = c(0.3, 0.5),
  median_angiogenesis_score = c(0.1, 0.2),
  median_lineage_or_tissue_PC1 = c(1.0, -1.0),
  median_lineage_or_tissue_PC2 = c(0.5, -0.5),
  missingness_proliferation_score = c(0, 0)
)

features <- env$build_level1_feature_matrix(
  projects,
  group_map,
  tmb,
  drivers,
  copy_number,
  expression,
  max_missing_fraction = 0.30
)

matrix <- features$feature_matrix
stopifnot(nrow(matrix) == 2)
stopifnot(all(c("project_id", "project_code", "broad_group", "mutation_prevalence_TP53") %in% names(matrix)))
stopifnot("median_purity" %in% names(matrix))
stopifnot("median_proliferation_score" %in% names(matrix))
stopifnot(!"median_fraction_genome_altered" %in% names(matrix))
stopifnot(!"median_copy_number_complexity_score" %in% names(matrix))
stopifnot(matrix[project_id == "TCGA-BRCA", tmb_proxy_median_nonsynonymous_count] == 15)
stopifnot(matrix[project_id == "TCGA-LUAD", mutation_prevalence_KRAS] == 0.25)
stopifnot(features$sample_support[project_id == "TCGA-BRCA", mutation_n_samples] == 2)
stopifnot(features$sample_support[project_id == "TCGA-LUAD", expression_n_samples] == 1)
stopifnot(features$feature_qc[feature_name == "median_fraction_genome_altered", retained] == FALSE)

env$validate_level1_outputs(matrix, projects)

cat("Synthetic Level 1 feature-matrix checks passed\n")
