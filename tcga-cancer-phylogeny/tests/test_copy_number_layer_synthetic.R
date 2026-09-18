#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "03_download_copy_number.R"), envir = env)
env$require_copy_number_packages()

tmp <- tempfile("copy_number_layer_synthetic_")
dir.create(file.path(tmp, "data", "raw", "purity_ploidy"), recursive = TRUE)
dir.create(file.path(tmp, "data", "raw", "copy_number"), recursive = TRUE)

projects <- data.table::data.table(
  project_id = c("TCGA-BRCA", "TCGA-LUAD"),
  project_code = c("BRCA", "LUAD")
)
sample_universe <- data.table::data.table(
  sample_barcode = c("TCGA-AB-1234-01A", "TCGA-CD-5678-01A", "TCGA-EF-9999-01A"),
  patient_barcode = c("TCGA-AB-1234", "TCGA-CD-5678", "TCGA-EF-9999"),
  project_id = c("TCGA-BRCA", "TCGA-BRCA", "TCGA-LUAD"),
  project_code = c("BRCA", "BRCA", "LUAD")
)

absolute_path <- file.path(tmp, "data", "raw", "purity_ploidy", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")
arm_path <- file.path(tmp, "data", "raw", "copy_number", "PANCAN_ArmCallsAndAneuploidyScore_092817.txt")

absolute <- data.frame(
  array = c("TCGA-AB-1234-01", "TCGA-CD-5678-01"),
  sample = c("TCGA-AB-1234-01A-01D-0000-01", "TCGA-CD-5678-01A-01D-0000-01"),
  `call status` = c("called", "called"),
  purity = c(0.70, 0.50),
  ploidy = c(2.8, 2.0),
  `Genome doublings` = c(1, 0),
  `Coverage for 80% power` = c(10, 8),
  `Cancer DNA fraction` = c(0.65, 0.45),
  `Subclonal genome fraction` = c(0.20, 0.10),
  solution = c("new", "new"),
  check.names = FALSE
)
arm_calls <- data.frame(
  Sample = c("TCGA-AB-1234-01", "TCGA-CD-5678-01", "TCGA-EF-9999-01"),
  Type = c("BRCA", "BRCA", "LUAD"),
  `Aneuploidy Score` = c(2, 1, 3),
  `1p` = c(-1, 0, 1),
  `1q` = c(1, 1, 1),
  `2p` = c(0, 0, -1),
  check.names = FALSE
)

write.table(absolute, absolute_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(arm_calls, arm_path, sep = "\t", quote = FALSE, row.names = FALSE)

retrieval <- env$retrieve_pancanatlas_resources(tmp, refresh = FALSE)
stopifnot(all(retrieval$status == "local"))

resources <- env$discover_local_resources(tmp)
features <- env$build_copy_number_features(tmp, projects, sample_universe, resources)
project_summary <- env$summarize_project_features(features$purity, features$aneuploidy, projects)

stopifnot(nrow(features$purity) == 3)
stopifnot(features$purity[sample_barcode == "TCGA-AB-1234-01A", purity] == 0.70)
stopifnot(features$purity[sample_barcode == "TCGA-CD-5678-01A", ploidy] == 2.0)
stopifnot(features$purity[sample_barcode == "TCGA-AB-1234-01A", whole_genome_doubling_status] == "WGD")
stopifnot(is.na(features$purity[sample_barcode == "TCGA-EF-9999-01A", purity]))

brca_a <- features$aneuploidy[sample_barcode == "TCGA-AB-1234-01A"]
stopifnot(brca_a$aneuploidy_score == 2)
stopifnot(brca_a$arm_gain_count == 1)
stopifnot(brca_a$arm_loss_count == 1)
stopifnot(brca_a$total_arm_alteration_count == 2)

luad <- features$aneuploidy[sample_barcode == "TCGA-EF-9999-01A"]
stopifnot(luad$aneuploidy_score == 3)
stopifnot(luad$arm_gain_count == 2)
stopifnot(luad$arm_loss_count == 1)
stopifnot(nrow(features$arm_calls) == 9)

brca_summary <- project_summary[project_id == "TCGA-BRCA"]
stopifnot(abs(brca_summary$median_purity - 0.60) < 1e-8)
stopifnot(abs(brca_summary$median_ploidy - 2.4) < 1e-8)
stopifnot(brca_summary$n_samples == 2)
stopifnot(abs(brca_summary$median_aneuploidy_score - 1.5) < 1e-8)

fallback_tmp <- tempfile("copy_number_layer_fallback_")
dir.create(file.path(fallback_tmp, "data", "raw", "purity_ploidy"), recursive = TRUE)
dir.create(file.path(fallback_tmp, "data", "raw", "copy_number"), recursive = TRUE)
fallback_features <- env$build_copy_number_features(fallback_tmp, projects, sample_universe, character())
stopifnot(all(is.na(fallback_features$purity$purity)))
stopifnot(all(is.na(fallback_features$aneuploidy$aneuploidy_score)))
stopifnot(all(fallback_features$purity$source == "no_source_available"))

cat("Synthetic copy-number-layer checks passed\n")
