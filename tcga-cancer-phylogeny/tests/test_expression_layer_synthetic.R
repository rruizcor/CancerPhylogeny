#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "04_download_expression.R"), envir = env)
env$require_expression_packages()

tmp <- tempfile("expression_layer_synthetic_")
dir.create(tmp, recursive = TRUE)

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
purity_samples <- sample_universe[1:2, .(sample_barcode)]
aneuploidy_samples <- sample_universe[, .(sample_barcode)]
gene_sets <- env$read_expression_gene_sets(file.path(root, "config", "expression_gene_sets.yaml"))

expression_path <- file.path(tmp, "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")
expression <- data.frame(
  gene_id = c("MKI67|4288", "TOP2A|7153", "PTPRC|5788", "CD8A|925", "VIM|7431", "BRCA1|672", "VEGFA|7422", "KRT8|3856"),
  `TCGA-AB-1234-01A-01R-0000-07` = c(10, 8, 1, 1, 2, 5, 3, 8),
  `TCGA-CD-5678-01A-01R-0000-07` = c(2, 1, 7, 10, 3, 6, 2, 9),
  `TCGA-EF-9999-01A-01R-0000-07` = c(5, 4, 3, 2, 8, 4, 9, 2),
  `TCGA-NO-0001-11A-01R-0000-07` = c(100, 100, 100, 100, 100, 100, 100, 100),
  check.names = FALSE
)
write.table(expression, expression_path, sep = "\t", quote = FALSE, row.names = FALSE)

subset <- env$read_expression_subset(expression_path, gene_sets, sample_universe, chunk_size = 2L)
stopifnot(ncol(subset$expression_matrix) == 3)
stopifnot(nrow(subset$expression_matrix) == 8)
stopifnot("MKI67" %in% rownames(subset$expression_matrix))
stopifnot(!any(subset$sample_metadata$sample_barcode == "TCGA-NO-0001-11A"))

normalized <- env$normalize_expression_matrix(subset$expression_matrix, env$detect_expression_data_type(expression_path))
scores <- env$score_gene_sets(normalized$matrix, subset$sample_metadata, gene_sets)
stopifnot(nrow(scores) == 3)
stopifnot(all(c("proliferation_score", "immune_inflammatory_score", "EMT_score") %in% names(scores)))
stopifnot(any(!is.na(scores$proliferation_score)))
stopifnot(any(!is.na(scores$immune_inflammatory_score)))

pca <- env$compute_expression_pca(normalized$matrix, subset$sample_metadata)
stopifnot(all(paste0("PC", 1:20) %in% names(pca$scores)))
stopifnot(nrow(pca$variance) == 20)

features <- env$build_expression_features(
  projects,
  sample_universe,
  purity_samples,
  aneuploidy_samples,
  gene_sets,
  expression_path,
  "local"
)
stopifnot(nrow(features$sample_scores) == 3)
stopifnot(features$project_scores[project_id == "TCGA-BRCA", n_samples] == 2)
stopifnot(features$project_scores[project_id == "TCGA-LUAD", n_samples] == 1)
stopifnot(features$gene_set_coverage[gene_set == "proliferation", n_found_genes] >= 2)
stopifnot(features$gene_set_coverage[gene_set == "proliferation", n_missing_genes] > 0)
stopifnot(features$pca_status == "selected_gene_expression_pca")

fallback <- env$missing_expression_outputs(projects, sample_universe, gene_sets, "synthetic_no_source")
stopifnot(nrow(fallback$sample_scores) == 3)
stopifnot(all(fallback$sample_scores$source == "no_source_available"))
stopifnot(all(fallback$project_scores$n_samples == 0))

cat("Synthetic expression-layer checks passed\n")
