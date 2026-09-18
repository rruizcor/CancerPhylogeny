#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), mustWork = TRUE)
env <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "02_download_mc3_maf.R"), envir = env)

env$require_mutation_packages()

tmp <- tempfile("mutation_layer_synthetic_")
dir.create(tmp, recursive = TRUE)
dir.create(file.path(tmp, "data", "processed", "mutations"), recursive = TRUE)
dir.create(file.path(tmp, "data", "processed", "features"), recursive = TRUE)
dir.create(file.path(tmp, "results", "tables"), recursive = TRUE)

projects <- data.frame(
  project_id = c("TCGA-BRCA", "TCGA-LUAD"),
  project_code = c("BRCA", "LUAD"),
  stringsAsFactors = FALSE
)
drivers <- data.frame(
  gene = c("TP53", "KRAS", "BRAF", "IDH1", "IDH2", "PIK3CA", "PTEN", "APC", "NF1", "BAP1", "RB1"),
  feature_group = c(
    "tumor_suppressor",
    "oncogene",
    "oncogene",
    "oncometabolite",
    "oncometabolite",
    "oncogene",
    "tumor_suppressor",
    "tumor_suppressor",
    "tumor_suppressor",
    "tumor_suppressor",
    "tumor_suppressor"
  ),
  rationale = rep("test", 11),
  stringsAsFactors = FALSE
)
maf <- data.frame(
  Tumor_Sample_Barcode = c(
    "TCGA-AB-1234-01A-01D-0000-01",
    "TCGA-AB-1234-01A-01D-0000-01",
    "TCGA-CD-5678-01A-01D-0000-01"
  ),
  project_code = c("BRCA", "BRCA", "LUAD"),
  Hugo_Symbol = c("TP53", "SILENT1", "KRAS"),
  Chromosome = c("17", "1", "12"),
  Start_Position = c(7673803, 100, 25398284),
  End_Position = c(7673803, 100, 25398284),
  Reference_Allele = c("C", "A", "G"),
  Tumor_Seq_Allele2 = c("T", "G", "A"),
  Variant_Classification = c("Missense_Mutation", "Silent", "Nonsense_Mutation"),
  Variant_Type = c("SNP", "SNP", "SNP"),
  t_ref_count = c(30, 20, 25),
  t_alt_count = c(12, 5, 15),
  HGVSp_Short = c("p.R175H", "", "p.G12*"),
  stringsAsFactors = FALSE
)

project_path <- file.path(tmp, "tcga_projects.tsv")
driver_path <- file.path(tmp, "driver_genes.csv")
maf_path <- file.path(tmp, "mc3_synthetic.maf")
mutation_output <- file.path(tmp, "data", "processed", "mutations", "mc3_somatic_mutations.parquet")
tmb_output <- file.path(tmp, "data", "processed", "features", "tmb_by_sample.tsv")
prevalence_output <- file.path(tmp, "data", "processed", "features", "mutation_prevalence_by_project.tsv")
driver_prevalence_output <- file.path(tmp, "data", "processed", "features", "mutation_prevalence_by_project_driver_only.tsv")
qc_output <- file.path(tmp, "results", "tables", "mutation_layer_qc_summary.tsv")

write.table(projects, project_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(drivers, driver_path, quote = FALSE, row.names = FALSE)
write.table(maf, maf_path, sep = "\t", quote = FALSE, row.names = FALSE)

standardized <- env$standardize_mc3_maf(maf_path, driver_path, project_path, mutation_output)
parsed <- env$validate_barcode_parsing(standardized$mutations)
tmb <- env$write_tmb_by_sample(standardized$mutations, tmb_output)
invisible(env$write_mutation_prevalence(
  standardized$mutations,
  standardized$expected_projects,
  prevalence_output,
  driver_prevalence_output,
  standardized$driver_genes
))
qc <- env$write_qc_summary(standardized$mutations, tmb, standardized$expected_projects, standardized$driver_genes, qc_output)
env$validate_outputs_nonempty(c(mutation_output, tmb_output, prevalence_output, driver_prevalence_output, qc_output))

stopifnot(nrow(standardized$mutations) == 3)
stopifnot(all(c("TCGA-BRCA", "TCGA-LUAD") %in% standardized$mutations$project_id))
stopifnot(all(parsed$patient_barcode == substr(parsed$Tumor_Sample_Barcode, 1, 12)))
stopifnot(qc$represented_project_count == 2)
stopifnot(qc$sample_count == 2)
stopifnot(isTRUE(qc$ref_alt_counts_present))

tmb_check <- read.delim(tmb_output, stringsAsFactors = FALSE)
stopifnot(all(tmb_check$metric_label == "mutation_count_proxy_no_callable_territory"))

cat("Synthetic mutation-layer checks passed\n")
