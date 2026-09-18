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

nonsynonymous_classes <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Splice_Site",
  "Translation_Start_Site"
)

required_maf_columns <- c(
  "Tumor_Sample_Barcode",
  "Hugo_Symbol",
  "Chromosome",
  "Start_Position",
  "End_Position",
  "Reference_Allele",
  "Tumor_Seq_Allele2",
  "Variant_Classification",
  "Variant_Type"
)

optional_maf_columns <- c("t_ref_count", "t_alt_count", "HGVSp_Short")

expected_driver_columns <- c("gene", "feature_group", "rationale")
required_driver_genes <- c("TP53", "KRAS", "BRAF", "IDH1", "IDH2", "PIK3CA", "PTEN", "APC", "NF1", "BAP1", "RB1")

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  list(refresh = "--refresh" %in% args)
}

require_mutation_packages <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    abort("Package 'data.table' is required for MC3 MAF processing")
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    log_warn("R package 'arrow' is not installed; will try Python pyarrow fallback for Parquet output")
  }
  suppressPackageStartupMessages(library(data.table))
}

write_parquet <- function(data, output_path) {
  ensure_dir(dirname(output_path))

  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(data, output_path)
    log_info(sprintf("Wrote %s", output_path))
    return(invisible(output_path))
  }

  tmp <- tempfile(pattern = "mutations_", fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  data.table::fwrite(data, tmp, sep = "\t", quote = FALSE, na = "")

  py_script <- tempfile(pattern = "write_parquet_", fileext = ".py")
  on.exit(unlink(py_script), add = TRUE)
  writeLines(
    c(
      "import sys",
      "import pyarrow.csv as pv",
      "import pyarrow.parquet as pq",
      "table = pv.read_csv(sys.argv[1], parse_options=pv.ParseOptions(delimiter='\\t'))",
      "pq.write_table(table, sys.argv[2])"
    ),
    py_script
  )

  result <- system2("python", c(py_script, tmp, output_path), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    abort(sprintf("Failed to write Parquet with Python pyarrow fallback: %s", paste(result, collapse = "\n")))
  }

  log_info(sprintf("Wrote %s with Python pyarrow fallback", output_path))
  invisible(output_path)
}

read_projects <- function(path) {
  require_file(path, "TCGA project metadata")
  projects <- data.table::fread(path, sep = "\t", data.table = FALSE)
  required <- c("project_id", "project_code")
  missing <- setdiff(required, names(projects))
  if (length(missing) > 0) {
    abort(sprintf("Project metadata is missing required columns: %s", paste(missing, collapse = ", ")))
  }
  projects
}

read_driver_genes <- function(path) {
  require_file(path, "driver gene configuration")
  warnings <- character()
  drivers <- withCallingHandlers(
    data.table::fread(
      path,
      sep = ",",
      header = TRUE,
      quote = "\"",
      fill = FALSE,
      data.table = FALSE,
      showProgress = FALSE
    ),
    warning = function(warn) {
      warnings <<- c(warnings, conditionMessage(warn))
      invokeRestart("muffleWarning")
    }
  )

  if (length(warnings) > 0) {
    abort(sprintf("Driver gene configuration emitted CSV parsing warnings: %s", paste(unique(warnings), collapse = " | ")))
  }

  if (!identical(names(drivers), expected_driver_columns)) {
    abort(sprintf(
      "Driver gene configuration must have exactly these columns: %s. Found: %s",
      paste(expected_driver_columns, collapse = ", "),
      paste(names(drivers), collapse = ", ")
    ))
  }

  malformed <- which(is.na(drivers$gene) | drivers$gene == "" | is.na(drivers$feature_group) | drivers$feature_group == "" | is.na(drivers$rationale) | drivers$rationale == "")
  if (length(malformed) > 0) {
    abort(sprintf("Driver gene configuration has malformed or blank rows: %s", paste(malformed + 1L, collapse = ", ")))
  }

  drivers$gene <- toupper(trimws(drivers$gene))
  drivers$feature_group <- trimws(drivers$feature_group)
  drivers$rationale <- trimws(drivers$rationale)

  duplicated_genes <- unique(drivers$gene[duplicated(drivers$gene)])
  if (length(duplicated_genes) > 0) {
    abort(sprintf("Driver gene configuration has duplicate gene symbols: %s", paste(duplicated_genes, collapse = ", ")))
  }

  missing_required <- setdiff(required_driver_genes, drivers$gene)
  if (length(missing_required) > 0) {
    abort(sprintf("Driver gene configuration is missing required driver genes: %s", paste(missing_required, collapse = ", ")))
  }

  drivers$gene
}

write_table_gz <- function(data, path) {
  ensure_dir(dirname(path))
  data <- data.table::as.data.table(data)
  tryCatch(
    {
      data.table::fwrite(data, path, sep = "\t", quote = FALSE, na = "", compress = "auto")
    },
    error = function(err) {
      log_warn(sprintf("data.table::fwrite gzip path failed, using gzfile fallback: %s", conditionMessage(err)))
      con <- gzfile(path, open = "wt")
      on.exit(close(con), add = TRUE)
      utils::write.table(data, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    }
  )
  log_info(sprintf("Wrote %s", path))
  invisible(path)
}

read_maf_file <- function(path, select = NULL) {
  require_file(path, "MC3 MAF file")
  log_info(sprintf("Reading MAF: %s", path))
  data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    na.strings = c("", "NA", "."),
    select = select,
    showProgress = interactive()
  )
}

find_local_mc3_maf <- function(raw_dir) {
  if (!dir.exists(raw_dir)) {
    return(NA_character_)
  }

  candidates <- list.files(
    raw_dir,
    pattern = "\\.(maf|maf\\.gz|tsv|tsv\\.gz|txt|txt\\.gz)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  candidates <- candidates[basename(candidates) != ".gitkeep"]
  if (length(candidates) == 0) {
    return(NA_character_)
  }

  preferred <- candidates[grepl("mc3", basename(candidates), ignore.case = TRUE)]
  if (length(preferred) > 0) {
    return(preferred[[1]])
  }
  candidates[[1]]
}

coerce_getmc3maf_result <- function(result, project_code = NA_character_) {
  if (is.character(result) && length(result) == 1 && file.exists(result)) {
    maf <- read_maf_file(result)
  } else if (inherits(result, "data.frame")) {
    maf <- data.table::as.data.table(result)
  } else {
    abort(sprintf("Unsupported getMC3MAF() return type: %s", paste(class(result), collapse = ", ")))
  }

  if (!is.na(project_code) && !"project_code" %in% names(maf)) {
    maf[, project_code := project_code]
  }
  maf
}

download_with_tcgabiolinks <- function(project_codes, output_path) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    log_warn("TCGAbiolinks is not installed; skipping getMC3MAF() retrieval")
    return(FALSE)
  }

  get_mc3_maf <- getExportedValue("TCGAbiolinks", "getMC3MAF")
  get_mc3_args <- names(formals(get_mc3_maf))

  log_info("Attempting MC3 MAF retrieval with TCGAbiolinks::getMC3MAF()")

  maf <- tryCatch(
    {
      if ("tumor" %in% get_mc3_args) {
        pieces <- vector("list", length(project_codes))
        names(pieces) <- project_codes
        for (code in project_codes) {
          log_info(sprintf("Retrieving MC3 MAF records for TCGA-%s", code))
          pieces[[code]] <- coerce_getmc3maf_result(get_mc3_maf(tumor = code), project_code = code)
        }
        data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE)
      } else {
        coerce_getmc3maf_result(get_mc3_maf())
      }
    },
    error = function(err) {
      log_warn(sprintf("TCGAbiolinks::getMC3MAF() failed: %s", conditionMessage(err)))
      NULL
    }
  )

  if (is.null(maf) || nrow(maf) == 0) {
    return(FALSE)
  }

  write_table_gz(maf, output_path)
  TRUE
}

locate_or_download_mc3_maf <- function(raw_dir, project_codes, refresh = FALSE) {
  ensure_dir(raw_dir)
  cached_path <- file.path(raw_dir, "mc3.maf.gz")

  if (file.exists(cached_path) && !refresh) {
    log_info(sprintf("Using cached MC3 MAF: %s", cached_path))
    return(cached_path)
  }

  downloaded <- download_with_tcgabiolinks(project_codes, cached_path)
  if (downloaded) {
    return(cached_path)
  }

  local_path <- find_local_mc3_maf(raw_dir)
  if (!is.na(local_path)) {
    log_warn(sprintf("Using local MC3 MAF fallback: %s", local_path))
    return(local_path)
  }

  abort(paste(
    "Could not retrieve MC3 MAF with TCGAbiolinks and no local MAF was found in data/raw/mc3.",
    "Install TCGAbiolinks in the project environment or place a .maf/.maf.gz file under data/raw/mc3/."
  ))
}

first_present_column <- function(columns, candidates) {
  found <- candidates[candidates %in% columns]
  if (length(found) == 0) {
    return(NA_character_)
  }
  found[[1]]
}

normalize_project_annotations <- function(maf, expected_projects) {
  expected_codes <- expected_projects$project_code
  expected_ids <- expected_projects$project_id

  project_id_col <- first_present_column(names(maf), c("project_id", "Project_ID", "GDC_Project", "gdc_project", "project"))
  project_code_col <- first_present_column(names(maf), c(
    "project_code", "Project_Code", "cancer", "Cancer_Type", "cancer_type",
    "Tumor_Type", "tumor_type", "cohort", "COHORT", "Study", "study"
  ))

  if (!is.na(project_id_col)) {
    values <- as.character(maf[[project_id_col]])
    if (all(grepl("^TCGA-[A-Z0-9]+$", values[!is.na(values) & values != ""]))) {
      maf[, project_id := values]
      maf[, project_code := sub("^TCGA-", "", project_id)]
    }
  }

  if (!is.na(project_code_col) && (!"project_id" %in% names(maf) || !"project_code" %in% names(maf))) {
    values <- as.character(maf[[project_code_col]])
    values <- sub("^TCGA-", "", values)
    values <- toupper(values)
    values[!values %in% expected_codes] <- NA_character_
    maf[, project_code := values]
    maf[, project_id := ifelse(is.na(project_code), NA_character_, paste0("TCGA-", project_code))]
  }

  if (!"project_code" %in% names(maf) || all(is.na(maf$project_code))) {
    abort(paste(
      "Could not identify TCGA project codes in the MC3 MAF.",
      "TCGA barcodes do not encode the cancer project code, so a local combined MAF must include a project/cohort column",
      "or TCGAbiolinks per-project retrieval must be available."
    ))
  }

  maf[!project_id %in% expected_ids, project_id := NA_character_]
  maf[!project_code %in% expected_codes, project_code := NA_character_]
  maf
}

add_missing_optional_columns <- function(maf) {
  for (column in optional_maf_columns) {
    if (!column %in% names(maf)) {
      maf[, (column) := NA]
    }
  }
  maf
}

standardize_mc3_maf <- function(maf_path, driver_gene_path, project_path, output_path) {
  expected_projects <- read_projects(project_path)
  driver_genes <- read_driver_genes(driver_gene_path)

  header <- names(data.table::fread(maf_path, sep = "\t", nrows = 0, quote = ""))
  missing_required <- setdiff(required_maf_columns, header)
  if (length(missing_required) > 0) {
    abort(sprintf("MC3 MAF is missing required columns: %s", paste(missing_required, collapse = ", ")))
  }

  annotation_candidates <- unique(c(
    "project_id", "Project_ID", "GDC_Project", "gdc_project", "project",
    "project_code", "Project_Code", "cancer", "Cancer_Type", "cancer_type",
    "Tumor_Type", "tumor_type", "cohort", "COHORT", "Study", "study"
  ))
  select <- unique(c(required_maf_columns, optional_maf_columns, annotation_candidates))
  select <- select[select %in% header]

  maf <- read_maf_file(maf_path, select = select)
  maf <- add_missing_optional_columns(maf)
  maf <- normalize_project_annotations(maf, expected_projects)

  maf[, Tumor_Sample_Barcode := as.character(Tumor_Sample_Barcode)]
  maf[, sample_barcode := substr(Tumor_Sample_Barcode, 1, 16)]
  maf[, patient_barcode := substr(Tumor_Sample_Barcode, 1, 12)]
  maf[, is_nonsynonymous := Variant_Classification %in% nonsynonymous_classes]
  maf[, is_driver_gene := Hugo_Symbol %in% driver_genes]

  keep <- c(
    "Tumor_Sample_Barcode",
    "project_id",
    "project_code",
    "sample_barcode",
    "patient_barcode",
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Variant_Classification",
    "Variant_Type",
    "t_ref_count",
    "t_alt_count",
    "HGVSp_Short",
    "is_nonsynonymous",
    "is_driver_gene"
  )
  maf <- maf[, ..keep]
  maf <- maf[!is.na(project_id) & !is.na(project_code)]

  write_parquet(maf, output_path)

  list(mutations = maf, driver_genes = driver_genes, expected_projects = expected_projects)
}

write_tmb_by_sample <- function(mutations, output_path) {
  ensure_dir(dirname(output_path))
  tmb <- mutations[, .(
    nonsynonymous_count = sum(is_nonsynonymous, na.rm = TRUE),
    total_mutation_count = .N
  ), by = .(sample_barcode, patient_barcode, project_id, project_code)]
  tmb[, metric_label := "mutation_count_proxy_no_callable_territory"]
  data.table::setorder(tmb, project_id, sample_barcode)
  data.table::fwrite(tmb, output_path, sep = "\t", quote = FALSE, na = "NA")
  log_info(sprintf("Wrote %s", output_path))
  tmb
}

write_mutation_prevalence <- function(mutations, expected_projects, output_path, driver_only_output_path, driver_genes) {
  ensure_dir(dirname(output_path))
  sample_counts <- unique(mutations[, .(project_id, project_code, sample_barcode)])[, .(
    n_samples = uniqueN(sample_barcode)
  ), by = .(project_id, project_code)]

  nonsyn <- mutations[is_nonsynonymous == TRUE & !is.na(Hugo_Symbol) & Hugo_Symbol != ""]
  gene_sample <- unique(nonsyn[, .(project_id, project_code, sample_barcode, Hugo_Symbol)])
  prevalence_long <- gene_sample[, .(n_mutated = uniqueN(sample_barcode)), by = .(project_id, project_code, Hugo_Symbol)]
  prevalence_long <- merge(prevalence_long, sample_counts, by = c("project_id", "project_code"), all.x = TRUE)
  prevalence_long[, prevalence := n_mutated / n_samples]

  if (nrow(prevalence_long) == 0) {
    abort("No nonsynonymous mutation records were available for project-level mutation prevalence")
  }

  prevalence_wide <- data.table::dcast(
    prevalence_long,
    project_id + project_code ~ Hugo_Symbol,
    value.var = "prevalence",
    fill = 0
  )
  prevalence_wide <- merge(
    expected_projects[, c("project_id", "project_code")],
    prevalence_wide,
    by = c("project_id", "project_code"),
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setDT(prevalence_wide)
  for (column in setdiff(names(prevalence_wide), c("project_id", "project_code"))) {
    data.table::set(prevalence_wide, which(is.na(prevalence_wide[[column]])), column, 0)
  }
  data.table::setorder(prevalence_wide, project_id)
  data.table::fwrite(prevalence_wide, output_path, sep = "\t", quote = FALSE, na = "NA")
  log_info(sprintf("Wrote %s", output_path))

  driver_long <- prevalence_long[Hugo_Symbol %in% driver_genes]
  if (nrow(driver_long) == 0) {
    log_warn("No configured driver genes had nonsynonymous mutation records in the MC3 MAF")
    driver_wide <- data.table::as.data.table(expected_projects[, c("project_id", "project_code")])
    for (gene in driver_genes) {
      driver_wide[, (gene) := 0]
    }
  } else {
    driver_wide <- data.table::dcast(
      driver_long,
      project_id + project_code ~ Hugo_Symbol,
      value.var = "prevalence",
      fill = 0
    )
    driver_wide <- merge(
      expected_projects[, c("project_id", "project_code")],
      driver_wide,
      by = c("project_id", "project_code"),
      all.x = TRUE,
      sort = FALSE
    )
    data.table::setDT(driver_wide)
    for (gene in setdiff(driver_genes, names(driver_wide))) {
      driver_wide[, (gene) := 0]
    }
    for (column in setdiff(names(driver_wide), c("project_id", "project_code"))) {
      data.table::set(driver_wide, which(is.na(driver_wide[[column]])), column, 0)
    }
    driver_wide <- driver_wide[, c("project_id", "project_code", driver_genes), with = FALSE]
  }

  data.table::setorder(driver_wide, project_id)
  data.table::fwrite(driver_wide, driver_only_output_path, sep = "\t", quote = FALSE, na = "NA")
  log_info(sprintf("Wrote %s", driver_only_output_path))

  prevalence_long
}

validate_barcode_parsing <- function(mutations) {
  examples <- unique(mutations$Tumor_Sample_Barcode)[seq_len(min(10, data.table::uniqueN(mutations$Tumor_Sample_Barcode)))]
  parsed <- data.table::data.table(
    Tumor_Sample_Barcode = examples,
    sample_barcode = substr(examples, 1, 16),
    patient_barcode = substr(examples, 1, 12)
  )
  invalid <- parsed[!grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", Tumor_Sample_Barcode)]
  if (nrow(invalid) > 0) {
    log_warn(sprintf("Barcode parsing examples include non-standard TCGA barcodes: %s", paste(invalid$Tumor_Sample_Barcode, collapse = ", ")))
  } else {
    log_info(sprintf("Barcode parsing examples passed: %s", paste(parsed$sample_barcode, collapse = ", ")))
  }
  parsed
}

write_qc_summary <- function(mutations, tmb, expected_projects, driver_genes, output_path) {
  ensure_dir(dirname(output_path))
  represented <- sort(unique(mutations$project_id))
  missing_projects <- setdiff(expected_projects$project_id, represented)
  absent_driver_genes <- setdiff(driver_genes, unique(mutations$Hugo_Symbol))

  if (length(missing_projects) > 0) {
    log_warn(sprintf("Expected TCGA projects with no mutation records: %s", paste(missing_projects, collapse = ", ")))
  }
  if (length(absent_driver_genes) > 0) {
    log_warn(sprintf("Configured driver genes absent from MAF: %s", paste(absent_driver_genes, collapse = ", ")))
  }

  ref_alt_present <- all(c("t_ref_count", "t_alt_count") %in% names(mutations)) &&
    any(!is.na(mutations$t_ref_count)) &&
    any(!is.na(mutations$t_alt_count))
  hgvsp_present <- "HGVSp_Short" %in% names(mutations) && any(!is.na(mutations$HGVSp_Short) & mutations$HGVSp_Short != "")

  global_rows <- data.table::data.table(
    qc_section = "global",
    project_id = "ALL",
    project_code = "ALL",
    metric = c(
      "total_mutation_records",
      "unique_tumor_samples",
      "unique_patients",
      "tcga_projects_detected",
      "expected_tcga_projects_missing",
      "missing_tcga_projects",
      "ref_alt_counts_present",
      "hgvsp_short_present",
      "configured_driver_genes_absent"
    ),
    value = c(
      as.character(nrow(mutations)),
      as.character(data.table::uniqueN(mutations$sample_barcode)),
      as.character(data.table::uniqueN(mutations$patient_barcode)),
      as.character(length(represented)),
      as.character(length(missing_projects)),
      ifelse(length(missing_projects) == 0, "none", paste(missing_projects, collapse = ";")),
      as.character(ref_alt_present),
      as.character(hgvsp_present),
      ifelse(length(absent_driver_genes) == 0, "none", paste(absent_driver_genes, collapse = ";"))
    )
  )

  sample_counts <- unique(mutations[, .(project_id, project_code, sample_barcode)])[, .(
    value = as.character(uniqueN(sample_barcode))
  ), by = .(project_id, project_code)]
  sample_counts[, `:=`(qc_section = "per_project", metric = "sample_count")]

  nonsyn_summary <- tmb[, .(
    nonsynonymous_mutation_count_min = as.numeric(min(nonsynonymous_count, na.rm = TRUE)),
    nonsynonymous_mutation_count_median = as.numeric(stats::median(nonsynonymous_count, na.rm = TRUE)),
    nonsynonymous_mutation_count_mean = mean(nonsynonymous_count, na.rm = TRUE),
    nonsynonymous_mutation_count_max = as.numeric(max(nonsynonymous_count, na.rm = TRUE))
  ), by = .(project_id, project_code)]
  nonsyn_long <- data.table::melt(
    nonsyn_summary,
    id.vars = c("project_id", "project_code"),
    variable.name = "metric",
    value.name = "value"
  )
  nonsyn_long[, value := as.character(round(as.numeric(value), 4))]
  nonsyn_long[, qc_section := "per_project"]

  qc <- data.table::rbindlist(
    list(
      global_rows,
      sample_counts[, .(qc_section, project_id, project_code, metric, value)],
      nonsyn_long[, .(qc_section, project_id, project_code, metric, value)]
    ),
    use.names = TRUE,
    fill = TRUE
  )
  data.table::setorder(qc, qc_section, project_id, metric)
  data.table::fwrite(qc, output_path, sep = "\t", quote = FALSE, na = "NA")
  log_info(sprintf("Wrote %s", output_path))

  list(
    represented_project_count = length(represented),
    sample_count = data.table::uniqueN(mutations$sample_barcode),
    ref_alt_counts_present = ref_alt_present,
    hgvsp_short_present = hgvsp_present,
    missing_projects = missing_projects,
    absent_driver_genes = absent_driver_genes
  )
}

validate_outputs_nonempty <- function(paths) {
  for (path in paths) {
    if (!file.exists(path) || file.info(path)$size == 0) {
      abort(sprintf("Expected non-empty output was not created: %s", path))
    }
  }
  invisible(TRUE)
}

main <- function() {
  args <- parse_args()
  require_mutation_packages()

  root <- find_project_root()
  raw_dir <- project_path("data", "raw", "mc3", root = root)
  projects_path <- project_path("data", "interim", "tcga_projects.tsv", root = root)
  driver_gene_path <- project_path("config", "driver_genes.csv", root = root)

  mutation_output <- project_path("data", "processed", "mutations", "mc3_somatic_mutations.parquet", root = root)
  tmb_output <- project_path("data", "processed", "features", "tmb_by_sample.tsv", root = root)
  prevalence_output <- project_path("data", "processed", "features", "mutation_prevalence_by_project.tsv", root = root)
  driver_prevalence_output <- project_path("data", "processed", "features", "mutation_prevalence_by_project_driver_only.tsv", root = root)
  qc_output <- project_path("results", "tables", "mutation_layer_qc_summary.tsv", root = root)

  projects <- read_projects(projects_path)
  log_info(sprintf("Preparing mutation layer for %d expected TCGA projects", nrow(projects)))

  maf_path <- locate_or_download_mc3_maf(raw_dir, projects$project_code, refresh = args$refresh)
  standardized <- standardize_mc3_maf(maf_path, driver_gene_path, projects_path, mutation_output)
  validate_barcode_parsing(standardized$mutations)
  tmb <- write_tmb_by_sample(standardized$mutations, tmb_output)
  write_mutation_prevalence(
    standardized$mutations,
    standardized$expected_projects,
    prevalence_output,
    driver_prevalence_output,
    standardized$driver_genes
  )
  qc <- write_qc_summary(standardized$mutations, tmb, standardized$expected_projects, standardized$driver_genes, qc_output)

  validate_outputs_nonempty(c(
    mutation_output,
    tmb_output,
    prevalence_output,
    driver_prevalence_output,
    qc_output
  ))

  log_info(sprintf("Mutation layer complete: %d projects represented, %d samples processed", qc$represented_project_count, qc$sample_count))
  log_info(sprintf("Ref/alt counts available: %s", qc$ref_alt_counts_present))
}

if (identical(environment(), globalenv())) {
  main()
}
