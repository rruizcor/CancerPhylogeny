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

copy_number_numeric_columns <- c(
  "purity",
  "ploidy",
  "aneuploidy_score",
  "fraction_genome_altered",
  "arm_gain_count",
  "arm_loss_count",
  "total_arm_alteration_count",
  "copy_number_complexity_score"
)

pancanatlas_resource_specs <- function(root) {
  data.table::data.table(
    resource_key = c("absolute_purity_ploidy", "arm_calls_aneuploidy"),
    expected_file = c("TCGA_mastercalls.abs_tables_JSedit.fixed.txt", "PANCAN_ArmCallsAndAneuploidyScore_092817.txt"),
    local_path = c(
      project_path("data", "raw", "purity_ploidy", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt", root = root),
      project_path("data", "raw", "copy_number", "PANCAN_ArmCallsAndAneuploidyScore_092817.txt", root = root)
    ),
    gdc_url = c(
      "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
      "https://api.gdc.cancer.gov/data/4c35f34f-b0f3-4891-8794-4840dd748aad"
    ),
    publication_page = c(
      "https://gdc.cancer.gov/about-data/publications/pancanatlas",
      "https://gdc.cancer.gov/about-data/publications/pancan-aneuploidy"
    )
  )
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  list(refresh = "--refresh" %in% args)
}

require_copy_number_packages <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    abort("Package 'data.table' is required for copy-number processing")
  }
  suppressPackageStartupMessages(library(data.table))
}

normalize_name <- function(x) {
  tolower(gsub("[^a-zA-Z0-9]+", "_", x))
}

find_column <- function(table, candidates) {
  normalized <- normalize_name(names(table))
  candidate_norm <- normalize_name(candidates)
  match_index <- match(candidate_norm, normalized, nomatch = 0L)
  match_index <- match_index[match_index > 0L]
  if (length(match_index) == 0) {
    return(NA_character_)
  }
  names(table)[match_index[[1]]]
}

coerce_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

read_projects <- function(path) {
  require_file(path, "TCGA project metadata")
  projects <- data.table::fread(path, sep = "\t", data.table = FALSE)
  required <- c("project_id", "project_code")
  missing <- setdiff(required, names(projects))
  if (length(missing) > 0) {
    abort(sprintf("Project metadata is missing required columns: %s", paste(missing, collapse = ", ")))
  }
  data.table::as.data.table(projects[, required])
}

read_parquet_with_fallback <- function(path) {
  if (requireNamespace("arrow", quietly = TRUE)) {
    return(data.table::as.data.table(arrow::read_parquet(path)))
  }

  tmp <- tempfile(pattern = "parquet_", fileext = ".tsv")
  py_script <- tempfile(pattern = "read_parquet_", fileext = ".py")
  on.exit(unlink(c(tmp, py_script)), add = TRUE)
  writeLines(
    c(
      "import sys",
      "import pyarrow.parquet as pq",
      "import pandas as pd",
      "table = pq.read_table(sys.argv[1])",
      "table.to_pandas().to_csv(sys.argv[2], sep='\\t', index=False)"
    ),
    py_script
  )
  result <- system2("python", c(py_script, path, tmp), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    abort(sprintf("Failed to read Parquet with Python pyarrow fallback: %s", paste(result, collapse = "\n")))
  }
  data.table::fread(tmp, sep = "\t", data.table = TRUE)
}

read_table_auto <- function(path) {
  require_file(path, "input table")
  log_info(sprintf("Reading local copy-number resource: %s", path))
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    return(read_parquet_with_fallback(path))
  }
  sep <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
  data.table::fread(path, sep = sep, quote = "\"", data.table = TRUE, showProgress = interactive())
}

write_tsv_table <- function(data, path) {
  ensure_dir(dirname(path))
  data.table::fwrite(data, path, sep = "\t", quote = FALSE, na = "NA")
  log_info(sprintf("Wrote %s", path))
  invisible(path)
}

write_parquet_or_tsv <- function(data, path) {
  ensure_dir(dirname(path))
  if (grepl("\\.parquet$", path, ignore.case = TRUE) && requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(data, path)
    log_info(sprintf("Wrote %s", path))
    return(invisible(path))
  }
  write_tsv_table(data, path)
}

read_mutation_sample_universe <- function(mutation_path) {
  require_file(mutation_path, "mutation-layer parquet")
  mutations <- read_parquet_with_fallback(mutation_path)
  required <- c("sample_barcode", "patient_barcode", "project_id", "project_code")
  missing <- setdiff(required, names(mutations))
  if (length(missing) > 0) {
    abort(sprintf("Mutation table is missing required sample columns: %s", paste(missing, collapse = ", ")))
  }
  universe <- unique(mutations[, ..required])
  data.table::setorder(universe, project_id, sample_barcode)
  universe
}

download_public_resource <- function(url, destination) {
  ensure_dir(dirname(destination))
  tmp <- tempfile(pattern = "pancanatlas_", fileext = ".download")
  on.exit(unlink(tmp), add = TRUE)
  tryCatch(
    {
      utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
      if (!file.exists(tmp) || file.info(tmp)$size == 0) {
        stop("downloaded file is empty", call. = FALSE)
      }
      file.copy(tmp, destination, overwrite = TRUE)
      TRUE
    },
    error = function(err) {
      log_warn(sprintf("Download failed for %s: %s", url, conditionMessage(err)))
      FALSE
    }
  )
}

manual_download_message <- function(specs) {
  paste(
    "Automatic download failed for one or more PanCanAtlas copy-number resources.",
    "Manual fallback:",
    paste(sprintf("- Download %s from %s and place it at %s", specs$expected_file, specs$publication_page, specs$local_path), collapse = "\n"),
    sep = "\n"
  )
}

retrieve_pancanatlas_resources <- function(root, refresh = FALSE) {
  specs <- pancanatlas_resource_specs(root)
  specs[, `:=`(status = "unavailable", path = local_path, file_size = NA_real_)]

  for (i in seq_len(nrow(specs))) {
    local_path <- specs$local_path[[i]]
    ensure_dir(dirname(local_path))

    if (file.exists(local_path) && file.info(local_path)$size > 0 && !refresh) {
      specs[i, `:=`(status = "local", file_size = file.info(local_path)$size)]
      log_info(sprintf("Using local PanCanAtlas resource: %s", local_path))
      next
    }

    log_info(sprintf("Attempting download of %s from GDC", specs$expected_file[[i]]))
    downloaded <- download_public_resource(specs$gdc_url[[i]], local_path)
    if (downloaded && file.exists(local_path) && file.info(local_path)$size > 0) {
      specs[i, `:=`(status = "downloaded", file_size = file.info(local_path)$size)]
      log_info(sprintf("Downloaded PanCanAtlas resource: %s", local_path))
    }
  }

  unavailable <- specs[status == "unavailable"]
  if (nrow(unavailable) > 0) {
    log_warn(manual_download_message(unavailable))
  }

  specs
}

discover_local_resources <- function(root) {
  raw_dirs <- c(
    project_path("data", "raw", "copy_number", root = root),
    project_path("data", "raw", "purity_ploidy", root = root)
  )
  for (raw_dir in raw_dirs) {
    ensure_dir(raw_dir)
  }
  patterns <- "\\.(tsv|tsv\\.gz|csv|csv\\.gz|txt|txt\\.gz|maf|maf\\.gz|parquet)$"
  files <- unlist(lapply(raw_dirs, function(raw_dir) {
    list.files(raw_dir, pattern = patterns, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  files <- files[basename(files) != ".gitkeep"]
  files
}

classify_resource <- function(path) {
  name <- tolower(basename(path))
  if (grepl("purity|ploidy|absolute|ascat", name)) {
    return("purity_ploidy")
  }
  if (grepl("arm|aneuploid", name)) {
    return("arm_or_aneuploidy")
  }
  if (grepl("segment|seg\\.|cna|copy", name)) {
    return("segments")
  }
  "unknown"
}

harmonize_sample_keys <- function(table, sample_universe) {
  data.table::setDT(table)
  sample_col <- find_column(table, c(
    "sample_barcode", "tumor_sample_barcode", "Tumor_Sample_Barcode", "sample",
    "sample_id", "aliquot_barcode", "barcode", "bcr_sample_barcode"
  ))
  patient_col <- find_column(table, c("patient_barcode", "participant_barcode", "case_id", "patient", "bcr_patient_barcode"))
  project_id_col <- find_column(table, c("project_id", "Project_ID", "gdc_project"))
  project_code_col <- find_column(table, c("project_code", "cancer_type", "tumor_type", "type", "cohort", "study"))

  if (!is.na(sample_col)) {
    table[, sample_barcode := substr(as.character(get(sample_col)), 1, 16)]
    table[, patient_barcode := substr(as.character(get(sample_col)), 1, 12)]
  } else if (!is.na(patient_col)) {
    table[, patient_barcode := substr(as.character(get(patient_col)), 1, 12)]
    table[, sample_barcode := NA_character_]
  } else {
    table[, `:=`(sample_barcode = NA_character_, patient_barcode = NA_character_)]
  }

  if (!is.na(project_id_col)) {
    table[, project_id := as.character(get(project_id_col))]
    table[, project_code := sub("^TCGA-", "", project_id)]
  } else if (!is.na(project_code_col)) {
    table[, project_code := toupper(sub("^TCGA-", "", as.character(get(project_code_col))))]
    table[, project_id := ifelse(is.na(project_code) | project_code == "", NA_character_, paste0("TCGA-", project_code))]
  } else {
    table[, `:=`(project_id = NA_character_, project_code = NA_character_)]
  }

  by_sample <- sample_universe[, .(sample_barcode, universe_patient_barcode = patient_barcode, universe_project_id = project_id, universe_project_code = project_code)]
  table <- merge(table, by_sample, by = "sample_barcode", all.x = TRUE, sort = FALSE)
  table[is.na(patient_barcode) | patient_barcode == "", patient_barcode := universe_patient_barcode]
  table[is.na(project_id) | project_id == "", project_id := universe_project_id]
  table[is.na(project_code) | project_code == "", project_code := universe_project_code]
  table[, c("universe_patient_barcode", "universe_project_id", "universe_project_code") := NULL]

  if (any(is.na(table$sample_barcode)) && !all(is.na(table$patient_barcode))) {
    patient_map <- sample_universe[, .N, by = patient_barcode][N == 1][sample_universe, on = "patient_barcode", nomatch = 0]
    patient_map <- patient_map[, .(patient_barcode, mapped_sample_barcode = sample_barcode, mapped_project_id = project_id, mapped_project_code = project_code)]
    table <- merge(table, patient_map, by = "patient_barcode", all.x = TRUE, sort = FALSE)
    table[is.na(sample_barcode) | sample_barcode == "", sample_barcode := mapped_sample_barcode]
    table[is.na(project_id) | project_id == "", project_id := mapped_project_id]
    table[is.na(project_code) | project_code == "", project_code := mapped_project_code]
    table[, c("mapped_sample_barcode", "mapped_project_id", "mapped_project_code") := NULL]
  }

  table
}

normalize_project_code <- function(project_code) {
  project_code <- toupper(sub("^TCGA-", "", as.character(project_code)))
  project_code[project_code == "OVCA"] <- "OV"
  project_code
}

align_to_sample_universe <- function(table, sample_universe) {
  data.table::setDT(table)
  table[, input_sample_barcode := sample_barcode]

  exact_map <- sample_universe[, .(
    input_sample_barcode = sample_barcode,
    mapped_sample_barcode = sample_barcode,
    mapped_patient_barcode = patient_barcode,
    mapped_project_id = project_id,
    mapped_project_code = project_code
  )]
  table <- merge(table, exact_map, by = "input_sample_barcode", all.x = TRUE, sort = FALSE)
  table[!is.na(mapped_sample_barcode), sample_barcode := mapped_sample_barcode]
  table[is.na(patient_barcode) | patient_barcode == "", patient_barcode := mapped_patient_barcode]
  table[is.na(project_id) | project_id == "", project_id := mapped_project_id]
  table[is.na(project_code) | project_code == "", project_code := mapped_project_code]
  table[, c("mapped_sample_barcode", "mapped_patient_barcode", "mapped_project_id", "mapped_project_code") := NULL]

  prefix_universe <- copy(sample_universe)
  prefix_universe[, sample_prefix15 := substr(sample_barcode, 1, 15)]
  unique_prefixes <- prefix_universe[, .N, by = sample_prefix15][N == 1]
  prefix_map <- prefix_universe[unique_prefixes, on = "sample_prefix15"][
    ,
    .(
      sample_prefix15,
      prefix_sample_barcode = sample_barcode,
      prefix_patient_barcode = patient_barcode,
      prefix_project_id = project_id,
      prefix_project_code = project_code
    )
  ]

  table[, sample_prefix15 := substr(input_sample_barcode, 1, 15)]
  table <- merge(table, prefix_map, by = "sample_prefix15", all.x = TRUE, sort = FALSE)
  table[(is.na(project_id) | project_id == "") & !is.na(prefix_project_id), project_id := prefix_project_id]
  table[(is.na(project_code) | project_code == "") & !is.na(prefix_project_code), project_code := prefix_project_code]
  table[(is.na(patient_barcode) | patient_barcode == "") & !is.na(prefix_patient_barcode), patient_barcode := prefix_patient_barcode]
  table[nchar(sample_barcode) < 16 & !is.na(prefix_sample_barcode), sample_barcode := prefix_sample_barcode]
  table[, c("input_sample_barcode", "sample_prefix15", "prefix_sample_barcode", "prefix_patient_barcode", "prefix_project_id", "prefix_project_code") := NULL]
  table
}

extract_absolute_purity_ploidy <- function(table, source_name, sample_universe) {
  sample_col <- find_column(table, c("sample", "Tumor_Sample_Barcode", "sample_barcode", "array"))
  purity_col <- find_column(table, c("purity"))
  ploidy_col <- find_column(table, c("ploidy"))
  genome_doubling_col <- find_column(table, c("Genome doublings", "genome_doublings", "genome doublings"))
  call_status_col <- find_column(table, c("call status", "call_status"))
  subclonal_col <- find_column(table, c("Subclonal genome fraction", "subclonal_genome_fraction"))

  if (is.na(sample_col) || (is.na(purity_col) && is.na(ploidy_col))) {
    return(NULL)
  }

  out <- data.table(
    sample_barcode = substr(as.character(table[[sample_col]]), 1, 16),
    patient_barcode = substr(as.character(table[[sample_col]]), 1, 12),
    purity = if (is.na(purity_col)) NA_real_ else coerce_numeric(table[[purity_col]]),
    ploidy = if (is.na(ploidy_col)) NA_real_ else coerce_numeric(table[[ploidy_col]]),
    genome_doublings = if (is.na(genome_doubling_col)) NA_real_ else coerce_numeric(table[[genome_doubling_col]]),
    subclonal_genome_fraction = if (is.na(subclonal_col)) NA_real_ else coerce_numeric(table[[subclonal_col]])
  )
  out[, whole_genome_doubling_status := fifelse(
    is.na(genome_doublings),
    NA_character_,
    fifelse(genome_doublings >= 1, "WGD", "no_WGD")
  )]
  out[, source := source_name]
  out[, notes := if (is.na(call_status_col)) "absolute_purity_ploidy" else paste0("absolute_call_status=", as.character(table[[call_status_col]]))]
  out[, `:=`(project_id = NA_character_, project_code = NA_character_)]
  out <- align_to_sample_universe(out, sample_universe)
  out[!is.na(project_id)]
}

arm_name_columns <- function(table) {
  excluded <- normalize_name(c(
    "sample", "type", "aneuploidy score", "aneuploidy_score", "fraction genome altered",
    "fraction_genome_altered", "fga", "sample_barcode", "project_code", "project_id"
  ))
  names(table)[!normalize_name(names(table)) %in% excluded & grepl("(^[0-9XY]+[pq]$)|(^[0-9]+ \\([0-9]+[pq]\\)$)", names(table))]
}

clean_arm_name <- function(x) {
  out <- as.character(x)
  parenthetical <- grepl("\\(", out)
  out[parenthetical] <- sub("^.*\\(([0-9XY]+[pq])\\).*$", "\\1", out[parenthetical])
  out
}

extract_wide_arm_calls <- function(table, source_name, sample_universe) {
  sample_col <- find_column(table, c("Sample", "sample", "sample_barcode", "Tumor_Sample_Barcode"))
  project_code_col <- find_column(table, c("Type", "type", "project_code", "cancer_type"))
  if (is.na(sample_col)) {
    return(NULL)
  }
  arm_cols <- arm_name_columns(table)
  if (length(arm_cols) == 0) {
    return(NULL)
  }

  keyed <- data.table(
    sample_barcode = substr(as.character(table[[sample_col]]), 1, 16),
    patient_barcode = substr(as.character(table[[sample_col]]), 1, 12),
    project_code = if (is.na(project_code_col)) NA_character_ else normalize_project_code(table[[project_code_col]])
  )
  keyed[, project_id := ifelse(is.na(project_code) | project_code == "", NA_character_, paste0("TCGA-", project_code))]
  keyed <- align_to_sample_universe(keyed, sample_universe)

  values <- copy(table[, ..arm_cols])
  values[, sample_barcode := keyed$sample_barcode]
  calls <- data.table::melt(
    values,
    id.vars = "sample_barcode",
    variable.name = "chromosome_arm",
    value.name = "arm_call"
  )
  calls <- merge(calls, keyed, by = "sample_barcode", all.x = TRUE, sort = FALSE)
  calls[, chromosome_arm := clean_arm_name(chromosome_arm)]
  calls[, arm_call := as.character(arm_call)]
  calls[, source := source_name]
  calls[, call_normalized := trimws(tolower(arm_call))]
  calls[, is_gain := call_normalized %in% c("1", "+1", "gain", "amp", "amplification")]
  calls[, is_loss := call_normalized %in% c("-1", "loss", "del", "deletion")]
  calls[!is.na(sample_barcode) & !is.na(arm_call) & arm_call != ""]
}

extract_pancan_aneuploidy_scores <- function(table, source_name, sample_universe) {
  sample_col <- find_column(table, c("Sample", "sample", "sample_barcode", "Tumor_Sample_Barcode"))
  project_code_col <- find_column(table, c("Type", "type", "project_code", "cancer_type"))
  aneuploidy_col <- find_column(table, c("Aneuploidy Score", "aneuploidy_score"))
  fga_col <- find_column(table, c("fraction_genome_altered", "fraction genome altered", "fga"))
  if (is.na(sample_col) || is.na(aneuploidy_col)) {
    return(NULL)
  }

  out <- data.table(
    sample_barcode = substr(as.character(table[[sample_col]]), 1, 16),
    patient_barcode = substr(as.character(table[[sample_col]]), 1, 12),
    project_code = if (is.na(project_code_col)) NA_character_ else normalize_project_code(table[[project_code_col]]),
    aneuploidy_score = coerce_numeric(table[[aneuploidy_col]]),
    fraction_genome_altered = if (is.na(fga_col)) NA_real_ else coerce_numeric(table[[fga_col]])
  )
  out[, project_id := ifelse(is.na(project_code) | project_code == "", NA_character_, paste0("TCGA-", project_code))]
  out <- align_to_sample_universe(out, sample_universe)
  out[, `:=`(
    arm_gain_count = NA_real_,
    arm_loss_count = NA_real_,
    total_arm_alteration_count = NA_real_,
    copy_number_complexity_score = NA_real_,
    source = source_name,
    notes = "pancanatlas_aneuploidy_score"
  )]
  out[!is.na(project_id)]
}

extract_purity_ploidy <- function(table, source_name, sample_universe) {
  table <- harmonize_sample_keys(table, sample_universe)
  purity_col <- find_column(table, c("purity", "tumor_purity", "absolute_purity", "purity_estimate", "cellularity"))
  ploidy_col <- find_column(table, c("ploidy", "absolute_ploidy", "tumor_ploidy", "ploidy_estimate"))
  genome_doubling_col <- find_column(table, c("Genome doublings", "genome_doublings", "whole_genome_doubling_status", "wgd"))
  subclonal_col <- find_column(table, c("Subclonal genome fraction", "subclonal_genome_fraction"))

  if (is.na(purity_col) && is.na(ploidy_col)) {
    return(NULL)
  }

  out <- table[, .(sample_barcode, patient_barcode, project_id, project_code)]
  out[, purity := if (is.na(purity_col)) NA_real_ else coerce_numeric(table[[purity_col]])]
  out[, ploidy := if (is.na(ploidy_col)) NA_real_ else coerce_numeric(table[[ploidy_col]])]
  out[, genome_doublings := if (is.na(genome_doubling_col)) NA_real_ else coerce_numeric(table[[genome_doubling_col]])]
  out[, whole_genome_doubling_status := fifelse(
    is.na(genome_doublings),
    NA_character_,
    fifelse(genome_doublings >= 1, "WGD", "no_WGD")
  )]
  out[, subclonal_genome_fraction := if (is.na(subclonal_col)) NA_real_ else coerce_numeric(table[[subclonal_col]])]
  out[, source := source_name]
  out[, notes := fifelse(is.na(purity) & is.na(ploidy), "purity_ploidy_missing", "direct_or_local_source")]
  out[!is.na(sample_barcode)]
}

extract_direct_aneuploidy <- function(table, source_name, sample_universe) {
  table <- harmonize_sample_keys(table, sample_universe)
  metric_map <- list(
    aneuploidy_score = c("aneuploidy_score", "aneuploidy", "aneuploidy_scores", "number_of_aneuploid_events"),
    fraction_genome_altered = c("fraction_genome_altered", "fga", "genome_instability", "fraction_altered"),
    arm_gain_count = c("arm_gain_count", "aneuploidy_gain", "gain_count"),
    arm_loss_count = c("arm_loss_count", "aneuploidy_loss", "loss_count"),
    total_arm_alteration_count = c("total_arm_alteration_count", "arm_alteration_count", "aneuploidy_count"),
    copy_number_complexity_score = c("copy_number_complexity_score", "cn_complexity", "number_of_segments", "segment_count")
  )
  found <- vapply(metric_map, function(candidates) find_column(table, candidates), character(1))
  if (all(is.na(found))) {
    return(NULL)
  }

  out <- table[, .(sample_barcode, patient_barcode, project_id, project_code)]
  for (metric in names(metric_map)) {
    column <- found[[metric]]
    out[, (metric) := if (is.na(column)) NA_real_ else coerce_numeric(table[[column]])]
  }
  out[is.na(total_arm_alteration_count) & (!is.na(arm_gain_count) | !is.na(arm_loss_count)),
      total_arm_alteration_count := fifelse(is.na(arm_gain_count), 0, arm_gain_count) + fifelse(is.na(arm_loss_count), 0, arm_loss_count)]
  out[, source := source_name]
  out[, notes := "direct_or_local_source"]
  out[!is.na(sample_barcode)]
}

extract_arm_calls <- function(table, source_name, sample_universe) {
  table <- harmonize_sample_keys(table, sample_universe)
  arm_col <- find_column(table, c("arm", "chromosome_arm", "chrom_arm", "cytoband_arm"))
  call_col <- find_column(table, c("call", "arm_call", "alteration", "status", "cn_call", "copy_number_status"))

  if (is.na(arm_col) || is.na(call_col)) {
    return(NULL)
  }

  calls <- table[, .(
    sample_barcode,
    patient_barcode,
    project_id,
    project_code,
    chromosome_arm = as.character(get(arm_col)),
    arm_call = as.character(get(call_col))
  )]
  calls[, source := source_name]
  calls[, call_normalized := tolower(trimws(arm_call))]
  calls[, is_gain := grepl("gain|amp|\\+|^1$", call_normalized)]
  calls[, is_loss := grepl("loss|del|-|^-1$", call_normalized)]
  calls[!is.na(sample_barcode)]
}

summarize_arm_calls <- function(arm_calls) {
  if (is.null(arm_calls) || nrow(arm_calls) == 0) {
    return(NULL)
  }
  arm_calls[, .(
    arm_gain_count = sum(is_gain, na.rm = TRUE),
    arm_loss_count = sum(is_loss, na.rm = TRUE),
    total_arm_alteration_count = sum(is_gain | is_loss, na.rm = TRUE)
  ), by = .(sample_barcode, patient_barcode, project_id, project_code)]
}

extract_segments <- function(table, source_name, sample_universe) {
  table <- harmonize_sample_keys(table, sample_universe)
  chr_col <- find_column(table, c("chromosome", "chrom", "chr"))
  start_col <- find_column(table, c("start", "start_position", "loc_start"))
  end_col <- find_column(table, c("end", "end_position", "loc_end"))
  value_col <- find_column(table, c("segment_mean", "seg_mean", "log2", "log2_copy_ratio", "copy_number", "cn"))

  if (is.na(chr_col) || is.na(start_col) || is.na(end_col)) {
    return(NULL)
  }

  segments <- table[, .(
    sample_barcode,
    patient_barcode,
    project_id,
    project_code,
    chromosome = as.character(get(chr_col)),
    start = coerce_numeric(get(start_col)),
    end = coerce_numeric(get(end_col))
  )]
  segments[, segment_value := if (is.na(value_col)) NA_real_ else coerce_numeric(table[[value_col]])]
  segments[, source := source_name]
  segments[!is.na(sample_barcode)]
}

summarize_segments <- function(segments) {
  if (is.null(segments) || nrow(segments) == 0) {
    return(NULL)
  }
  segments[, segment_length := pmax(0, end - start + 1)]
  segments[, altered_segment := ifelse(is.na(segment_value), NA, abs(segment_value) >= 0.2)]
  segments[, .(
    number_of_segments = .N,
    fraction_segments_altered = mean(altered_segment, na.rm = TRUE),
    copy_number_complexity_score = .N,
    fraction_genome_altered = if (all(is.na(altered_segment)) || sum(segment_length, na.rm = TRUE) == 0) NA_real_ else
      sum(segment_length[altered_segment %in% TRUE], na.rm = TRUE) / sum(segment_length, na.rm = TRUE)
  ), by = .(sample_barcode, patient_barcode, project_id, project_code)]
}

combine_first_nonmissing <- function(base, update, value_columns, source_label) {
  if (is.null(update) || nrow(update) == 0) {
    return(base)
  }
  update <- update[!is.na(sample_barcode)]
  update <- update[order(sample_barcode)]
  update <- update[, lapply(.SD, function(x) x[which.max(!is.na(x))[1]]), by = sample_barcode, .SDcols = setdiff(names(update), "sample_barcode")]
  merged <- merge(base, update, by = "sample_barcode", all.x = TRUE, suffixes = c("", ".new"), sort = FALSE)
  for (column in value_columns) {
    new_col <- paste0(column, ".new")
    if (new_col %in% names(merged)) {
      merged[is.na(get(column)) & !is.na(get(new_col)), (column) := get(new_col)]
    }
  }
  for (column in c("patient_barcode", "project_id", "project_code")) {
    new_col <- paste0(column, ".new")
    if (new_col %in% names(merged)) {
      merged[is.na(get(column)) & !is.na(get(new_col)), (column) := get(new_col)]
    }
  }
  if ("source.new" %in% names(merged)) {
    merged[!is.na(source.new) & source.new != "", source := fifelse(source == "no_source_available", source.new, paste(source, source.new, sep = ";"))]
  } else if (!is.null(source_label)) {
    merged[rowSums(!is.na(.SD)) > 0, source := source_label, .SDcols = value_columns]
  }
  if ("notes.new" %in% names(merged)) {
    merged[!is.na(notes.new) & notes.new != "", notes := fifelse(
      grepl("^no_.*_source_found$", notes),
      notes.new,
      paste(notes, notes.new, sep = ";")
    )]
  }
  drop_cols <- grep("\\.new$", names(merged), value = TRUE)
  if (length(drop_cols) > 0) {
    merged[, (drop_cols) := NULL]
  }
  merged
}

build_copy_number_features <- function(root, projects, sample_universe, resource_files) {
  purity <- copy(sample_universe)
  purity[, `:=`(
    purity = NA_real_,
    ploidy = NA_real_,
    genome_doublings = NA_real_,
    whole_genome_doubling_status = NA_character_,
    subclonal_genome_fraction = NA_real_,
    source = "no_source_available",
    notes = "no_purity_ploidy_source_found"
  )]

  aneuploidy <- copy(sample_universe)
  for (column in setdiff(copy_number_numeric_columns, c("purity", "ploidy"))) {
    aneuploidy[, (column) := NA_real_]
  }
  aneuploidy[, `:=`(source = "no_source_available", notes = "no_copy_number_source_found")]

  arm_calls_all <- data.table()
  segments_all <- data.table()
  source_records <- data.table(source = character(), resource_type = character(), rows = integer())
  parsing_warnings <- data.table(source = character(), warning = character())

  for (path in resource_files) {
    warnings <- character()
    table <- tryCatch(
      withCallingHandlers(
        read_table_auto(path),
        warning = function(warn) {
          warnings <<- c(warnings, conditionMessage(warn))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(err) {
        log_warn(sprintf("Skipping unreadable copy-number resource %s: %s", path, conditionMessage(err)))
        NULL
      }
    )
    if (length(warnings) > 0) {
      parsing_warnings <- rbind(
        parsing_warnings,
        data.table(source = basename(path), warning = unique(warnings)),
        use.names = TRUE,
        fill = TRUE
      )
    }
    if (is.null(table) || nrow(table) == 0) {
      next
    }
    source_name <- basename(path)
    type <- classify_resource(path)

    purity_piece <- if (identical(source_name, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")) {
      extract_absolute_purity_ploidy(table, source_name, sample_universe)
    } else {
      extract_purity_ploidy(table, source_name, sample_universe)
    }
    if (!is.null(purity_piece) && nrow(purity_piece) > 0) {
      purity <- combine_first_nonmissing(
        purity,
        purity_piece,
        c("purity", "ploidy", "genome_doublings", "whole_genome_doubling_status", "subclonal_genome_fraction"),
        source_name
      )
      source_records <- rbind(source_records, data.table(source = source_name, resource_type = "purity_ploidy", rows = nrow(purity_piece)))
    }

    direct_piece <- if (identical(source_name, "PANCAN_ArmCallsAndAneuploidyScore_092817.txt")) {
      extract_pancan_aneuploidy_scores(table, source_name, sample_universe)
    } else {
      extract_direct_aneuploidy(table, source_name, sample_universe)
    }
    if (!is.null(direct_piece) && nrow(direct_piece) > 0) {
      aneuploidy <- combine_first_nonmissing(
        aneuploidy,
        direct_piece,
        setdiff(copy_number_numeric_columns, c("purity", "ploidy")),
        source_name
      )
      source_records <- rbind(source_records, data.table(source = source_name, resource_type = "direct_aneuploidy_or_fga", rows = nrow(direct_piece)))
    }

    arm_piece <- if (identical(source_name, "PANCAN_ArmCallsAndAneuploidyScore_092817.txt")) {
      extract_wide_arm_calls(table, source_name, sample_universe)
    } else {
      extract_arm_calls(table, source_name, sample_universe)
    }
    if (!is.null(arm_piece) && nrow(arm_piece) > 0) {
      arm_calls_all <- data.table::rbindlist(list(arm_calls_all, arm_piece), use.names = TRUE, fill = TRUE)
      source_records <- rbind(source_records, data.table(source = source_name, resource_type = "arm_level_calls", rows = nrow(arm_piece)))
    }

    segment_piece <- extract_segments(table, source_name, sample_universe)
    if (!is.null(segment_piece) && nrow(segment_piece) > 0) {
      segments_all <- data.table::rbindlist(list(segments_all, segment_piece), use.names = TRUE, fill = TRUE)
      source_records <- rbind(source_records, data.table(source = source_name, resource_type = "segments", rows = nrow(segment_piece)))
    }

    if (identical(type, "unknown") && is.null(purity_piece) && is.null(direct_piece) && is.null(arm_piece) && is.null(segment_piece)) {
      log_warn(sprintf("Local resource did not match supported copy-number schemas: %s", path))
    }
  }

  arm_summary <- summarize_arm_calls(arm_calls_all)
  if (!is.null(arm_summary)) {
    aneuploidy <- combine_first_nonmissing(
      aneuploidy,
      arm_summary[, `:=`(source = "derived_from_arm_level_calls", notes = "arm_counts_derived_from_local_arm_calls")],
      c("arm_gain_count", "arm_loss_count", "total_arm_alteration_count"),
      "derived_from_arm_level_calls"
    )
  }

  segment_summary <- summarize_segments(segments_all)
  if (!is.null(segment_summary)) {
    aneuploidy <- combine_first_nonmissing(
      aneuploidy,
      segment_summary[, `:=`(source = "derived_from_segments", notes = "copy_number_complexity_and_fga_approximate_from_segments")],
      c("fraction_genome_altered", "copy_number_complexity_score"),
      "derived_from_segments"
    )
  }

  aneuploidy[is.na(total_arm_alteration_count) & (!is.na(arm_gain_count) | !is.na(arm_loss_count)),
             total_arm_alteration_count := fifelse(is.na(arm_gain_count), 0, arm_gain_count) + fifelse(is.na(arm_loss_count), 0, arm_loss_count)]

  if (length(resource_files) == 0) {
    log_warn("No local copy-number, purity, ploidy, arm-level, or segment resources were found. Writing mutation-sample-aligned tables with missing feature values.")
  }

  list(
    purity = purity,
    aneuploidy = aneuploidy,
    arm_calls = arm_calls_all,
    segments = segments_all,
    source_records = source_records,
    parsing_warnings = parsing_warnings
  )
}

summarize_project_features <- function(purity, aneuploidy, projects) {
  combined <- merge(
    purity[, .(sample_barcode, project_id, project_code, purity, ploidy)],
    aneuploidy[, .(
      sample_barcode,
      aneuploidy_score,
      fraction_genome_altered,
      arm_gain_count,
      arm_loss_count,
      total_arm_alteration_count,
      copy_number_complexity_score
    )],
    by = "sample_barcode",
    all = TRUE,
    sort = FALSE
  )
  metrics <- copy_number_numeric_columns
  summary <- combined[, c(
    list(n_samples = .N),
    lapply(.SD, function(x) stats::median(x, na.rm = TRUE)),
    setNames(lapply(.SD, function(x) mean(is.na(x))), paste0("missingness_", metrics))
  ), by = .(project_id, project_code), .SDcols = metrics]
  for (metric in metrics) {
    data.table::setnames(summary, metric, paste0("median_", metric))
    summary[is.infinite(get(paste0("median_", metric))), (paste0("median_", metric)) := NA_real_]
  }
  summary <- merge(projects, summary, by = c("project_id", "project_code"), all.x = TRUE, sort = FALSE)
  summary[is.na(n_samples), n_samples := 0L]
  data.table::setorder(summary, project_id)
  summary
}

write_qc_summary <- function(purity, aneuploidy, segments, project_summary, mutation_samples, projects, source_records, retrieval_records, parsing_warnings, output_path) {
  cn_metric_columns <- setdiff(copy_number_numeric_columns, c("purity", "ploidy"))
  aneuploidy_has_data <- rowSums(!is.na(aneuploidy[, ..cn_metric_columns])) > 0
  represented_projects <- sort(unique(c(
    purity[!is.na(purity) | !is.na(ploidy), project_id],
    aneuploidy[aneuploidy_has_data, project_id],
    if (nrow(segments) == 0) character() else segments[, project_id]
  )))
  represented_projects <- represented_projects[!is.na(represented_projects)]
  missing_projects <- setdiff(projects$project_id, represented_projects)

  purity_samples <- purity[!is.na(purity), uniqueN(sample_barcode)]
  ploidy_samples <- purity[!is.na(ploidy), uniqueN(sample_barcode)]
  pp_samples <- purity[!is.na(purity) | !is.na(ploidy), uniqueN(sample_barcode)]
  aneuploidy_score_samples <- aneuploidy[!is.na(aneuploidy_score), uniqueN(sample_barcode)]
  cn_samples <- aneuploidy[aneuploidy_has_data, uniqueN(sample_barcode)]
  segment_samples <- if (nrow(segments) == 0) 0L else segments[, uniqueN(sample_barcode)]
  arm_call_samples <- source_records[resource_type == "arm_level_calls", {
    if (nrow(source_records[resource_type == "arm_level_calls"]) == 0) 0L else NA_integer_
  }]
  arm_call_samples <- if (nrow(source_records[resource_type == "arm_level_calls"]) == 0) {
    0L
  } else {
    aneuploidy[!is.na(arm_gain_count) | !is.na(arm_loss_count) | !is.na(total_arm_alteration_count), uniqueN(sample_barcode)]
  }
  mutation_sample_count <- mutation_samples[, uniqueN(sample_barcode)]
  pp_overlap <- purity[(!is.na(purity) | !is.na(ploidy)) & sample_barcode %in% mutation_samples$sample_barcode, uniqueN(sample_barcode)]
  cn_overlap <- aneuploidy[aneuploidy_has_data & sample_barcode %in% mutation_samples$sample_barcode, uniqueN(sample_barcode)]
  parsing_warning_value <- if (nrow(parsing_warnings) == 0) "none" else paste(sprintf("%s: %s", parsing_warnings$source, parsing_warnings$warning), collapse = "; ")

  absolute_status <- retrieval_records[resource_key == "absolute_purity_ploidy", status]
  arm_status <- retrieval_records[resource_key == "arm_calls_aneuploidy", status]

  global <- data.table(
    qc_section = "global",
    project_id = "ALL",
    project_code = "ALL",
    metric = c(
      "absolute_purity_ploidy_file_status",
      "arm_calls_aneuploidy_file_status",
      "samples_with_purity",
      "samples_with_ploidy",
      "samples_with_aneuploidy_score",
      "samples_with_arm_level_calls",
      "samples_with_purity_ploidy_data",
      "samples_with_aneuploidy_or_copy_number_data",
      "samples_with_copy_number_segment_data",
      "tcga_projects_with_copy_number_or_purity_data",
      "expected_tcga_projects_missing_copy_number_or_purity",
      "missing_tcga_projects_copy_number_or_purity",
      "mutation_layer_samples",
      "mutation_layer_samples_overlapping_purity_ploidy",
      "mutation_layer_samples_overlapping_copy_number_features",
      "percent_mutation_layer_samples_with_purity_ploidy",
      "percent_mutation_layer_samples_with_copy_number_features",
      "local_or_downloaded_source_files_used",
      "parsing_warnings"
    ),
    value = c(
      ifelse(length(absolute_status) == 0, "unavailable", absolute_status),
      ifelse(length(arm_status) == 0, "unavailable", arm_status),
      as.character(purity_samples),
      as.character(ploidy_samples),
      as.character(aneuploidy_score_samples),
      as.character(arm_call_samples),
      as.character(pp_samples),
      as.character(cn_samples),
      as.character(segment_samples),
      as.character(length(represented_projects)),
      as.character(length(missing_projects)),
      ifelse(length(missing_projects) == 0, "none", paste(missing_projects, collapse = ";")),
      as.character(mutation_sample_count),
      as.character(pp_overlap),
      as.character(cn_overlap),
      sprintf("%.4f", ifelse(mutation_sample_count == 0, NA_real_, pp_overlap / mutation_sample_count * 100)),
      sprintf("%.4f", ifelse(mutation_sample_count == 0, NA_real_, cn_overlap / mutation_sample_count * 100)),
      ifelse(nrow(source_records) == 0, "none", paste(unique(source_records$source), collapse = ";")),
      parsing_warning_value
    )
  )

  per_project_metrics <- c(
    "n_samples",
    "median_purity",
    "median_ploidy",
    "median_aneuploidy_score",
    "median_fraction_genome_altered",
    "missingness_purity",
    "missingness_ploidy",
    "missingness_aneuploidy_score",
    "missingness_fraction_genome_altered",
    "missingness_arm_gain_count",
    "missingness_arm_loss_count",
    "missingness_total_arm_alteration_count"
  )
  per_project_input <- copy(project_summary[, c("project_id", "project_code", per_project_metrics), with = FALSE])
  for (metric in per_project_metrics) {
    per_project_input[, (metric) := as.character(get(metric))]
  }
  per_project <- data.table::melt(
    per_project_input,
    id.vars = c("project_id", "project_code"),
    variable.name = "metric",
    value.name = "value"
  )
  per_project[, qc_section := "per_project"]

  retrieval_qc <- retrieval_records[, .(
    qc_section = "sources",
    project_id = "ALL",
    project_code = "ALL",
    metric = paste0(resource_key, "_status"),
    value = paste(status, path, sep = ":")
  )]

  source_qc <- if (nrow(source_records) == 0) {
    data.table(qc_section = "sources", project_id = "ALL", project_code = "ALL", metric = "parsed_source_files", value = "none")
  } else {
    source_records[, .(qc_section = "sources", project_id = "ALL", project_code = "ALL", metric = paste(resource_type, source, sep = ":"), value = as.character(rows))]
  }

  qc <- data.table::rbindlist(
    list(global, per_project[, .(qc_section, project_id, project_code, metric, value)], retrieval_qc, source_qc),
    use.names = TRUE,
    fill = TRUE
  )
  data.table::setorder(qc, qc_section, project_id, metric)
  write_tsv_table(qc, output_path)

  list(
    pp_samples = pp_samples,
    purity_samples = purity_samples,
    ploidy_samples = ploidy_samples,
    aneuploidy_score_samples = aneuploidy_score_samples,
    arm_call_samples = arm_call_samples,
    cn_samples = cn_samples,
    segment_samples = segment_samples,
    represented_projects = length(represented_projects),
    missing_projects = missing_projects,
    pp_overlap = pp_overlap,
    cn_overlap = cn_overlap,
    mutation_sample_count = mutation_sample_count
  )
}

write_qc_figure <- function(project_summary, output_path) {
  ensure_dir(dirname(output_path))
  grDevices::pdf(output_path, width = 10, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(8, 4, 3, 1))
  values <- 1 - project_summary$missingness_purity
  values[is.na(values)] <- 0
  names(values) <- project_summary$project_code
  graphics::barplot(
    values,
    las = 2,
    ylim = c(0, 1),
    ylab = "Fraction with purity data",
    main = "Copy-number layer QC: purity/ploidy coverage"
  )
  log_info(sprintf("Wrote %s", output_path))
}

validate_outputs_nonempty <- function(paths) {
  for (path in paths) {
    if (!file.exists(path) || file.info(path)$size == 0) {
      abort(sprintf("Expected non-empty output was not created: %s", path))
    }
  }
}

main <- function() {
  args <- parse_args()
  require_copy_number_packages()
  root <- find_project_root()

  projects_path <- project_path("data", "interim", "tcga_projects.tsv", root = root)
  mutation_path <- project_path("data", "processed", "mutations", "mc3_somatic_mutations.parquet", root = root)
  projects <- read_projects(projects_path)
  mutation_samples <- read_mutation_sample_universe(mutation_path)
  log_info(sprintf("Loaded %d mutation-layer samples as copy-number sample universe", nrow(mutation_samples)))

  retrieval_records <- retrieve_pancanatlas_resources(root, refresh = args$refresh)
  resource_files <- discover_local_resources(root)
  log_info(sprintf("Discovered %d local copy-number/purity resource files", length(resource_files)))

  features <- build_copy_number_features(root, projects, mutation_samples, resource_files)
  project_summary <- summarize_project_features(features$purity, features$aneuploidy, projects)

  purity_output <- project_path("data", "processed", "features", "purity_ploidy_by_sample.tsv", root = root)
  aneuploidy_output <- project_path("data", "processed", "features", "aneuploidy_by_sample.tsv", root = root)
  project_output <- project_path("data", "processed", "features", "aneuploidy_by_project.tsv", root = root)
  arm_output <- project_path("data", "processed", "copy_number", "arm_level_calls_by_sample.tsv", root = root)
  segment_output <- project_path("data", "processed", "copy_number", "segments_by_sample.tsv", root = root)
  qc_output <- project_path("results", "tables", "copy_number_layer_qc_summary.tsv", root = root)
  figure_output <- project_path("results", "figures", "copy_number_layer_qc_overview.pdf", root = root)

  write_tsv_table(features$purity, purity_output)
  write_tsv_table(features$aneuploidy, aneuploidy_output)
  write_tsv_table(project_summary, project_output)

  arm_out <- features$arm_calls
  if (nrow(arm_out) == 0) {
    arm_out <- data.table(
      sample_barcode = character(),
      patient_barcode = character(),
      project_id = character(),
      project_code = character(),
      chromosome_arm = character(),
      arm_call = character(),
      source = character()
    )
  }
  write_tsv_table(arm_out, arm_output)

  segment_out <- features$segments
  if (nrow(segment_out) == 0) {
    segment_out <- data.table(
      sample_barcode = character(),
      patient_barcode = character(),
      project_id = character(),
      project_code = character(),
      chromosome = character(),
      start = numeric(),
      end = numeric(),
      segment_value = numeric(),
      source = character()
    )
  }
  write_tsv_table(segment_out, segment_output)

  qc <- write_qc_summary(
    features$purity,
    features$aneuploidy,
    features$segments,
    project_summary,
    mutation_samples,
    projects,
    features$source_records,
    retrieval_records,
    features$parsing_warnings,
    qc_output
  )
  write_qc_figure(project_summary, figure_output)
  validate_outputs_nonempty(c(purity_output, aneuploidy_output, project_output, qc_output, figure_output))

  log_info(sprintf("Copy-number layer complete: %d projects with any sourced purity/CN data", qc$represented_projects))
  log_info(sprintf("Samples with purity/ploidy data: %d", qc$pp_samples))
  log_info(sprintf("Samples with aneuploidy/copy-number data: %d", qc$cn_samples))
  log_info(sprintf("Overlap with mutation-layer samples: purity/ploidy=%d, copy-number=%d of %d", qc$pp_overlap, qc$cn_overlap, qc$mutation_sample_count))
}

if (identical(environment(), globalenv())) {
  main()
}
