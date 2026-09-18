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

expected_expression_gene_sets <- c(
  "proliferation",
  "immune_inflammatory",
  "interferon_gamma",
  "cytotoxic_t_cell",
  "stromal",
  "epithelial",
  "EMT",
  "hypoxia",
  "cell_cycle",
  "DNA_repair",
  "angiogenesis"
)

expression_score_map <- data.table::data.table(
  gene_set = expected_expression_gene_sets,
  score_column = c(
    "proliferation_score",
    "immune_inflammatory_score",
    "interferon_gamma_score",
    "cytotoxic_t_cell_score",
    "stromal_score",
    "epithelial_score",
    "EMT_score",
    "hypoxia_score",
    "cell_cycle_score",
    "DNA_repair_score",
    "angiogenesis_score"
  )
)

expression_score_columns <- expression_score_map$score_column
lineage_pc_columns <- c("lineage_or_tissue_PC1", "lineage_or_tissue_PC2")
pca_columns <- paste0("PC", seq_len(20))
tumor_sample_type_codes <- c("01", "02", "03", "05", "06", "07", "08", "09")

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  list(refresh = "--refresh" %in% args)
}

require_expression_packages <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    abort("Package 'data.table' is required for expression processing")
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    abort("Package 'yaml' is required for expression gene-set configuration")
  }
  suppressPackageStartupMessages(library(data.table))
}

write_tsv_table <- function(data, path) {
  ensure_dir(dirname(path))
  data.table::fwrite(data, path, sep = "\t", quote = FALSE, na = "NA")
  log_info(sprintf("Wrote %s", path))
  invisible(path)
}

write_tsv_gz_table <- function(data, path) {
  ensure_dir(dirname(path))
  data.table::fwrite(data, path, sep = "\t", quote = FALSE, na = "NA", compress = "auto")
  log_info(sprintf("Wrote %s", path))
  invisible(path)
}

read_projects <- function(path) {
  require_file(path, "TCGA project metadata")
  projects <- data.table::fread(path, sep = "\t", data.table = TRUE)
  required <- c("project_id", "project_code")
  missing <- setdiff(required, names(projects))
  if (length(missing) > 0) {
    abort(sprintf("Project metadata is missing required columns: %s", paste(missing, collapse = ", ")))
  }
  unique(projects[, ..required])
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

read_mutation_sample_universe <- function(mutation_path, tmb_path = NA_character_) {
  if (!is.na(tmb_path) && file.exists(tmb_path)) {
    table <- data.table::fread(tmb_path, sep = "\t", data.table = TRUE)
  } else {
    require_file(mutation_path, "mutation-layer parquet")
    table <- read_parquet_with_fallback(mutation_path)
  }

  required <- c("sample_barcode", "patient_barcode", "project_id", "project_code")
  missing <- setdiff(required, names(table))
  if (length(missing) > 0) {
    abort(sprintf("Mutation sample universe is missing required columns: %s", paste(missing, collapse = ", ")))
  }

  universe <- unique(table[, ..required])
  data.table::setorder(universe, project_id, sample_barcode)
  universe
}

read_sample_table_optional <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.table(sample_barcode = character()))
  }
  table <- data.table::fread(path, sep = "\t", data.table = TRUE)
  if (!"sample_barcode" %in% names(table)) {
    return(data.table(sample_barcode = character()))
  }
  unique(table[!is.na(sample_barcode), .(sample_barcode)])
}

pancanatlas_expression_spec <- function(root) {
  data.table::data.table(
    resource_key = "pancanatlas_ebplusplus_rnaseqv2",
    expected_file = "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",
    local_path = project_path("data", "raw", "expression", "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", root = root),
    gdc_url = "https://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611",
    publication_page = "https://gdc.cancer.gov/about-data/publications/pancanatlas"
  )
}

download_public_resource <- function(url, destination) {
  ensure_dir(dirname(destination))
  tmp <- tempfile(pattern = "pancanatlas_expression_", fileext = ".download")
  on.exit(unlink(tmp), add = TRUE)
  old_timeout <- getOption("timeout")
  options(timeout = max(3600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  tryCatch(
    {
      utils::download.file(url, tmp, mode = "wb", quiet = FALSE, method = "auto")
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

manual_expression_download_message <- function(spec) {
  paste(
    "Automatic download failed for the PanCanAtlas expression resource.",
    sprintf("Manual fallback: download %s from %s and place it at %s", spec$expected_file, spec$publication_page, spec$local_path),
    sep = "\n"
  )
}

retrieve_expression_resource <- function(root, refresh = FALSE) {
  spec <- pancanatlas_expression_spec(root)
  ensure_dir(dirname(spec$local_path))
  spec[, `:=`(status = "unavailable", path = local_path, file_size = NA_real_)]

  if (file.exists(spec$local_path) && file.info(spec$local_path)$size > 0 && !refresh) {
    spec[, `:=`(status = "local", file_size = file.info(local_path)$size)]
    log_info(sprintf("Using local PanCanAtlas expression resource: %s", spec$local_path))
    return(spec)
  }

  log_info(sprintf("Attempting download of %s from GDC", spec$expected_file))
  downloaded <- download_public_resource(spec$gdc_url, spec$local_path)
  if (downloaded && file.exists(spec$local_path) && file.info(spec$local_path)$size > 0) {
    spec[, `:=`(status = "downloaded", file_size = file.info(local_path)$size)]
    log_info(sprintf("Downloaded PanCanAtlas expression resource: %s", spec$local_path))
  } else {
    log_warn(manual_expression_download_message(spec))
  }
  spec
}

discover_local_expression_resource <- function(root, retrieval_record) {
  raw_dir <- project_path("data", "raw", "expression", root = root)
  ensure_dir(raw_dir)
  expected <- retrieval_record$local_path[[1]]
  if (file.exists(expected) && file.info(expected)$size > 0) {
    return(expected)
  }

  candidates <- list.files(
    raw_dir,
    pattern = "\\.(tsv|tsv\\.gz|txt|txt\\.gz|csv|csv\\.gz)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  candidates <- candidates[basename(candidates) != ".gitkeep"]
  if (length(candidates) == 0) {
    return(NA_character_)
  }

  preferred <- candidates[grepl("geneexp|expression|rnaseq|fpkm|tpm|count", basename(candidates), ignore.case = TRUE)]
  if (length(preferred) > 0) {
    return(preferred[[1]])
  }
  candidates[[1]]
}

default_gene_set_config_lines <- function() {
  c(
    "gene_sets:",
    "  proliferation:",
    "    genes: [MKI67, TOP2A, PCNA, MCM2, MCM5, CCNB1, CCNA2, CDK1]",
    "  immune_inflammatory:",
    "    genes: [PTPRC, CD3D, CD3E, CD8A, CD8B, GZMB, PRF1, IFNG, CXCL9, CXCL10]",
    "  interferon_gamma:",
    "    genes: [IFNG, STAT1, IRF1, CXCL9, CXCL10, IDO1, GBP1, HLA-DRA]",
    "  cytotoxic_t_cell:",
    "    genes: [CD8A, CD8B, GZMB, PRF1, NKG7, GNLY, GZMA]",
    "  stromal:",
    "    genes: [COL1A1, COL1A2, COL3A1, ACTA2, TAGLN, FAP, PDGFRB]",
    "  epithelial:",
    "    genes: [EPCAM, KRT8, KRT18, KRT19, CDH1]",
    "  EMT:",
    "    genes: [VIM, ZEB1, ZEB2, SNAI1, SNAI2, TWIST1, FN1, ITGA5]",
    "  hypoxia:",
    "    genes: [VEGFA, CA9, SLC2A1, LDHA, PGK1, ENO1]",
    "  cell_cycle:",
    "    genes: [MKI67, TOP2A, PCNA, MCM2, MCM5, CCNB1, CCNA2, CDK1, BUB1, CDC20]",
    "  DNA_repair:",
    "    genes: [BRCA1, BRCA2, RAD51, ATM, ATR, CHEK1, CHEK2, FANCD2]",
    "  angiogenesis:",
    "    genes: [VEGFA, KDR, FLT1, ANGPT2, TEK, PECAM1]"
  )
}

ensure_expression_gene_set_config <- function(path) {
  if (!file.exists(path)) {
    ensure_dir(dirname(path))
    writeLines(default_gene_set_config_lines(), path)
    log_warn(sprintf("Created starter expression gene-set config: %s", path))
  }
  invisible(path)
}

extract_genes_from_gene_set_entry <- function(entry) {
  if (is.null(entry)) {
    return(character())
  }
  genes <- if (is.list(entry) && !is.null(entry$genes)) entry$genes else entry
  genes <- toupper(trimws(as.character(unlist(genes, use.names = FALSE))))
  unique(genes[!is.na(genes) & genes != ""])
}

read_expression_gene_sets <- function(path) {
  ensure_expression_gene_set_config(path)
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg$gene_sets)) {
    abort(sprintf("Expression gene-set config must contain top-level key 'gene_sets': %s", path))
  }

  gene_sets <- lapply(expected_expression_gene_sets, function(set_name) {
    extract_genes_from_gene_set_entry(cfg$gene_sets[[set_name]])
  })
  names(gene_sets) <- expected_expression_gene_sets

  missing_sets <- expected_expression_gene_sets[vapply(gene_sets, length, integer(1)) == 0L]
  if (length(missing_sets) > 0) {
    log_warn(sprintf("Expression gene-set config has empty or missing gene sets: %s", paste(missing_sets, collapse = ", ")))
  }
  gene_sets
}

strip_quotes <- function(x) {
  x <- sub("^\ufeff", "", x)
  gsub('^"|"$', "", x)
}

delimiter_for_file <- function(path) {
  if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
}

open_text_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(gzfile(path, open = "rt"))
  }
  file(path, open = "rt")
}

clean_expression_gene_symbol <- function(gene_id) {
  gene_id <- strip_quotes(trimws(as.character(gene_id)))
  gene_id <- sub("\\.[0-9]+$", "", gene_id)
  before_pipe <- sub("\\|.*$", "", gene_id)
  after_pipe <- sub("^.*\\|", "", gene_id)
  symbol <- ifelse(before_pipe %in% c("", "?", "NA"), after_pipe, before_pipe)
  toupper(trimws(symbol))
}

normalize_tcga_barcode_text <- function(x) {
  x <- strip_quotes(trimws(as.character(x)))
  x <- gsub("\\.", "-", x)
  x
}

harmonize_expression_sample_metadata <- function(sample_header, sample_universe) {
  metadata <- data.table::as.data.table(sample_header)
  metadata[, raw_barcode := normalize_tcga_barcode_text(column_name)]
  metadata[, sample_barcode := substr(raw_barcode, 1, 16)]
  metadata[, patient_barcode := substr(raw_barcode, 1, 12)]
  metadata[, sample_type_code := substr(sample_barcode, 14, 15)]
  metadata[, is_tumor_sample := sample_type_code %in% tumor_sample_type_codes]

  universe_exact <- sample_universe[, .(
    sample_barcode,
    mapped_patient_barcode = patient_barcode,
    mapped_project_id = project_id,
    mapped_project_code = project_code
  )]
  metadata <- merge(metadata, universe_exact, by = "sample_barcode", all.x = TRUE, sort = FALSE)
  metadata[!is.na(mapped_patient_barcode), patient_barcode := mapped_patient_barcode]
  metadata[, `:=`(project_id = mapped_project_id, project_code = mapped_project_code)]
  metadata[, c("mapped_patient_barcode", "mapped_project_id", "mapped_project_code") := NULL]

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

  metadata[, sample_prefix15 := substr(sample_barcode, 1, 15)]
  metadata <- merge(metadata, prefix_map, by = "sample_prefix15", all.x = TRUE, sort = FALSE)
  metadata[(is.na(project_id) | !(sample_barcode %in% sample_universe$sample_barcode)) & !is.na(prefix_sample_barcode),
           sample_barcode := prefix_sample_barcode]
  metadata[is.na(project_id) & !is.na(prefix_project_id), project_id := prefix_project_id]
  metadata[is.na(project_code) & !is.na(prefix_project_code), project_code := prefix_project_code]
  metadata[(is.na(patient_barcode) | patient_barcode == "") & !is.na(prefix_patient_barcode), patient_barcode := prefix_patient_barcode]
  metadata[, c("sample_prefix15", "prefix_sample_barcode", "prefix_patient_barcode", "prefix_project_id", "prefix_project_code") := NULL]

  metadata[is_tumor_sample %in% TRUE & !is.na(project_id) & !is.na(sample_barcode)]
}

collapse_expression_columns_by_sample <- function(matrix_values, metadata) {
  if (ncol(matrix_values) == 0 || nrow(metadata) == 0) {
    return(list(
      matrix = matrix(numeric(), nrow = nrow(matrix_values), ncol = 0, dimnames = list(rownames(matrix_values), character())),
      metadata = metadata[0]
    ))
  }

  groups <- split(seq_len(ncol(matrix_values)), metadata$sample_barcode)
  sample_order <- names(groups)
  collapsed <- vapply(groups, function(idx) {
    values <- rowMeans(matrix_values[, idx, drop = FALSE], na.rm = TRUE)
    values[is.nan(values)] <- NA_real_
    values
  }, numeric(nrow(matrix_values)))
  if (nrow(matrix_values) == 1L) {
    collapsed <- matrix(collapsed, nrow = 1, dimnames = list(rownames(matrix_values), sample_order))
  }
  colnames(collapsed) <- sample_order

  metadata_unique <- metadata[match(sample_order, sample_barcode)]
  metadata_unique <- metadata_unique[, .(sample_barcode, patient_barcode, project_id, project_code, sample_type_code)]
  list(matrix = collapsed, metadata = metadata_unique)
}

collapse_duplicate_genes <- function(matrix_values) {
  if (nrow(matrix_values) == 0) {
    return(matrix_values)
  }
  gene_groups <- split(seq_len(nrow(matrix_values)), rownames(matrix_values))
  collapsed <- t(vapply(gene_groups, function(idx) {
    values <- colMeans(matrix_values[idx, , drop = FALSE], na.rm = TRUE)
    values[is.nan(values)] <- NA_real_
    values
  }, numeric(ncol(matrix_values))))
  rownames(collapsed) <- names(gene_groups)
  colnames(collapsed) <- colnames(matrix_values)
  collapsed
}

read_expression_subset <- function(path, gene_sets, sample_universe, chunk_size = 1000L) {
  require_file(path, "expression matrix")
  delimiter <- delimiter_for_file(path)
  gene_universe <- unique(unlist(gene_sets, use.names = FALSE))
  gene_universe <- toupper(gene_universe[!is.na(gene_universe) & gene_universe != ""])

  con <- open_text_connection(path)
  on.exit(close(con), add = TRUE)
  header_line <- readLines(con, n = 1L, warn = FALSE)
  if (length(header_line) == 0) {
    abort(sprintf("Expression matrix has no header: %s", path))
  }
  header <- strip_quotes(strsplit(header_line, delimiter, fixed = TRUE)[[1]])
  gene_col <- match(TRUE, tolower(header) %in% c("gene_id", "gene", "gene_symbol", "hugo_symbol", "hugo_symbol_or_gene_id"), nomatch = 0L)
  if (gene_col == 0L) {
    gene_col <- 1L
  }

  sample_positions <- setdiff(seq_along(header), gene_col)
  sample_header <- data.table::data.table(header_position = sample_positions, column_name = header[sample_positions])
  metadata <- harmonize_expression_sample_metadata(sample_header, sample_universe)
  metadata <- metadata[!duplicated(paste(header_position, sample_barcode))]
  selected_positions <- metadata$header_position

  if (length(selected_positions) == 0) {
    log_warn("No expression columns matched tumor samples in the mutation-layer sample universe")
    return(list(
      expression_matrix = matrix(numeric(), nrow = 0, ncol = 0, dimnames = list(character(), character())),
      sample_metadata = data.table(sample_barcode = character(), patient_barcode = character(), project_id = character(), project_code = character(), sample_type_code = character()),
      raw_expression_sample_count = length(sample_positions),
      tumor_expression_column_count = 0L,
      parsing_warnings = "no_matching_tumor_samples"
    ))
  }

  log_info(sprintf("Scanning expression matrix for %d configured genes across %d tumor sample columns", length(gene_universe), length(selected_positions)))
  row_values <- list()
  row_genes <- character()
  parsing_warnings <- character()
  rows_scanned <- 0L

  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (length(lines) == 0) {
      break
    }
    for (line in lines) {
      rows_scanned <- rows_scanned + 1L
      fields <- strsplit(line, delimiter, fixed = TRUE)[[1]]
      if (length(fields) < max(selected_positions)) {
        parsing_warnings <- c(parsing_warnings, sprintf("short_row_%d", rows_scanned))
        next
      }
      gene_symbol <- clean_expression_gene_symbol(fields[[gene_col]])
      if (!gene_symbol %in% gene_universe) {
        next
      }
      values <- suppressWarnings(as.numeric(strip_quotes(fields[selected_positions])))
      row_values[[length(row_values) + 1L]] <- values
      row_genes <- c(row_genes, gene_symbol)
    }
    if (rows_scanned %% 5000L == 0L) {
      log_info(sprintf("Scanned %d expression rows; found %d configured gene rows", rows_scanned, length(row_values)))
    }
  }

  if (length(row_values) == 0) {
    log_warn("No configured expression gene-set genes were found in the expression matrix")
    return(list(
      expression_matrix = matrix(numeric(), nrow = 0, ncol = 0, dimnames = list(character(), character())),
      sample_metadata = metadata[, .(sample_barcode, patient_barcode, project_id, project_code, sample_type_code)][0],
      raw_expression_sample_count = length(sample_positions),
      tumor_expression_column_count = nrow(metadata),
      parsing_warnings = ifelse(length(parsing_warnings) == 0, "no_configured_genes_found", paste(unique(parsing_warnings), collapse = ";"))
    ))
  }

  matrix_values <- do.call(rbind, row_values)
  rownames(matrix_values) <- row_genes
  colnames(matrix_values) <- metadata$column_name
  matrix_values <- collapse_duplicate_genes(matrix_values)
  collapsed <- collapse_expression_columns_by_sample(matrix_values, metadata)
  collapsed$matrix <- collapse_duplicate_genes(collapsed$matrix)

  list(
    expression_matrix = collapsed$matrix,
    sample_metadata = collapsed$metadata,
    raw_expression_sample_count = length(sample_positions),
    tumor_expression_column_count = nrow(metadata),
    parsing_warnings = ifelse(length(parsing_warnings) == 0, "none", paste(unique(parsing_warnings), collapse = ";"))
  )
}

detect_expression_data_type <- function(path) {
  name <- tolower(basename(path))
  if (grepl("count|htseq", name)) {
    return("counts")
  }
  if (grepl("tpm", name)) {
    return("TPM")
  }
  if (grepl("fpkm", name)) {
    return("FPKM")
  }
  if (grepl("ebplusplus|rnaseqv2|geneexp", name)) {
    return("pancanatlas_ebplusplus_adjusted_rnaseqv2")
  }
  "unknown_normalized_or_transformed"
}

normalize_expression_matrix <- function(expression_matrix, data_type) {
  if (nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    return(list(matrix = expression_matrix, method = "not_applied_no_expression_values"))
  }

  if (identical(data_type, "counts")) {
    library_sizes <- colSums(expression_matrix, na.rm = TRUE)
    library_sizes[library_sizes <= 0 | is.na(library_sizes)] <- NA_real_
    normalized <- sweep(expression_matrix, 2, library_sizes, "/") * 1e6
    normalized <- log2(normalized + 1)
    return(list(matrix = normalized, method = "library_size_normalized_log2_cpm_plus_1"))
  }

  if (data_type %in% c("TPM", "FPKM")) {
    return(list(matrix = log2(expression_matrix + 1), method = sprintf("log2_%s_plus_1", tolower(data_type))))
  }

  list(matrix = expression_matrix, method = "as_provided_gene_level_z_scoring_for_pathway_scores")
}

zscore_gene_matrix <- function(expression_matrix) {
  if (nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    return(expression_matrix)
  }
  centers <- rowMeans(expression_matrix, na.rm = TRUE)
  sds <- apply(expression_matrix, 1, stats::sd, na.rm = TRUE)
  z <- sweep(expression_matrix, 1, centers, "-")
  z <- sweep(z, 1, sds, "/")
  z[!is.finite(z)] <- NA_real_
  z[sds == 0 | is.na(sds), ] <- NA_real_
  z
}

score_gene_sets <- function(normalized_matrix, sample_metadata, gene_sets) {
  scores <- copy(sample_metadata[, .(sample_barcode, patient_barcode, project_id, project_code)])
  z_matrix <- zscore_gene_matrix(normalized_matrix)

  for (i in seq_len(nrow(expression_score_map))) {
    gene_set <- expression_score_map$gene_set[[i]]
    score_column <- expression_score_map$score_column[[i]]
    found <- intersect(gene_sets[[gene_set]], rownames(z_matrix))
    if (length(found) == 0 || ncol(z_matrix) == 0) {
      scores[, (score_column) := NA_real_]
      next
    }
    values <- colMeans(z_matrix[found, , drop = FALSE], na.rm = TRUE)
    values[is.nan(values)] <- NA_real_
    scores[, (score_column) := values[match(sample_barcode, names(values))]]
  }
  scores
}

gene_set_coverage <- function(gene_sets, expression_matrix) {
  available_genes <- rownames(expression_matrix)
  coverage <- expression_score_map[, .(gene_set, score_column)]
  coverage[, configured_genes := vapply(gene_set, function(set_name) paste(gene_sets[[set_name]], collapse = ";"), character(1))]
  coverage[, n_configured_genes := vapply(gene_set, function(set_name) length(gene_sets[[set_name]]), integer(1))]
  coverage[, found_genes := vapply(gene_set, function(set_name) paste(intersect(gene_sets[[set_name]], available_genes), collapse = ";"), character(1))]
  coverage[, missing_genes := vapply(gene_set, function(set_name) paste(setdiff(gene_sets[[set_name]], available_genes), collapse = ";"), character(1))]
  coverage[, n_found_genes := vapply(gene_set, function(set_name) length(intersect(gene_sets[[set_name]], available_genes)), integer(1))]
  coverage[, n_missing_genes := n_configured_genes - n_found_genes]
  coverage
}

impute_feature_matrix <- function(feature_matrix) {
  keep <- vapply(seq_len(ncol(feature_matrix)), function(j) {
    values <- feature_matrix[, j]
    sum(!is.na(values)) >= 2L && stats::sd(values, na.rm = TRUE) > 0
  }, logical(1))
  feature_matrix <- feature_matrix[, keep, drop = FALSE]
  if (ncol(feature_matrix) == 0) {
    return(feature_matrix)
  }
  for (j in seq_len(ncol(feature_matrix))) {
    values <- feature_matrix[, j]
    replacement <- stats::median(values, na.rm = TRUE)
    if (!is.finite(replacement)) {
      replacement <- 0
    }
    values[is.na(values)] <- replacement
    feature_matrix[, j] <- values
  }
  feature_matrix
}

empty_pca_outputs <- function(sample_metadata, status) {
  pca <- copy(sample_metadata[, .(sample_barcode, patient_barcode, project_id, project_code)])
  for (column in pca_columns) {
    pca[, (column) := NA_real_]
  }
  variance <- data.table::data.table(
    PC = pca_columns,
    percent_variance_explained = NA_real_,
    pca_basis = status
  )
  list(scores = pca, variance = variance, status = status)
}

compute_expression_pca <- function(normalized_matrix, sample_metadata) {
  if (nrow(normalized_matrix) < 2 || ncol(normalized_matrix) < 2 || nrow(sample_metadata) < 2) {
    return(empty_pca_outputs(sample_metadata, "not_run_insufficient_expression_features_or_samples"))
  }

  z_matrix <- zscore_gene_matrix(normalized_matrix)
  feature_matrix <- t(z_matrix)
  rownames(feature_matrix) <- colnames(z_matrix)
  feature_matrix <- feature_matrix[sample_metadata$sample_barcode, , drop = FALSE]
  feature_matrix <- impute_feature_matrix(feature_matrix)
  if (nrow(feature_matrix) < 2 || ncol(feature_matrix) < 2) {
    return(empty_pca_outputs(sample_metadata, "not_run_insufficient_nonconstant_expression_features"))
  }

  pca_fit <- stats::prcomp(feature_matrix, center = TRUE, scale. = FALSE)
  scores <- copy(sample_metadata[, .(sample_barcode, patient_barcode, project_id, project_code)])
  for (column in pca_columns) {
    scores[, (column) := NA_real_]
  }

  available_pcs <- min(ncol(pca_fit$x), length(pca_columns))
  for (i in seq_len(available_pcs)) {
    scores[, (pca_columns[[i]]) := pca_fit$x[, i]]
  }

  variance <- data.table::data.table(
    PC = pca_columns,
    percent_variance_explained = NA_real_,
    pca_basis = "selected_gene_expression_pca"
  )
  variance$percent_variance_explained[seq_len(available_pcs)] <- (pca_fit$sdev[seq_len(available_pcs)]^2 / sum(pca_fit$sdev^2)) * 100

  list(scores = scores, variance = variance, status = "selected_gene_expression_pca")
}

add_missing_score_columns <- function(table) {
  for (column in c(expression_score_columns, lineage_pc_columns)) {
    if (!column %in% names(table)) {
      table[, (column) := NA_real_]
    }
  }
  table
}

missing_expression_outputs <- function(projects, sample_universe, gene_sets, reason) {
  sample_scores <- copy(sample_universe)
  sample_scores <- add_missing_score_columns(sample_scores)
  sample_scores[, `:=`(source = "no_source_available", notes = reason)]

  pca <- copy(sample_universe)
  for (column in pca_columns) {
    pca[, (column) := NA_real_]
  }
  variance <- data.table::data.table(
    PC = pca_columns,
    percent_variance_explained = NA_real_,
    pca_basis = "not_run_no_expression_source"
  )

  coverage <- expression_score_map[, .(gene_set, score_column)]
  coverage[, configured_genes := vapply(gene_set, function(set_name) paste(gene_sets[[set_name]], collapse = ";"), character(1))]
  coverage[, n_configured_genes := vapply(gene_set, function(set_name) length(gene_sets[[set_name]]), integer(1))]
  coverage[, `:=`(found_genes = "", missing_genes = configured_genes, n_found_genes = 0L, n_missing_genes = n_configured_genes)]

  list(
    expression_matrix = matrix(numeric(), nrow = 0, ncol = 0, dimnames = list(character(), character())),
    sample_metadata = sample_universe[0],
    sample_scores = sample_scores,
    project_scores = summarize_project_expression(sample_scores, projects),
    pca_scores = pca,
    pca_variance = variance,
    gene_set_coverage = coverage,
    source = "no_source_available",
    data_type = "not_available",
    normalization_method = "not_applied_no_expression_source",
    raw_expression_sample_count = 0L,
    tumor_expression_column_count = 0L,
    pca_status = "not_run_no_expression_source",
    parsing_warnings = reason
  )
}

summarize_project_expression <- function(sample_scores, projects) {
  score_cols <- c(expression_score_columns, lineage_pc_columns)
  sourced <- sample_scores[source != "no_source_available"]
  if (nrow(sourced) == 0) {
    summary <- copy(projects)
    summary[, n_samples := 0L]
    for (column in score_cols) {
      summary[, (paste0("median_", column)) := NA_real_]
      summary[, (paste0("missingness_", column)) := NA_real_]
    }
    data.table::setorder(summary, project_id)
    return(summary)
  }

  summary <- sourced[, c(
    list(n_samples = .N),
    lapply(.SD, function(x) stats::median(x, na.rm = TRUE)),
    setNames(lapply(.SD, function(x) mean(is.na(x))), paste0("missingness_", score_cols))
  ), by = .(project_id, project_code), .SDcols = score_cols]

  for (column in score_cols) {
    median_col <- paste0("median_", column)
    data.table::setnames(summary, column, median_col)
    summary[is.infinite(get(median_col)) | is.nan(get(median_col)), (median_col) := NA_real_]
  }

  summary <- merge(projects, summary, by = c("project_id", "project_code"), all.x = TRUE, sort = FALSE)
  summary[is.na(n_samples), n_samples := 0L]
  data.table::setorder(summary, project_id)
  summary
}

build_expression_features <- function(projects, sample_universe, purity_samples, aneuploidy_samples, gene_sets, expression_path, source_status) {
  if (is.na(expression_path) || !file.exists(expression_path) || file.info(expression_path)$size == 0) {
    log_warn("No expression source file was available; writing mutation-sample-aligned missing-feature tables")
    return(missing_expression_outputs(projects, sample_universe, gene_sets, "no_expression_source_found"))
  }

  data_type <- detect_expression_data_type(expression_path)
  subset <- read_expression_subset(expression_path, gene_sets, sample_universe)
  normalized <- normalize_expression_matrix(subset$expression_matrix, data_type)

  if (nrow(subset$sample_metadata) == 0) {
    missing <- missing_expression_outputs(projects, sample_universe, gene_sets, "expression_source_had_no_mutation_matched_tumor_samples")
    missing$source <- basename(expression_path)
    missing$data_type <- data_type
    missing$normalization_method <- normalized$method
    missing$raw_expression_sample_count <- subset$raw_expression_sample_count
    missing$tumor_expression_column_count <- subset$tumor_expression_column_count
    missing$parsing_warnings <- subset$parsing_warnings
    return(missing)
  }

  scores <- score_gene_sets(normalized$matrix, subset$sample_metadata, gene_sets)
  pca <- compute_expression_pca(normalized$matrix, subset$sample_metadata)
  scores[, lineage_or_tissue_PC1 := pca$scores$PC1[match(sample_barcode, pca$scores$sample_barcode)]]
  scores[, lineage_or_tissue_PC2 := pca$scores$PC2[match(sample_barcode, pca$scores$sample_barcode)]]
  scores[, source := basename(expression_path)]
  sample_type_note <- if (any(subset$sample_metadata$sample_type_code == "06")) {
    "tumor_samples_include_metastatic_sample_type_06_when_matched_to_mutation_layer"
  } else {
    "tumor_samples_exclude_normal_samples"
  }
  scores[, notes := paste(
    c(
      paste0("source_status=", source_status),
      paste0("data_type=", data_type),
      paste0("normalization=", normalized$method),
      sample_type_note,
      "pathway_scores_are_mean_z_scores_from_configured_gene_sets"
    ),
    collapse = ";"
  )]
  scores <- add_missing_score_columns(scores)
  data.table::setcolorder(scores, c("sample_barcode", "patient_barcode", "project_id", "project_code", expression_score_columns, lineage_pc_columns, "source", "notes"))

  project_scores <- summarize_project_expression(scores, projects)
  coverage <- gene_set_coverage(gene_sets, normalized$matrix)

  list(
    expression_matrix = normalized$matrix,
    sample_metadata = subset$sample_metadata,
    sample_scores = scores,
    project_scores = project_scores,
    pca_scores = pca$scores,
    pca_variance = pca$variance,
    gene_set_coverage = coverage,
    source = basename(expression_path),
    data_type = data_type,
    normalization_method = normalized$method,
    raw_expression_sample_count = subset$raw_expression_sample_count,
    tumor_expression_column_count = subset$tumor_expression_column_count,
    pca_status = pca$status,
    parsing_warnings = subset$parsing_warnings,
    purity_overlap = sum(scores$sample_barcode %in% purity_samples$sample_barcode),
    aneuploidy_overlap = sum(scores$sample_barcode %in% aneuploidy_samples$sample_barcode)
  )
}

write_expression_matrix <- function(expression_matrix, path) {
  if (nrow(expression_matrix) == 0 || ncol(expression_matrix) == 0) {
    write_tsv_gz_table(data.table::data.table(gene_symbol = character()), path)
    return(invisible(path))
  }
  table <- data.table::data.table(gene_symbol = rownames(expression_matrix))
  table <- cbind(table, data.table::as.data.table(expression_matrix))
  write_tsv_gz_table(table, path)
}

write_qc_summary <- function(features, projects, sample_universe, purity_samples, aneuploidy_samples, retrieval_record, output_path) {
  expressed_samples <- features$sample_scores[source != "no_source_available"]
  represented_projects <- sort(unique(expressed_samples$project_id))
  represented_projects <- represented_projects[!is.na(represented_projects)]
  missing_projects <- setdiff(projects$project_id, represented_projects)

  mutation_sample_count <- sample_universe[, uniqueN(sample_barcode)]
  expression_sample_count <- expressed_samples[, uniqueN(sample_barcode)]
  purity_overlap <- expressed_samples[sample_barcode %in% purity_samples$sample_barcode, uniqueN(sample_barcode)]
  aneuploidy_overlap <- expressed_samples[sample_barcode %in% aneuploidy_samples$sample_barcode, uniqueN(sample_barcode)]

  global <- data.table::data.table(
    qc_section = "global",
    project_id = "ALL",
    project_code = "ALL",
    metric = c(
      "expression_source_file_status",
      "expression_source_file",
      "expression_data_source",
      "expression_data_type",
      "normalization_method",
      "genes_available_in_scoring_matrix",
      "raw_expression_sample_columns",
      "tumor_expression_sample_columns_matched_to_mutation_layer",
      "expression_samples_processed",
      "tcga_projects_with_expression_data",
      "expected_tcga_projects_missing_expression_data",
      "missing_tcga_projects_expression_data",
      "mutation_layer_samples",
      "samples_overlapping_mutation_layer",
      "samples_overlapping_purity_ploidy_layer",
      "samples_overlapping_aneuploidy_layer",
      "percent_mutation_layer_samples_with_expression",
      "pca_status",
      "expression_matrix_output",
      "parsing_warnings",
      "limitations"
    ),
    value = c(
      retrieval_record$status[[1]],
      retrieval_record$path[[1]],
      features$source,
      features$data_type,
      features$normalization_method,
      as.character(nrow(features$expression_matrix)),
      as.character(features$raw_expression_sample_count),
      as.character(features$tumor_expression_column_count),
      as.character(expression_sample_count),
      as.character(length(represented_projects)),
      as.character(length(missing_projects)),
      ifelse(length(missing_projects) == 0, "none", paste(missing_projects, collapse = ";")),
      as.character(mutation_sample_count),
      as.character(expression_sample_count),
      as.character(purity_overlap),
      as.character(aneuploidy_overlap),
      sprintf("%.4f", ifelse(mutation_sample_count == 0, NA_real_, expression_sample_count / mutation_sample_count * 100)),
      features$pca_status,
      "data/processed/expression/expression_matrix_by_sample.tsv.gz",
      features$parsing_warnings,
      "PCA is computed from the configured expression-gene subset, not the full transcriptome, unless the local source is pre-filtered to a broader set."
    )
  )

  project_metrics <- c(
    "n_samples",
    paste0("median_", c(expression_score_columns, lineage_pc_columns)),
    paste0("missingness_", c(expression_score_columns, lineage_pc_columns))
  )
  per_project_input <- copy(features$project_scores[, c("project_id", "project_code", project_metrics), with = FALSE])
  for (metric in setdiff(project_metrics, c("project_id", "project_code"))) {
    per_project_input[, (metric) := as.character(get(metric))]
  }
  per_project <- data.table::melt(
    per_project_input,
    id.vars = c("project_id", "project_code"),
    variable.name = "metric",
    value.name = "value"
  )
  per_project[, qc_section := "per_project"]

  gene_set_qc <- data.table::rbindlist(lapply(seq_len(nrow(features$gene_set_coverage)), function(i) {
    row <- features$gene_set_coverage[i]
    data.table::data.table(
      qc_section = "gene_set",
      project_id = "ALL",
      project_code = "ALL",
      metric = paste(row$gene_set, c("n_configured_genes", "n_found_genes", "n_missing_genes", "missing_genes"), sep = ":"),
      value = c(as.character(row$n_configured_genes), as.character(row$n_found_genes), as.character(row$n_missing_genes), row$missing_genes)
    )
  }), use.names = TRUE)

  qc <- data.table::rbindlist(
    list(global, per_project[, .(qc_section, project_id, project_code, metric, value)], gene_set_qc),
    use.names = TRUE,
    fill = TRUE
  )
  data.table::setorder(qc, qc_section, project_id, metric)
  write_tsv_table(qc, output_path)

  list(
    expression_sample_count = expression_sample_count,
    represented_projects = length(represented_projects),
    missing_projects = missing_projects,
    mutation_sample_count = mutation_sample_count,
    purity_overlap = purity_overlap,
    aneuploidy_overlap = aneuploidy_overlap
  )
}

write_qc_figure <- function(project_scores, output_path) {
  ensure_dir(dirname(output_path))
  grDevices::pdf(output_path, width = 10, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(8, 4, 3, 1))
  values <- project_scores$n_samples
  values[is.na(values)] <- 0
  names(values) <- project_scores$project_code
  graphics::barplot(
    values,
    las = 2,
    ylab = "Expression samples",
    main = "Expression layer QC: mutation-matched tumor samples"
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
  require_expression_packages()
  root <- find_project_root()

  projects <- read_projects(project_path("data", "interim", "tcga_projects.tsv", root = root))
  mutation_path <- project_path("data", "processed", "mutations", "mc3_somatic_mutations.parquet", root = root)
  tmb_path <- project_path("data", "processed", "features", "tmb_by_sample.tsv", root = root)
  purity_path <- project_path("data", "processed", "features", "purity_ploidy_by_sample.tsv", root = root)
  aneuploidy_path <- project_path("data", "processed", "features", "aneuploidy_by_sample.tsv", root = root)
  gene_set_path <- project_path("config", "expression_gene_sets.yaml", root = root)

  sample_universe <- read_mutation_sample_universe(mutation_path, tmb_path)
  purity_samples <- read_sample_table_optional(purity_path)
  aneuploidy_samples <- read_sample_table_optional(aneuploidy_path)
  gene_sets <- read_expression_gene_sets(gene_set_path)
  log_info(sprintf("Loaded %d mutation-layer samples as expression sample universe", nrow(sample_universe)))

  retrieval_record <- retrieve_expression_resource(root, refresh = args$refresh)
  expression_path <- discover_local_expression_resource(root, retrieval_record)
  if (is.na(expression_path)) {
    log_warn("No local expression matrix was found under data/raw/expression/")
  } else {
    log_info(sprintf("Using expression matrix: %s", expression_path))
  }

  features <- build_expression_features(
    projects,
    sample_universe,
    purity_samples,
    aneuploidy_samples,
    gene_sets,
    expression_path,
    retrieval_record$status[[1]]
  )

  matrix_output <- project_path("data", "processed", "expression", "expression_matrix_by_sample.tsv.gz", root = root)
  sample_output <- project_path("data", "processed", "features", "expression_pathway_scores_by_sample.tsv", root = root)
  project_output <- project_path("data", "processed", "features", "expression_pathway_scores_by_project.tsv", root = root)
  pca_output <- project_path("data", "processed", "features", "expression_pca_by_sample.tsv", root = root)
  variance_output <- project_path("data", "processed", "features", "expression_pca_variance.tsv", root = root)
  qc_output <- project_path("results", "tables", "expression_layer_qc_summary.tsv", root = root)
  figure_output <- project_path("results", "figures", "expression_layer_qc_overview.pdf", root = root)

  write_expression_matrix(features$expression_matrix, matrix_output)
  write_tsv_table(features$sample_scores, sample_output)
  write_tsv_table(features$project_scores, project_output)
  write_tsv_table(features$pca_scores, pca_output)
  write_tsv_table(features$pca_variance, variance_output)
  qc <- write_qc_summary(features, projects, sample_universe, purity_samples, aneuploidy_samples, retrieval_record, qc_output)
  write_qc_figure(features$project_scores, figure_output)

  validate_outputs_nonempty(c(matrix_output, sample_output, project_output, pca_output, variance_output, qc_output, figure_output))
  log_info(sprintf("Expression layer complete: %d projects represented", qc$represented_projects))
  log_info(sprintf("Expression samples processed: %d of %d mutation-layer samples", qc$expression_sample_count, qc$mutation_sample_count))
  log_info(sprintf("Overlaps: purity/ploidy=%d, aneuploidy=%d", qc$purity_overlap, qc$aneuploidy_overlap))
}

if (identical(environment(), globalenv())) {
  main()
}
