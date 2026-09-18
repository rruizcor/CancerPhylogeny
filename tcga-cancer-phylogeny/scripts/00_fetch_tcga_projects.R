#!/usr/bin/env Rscript

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "lib", "common.R"))

gdc_projects_url <- function() {
  paste0(
    "https://api.gdc.cancer.gov/projects?",
    paste(
      c(
        "format=JSON",
        "size=1000",
        "sort=project_id:asc",
        "fields=project_id,name,disease_type,primary_site,dbgap_accession_number,released,state,program.name,summary.case_count,summary.file_count"
      ),
      collapse = "&"
    )
  )
}

download_gdc_projects_json <- function(cache_path, refresh = FALSE) {
  ensure_dir(dirname(cache_path))

  if (file.exists(cache_path) && !refresh) {
    log_info(sprintf("Using cached GDC project metadata: %s", cache_path))
    return(cache_path)
  }

  url <- gdc_projects_url()
  tmp <- tempfile(pattern = "gdc_projects_", fileext = ".json")
  log_info(sprintf("Fetching GDC project metadata from %s", url))

  tryCatch(
    utils::download.file(url, tmp, mode = "wb", quiet = TRUE),
    error = function(err) {
      abort(sprintf("Failed to download GDC project metadata: %s", conditionMessage(err)))
    }
  )

  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    abort("Downloaded GDC project metadata file is empty")
  }

  file.copy(tmp, cache_path, overwrite = TRUE)
  unlink(tmp)
  log_info(sprintf("Cached GDC project metadata: %s", cache_path))
  cache_path
}

read_gdc_project_hits <- function(json_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    abort("Package 'jsonlite' is required")
  }

  require_file(json_path, "cached GDC project metadata JSON")
  payload <- jsonlite::fromJSON(json_path, flatten = TRUE)

  if (is.null(payload$data$hits) || nrow(payload$data$hits) == 0) {
    abort("GDC project metadata response did not contain project hits")
  }

  payload$data$hits
}

standardize_gdc_projects <- function(hits) {
  collapse_values <- function(value) {
    if (is.null(value) || length(value) == 0) {
      return(NA_character_)
    }
    value <- as.character(value)
    value <- value[!is.na(value) & value != ""]
    if (length(value) == 0) {
      return(NA_character_)
    }
    paste(unique(value), collapse = ";")
  }

  optional_column <- function(name, default = NA) {
    if (!name %in% names(hits)) {
      return(rep(default, nrow(hits)))
    }

    value <- hits[[name]]
    if (is.list(value)) {
      return(vapply(value, collapse_values, character(1)))
    }

    value
  }

  data.frame(
    project_id = optional_column("project_id"),
    disease_type = optional_column("disease_type"),
    primary_site = optional_column("primary_site"),
    name = optional_column("name"),
    program = optional_column("program.name"),
    dbgap_accession_number = optional_column("dbgap_accession_number"),
    released = optional_column("released"),
    state = optional_column("state"),
    case_count = optional_column("summary.case_count"),
    file_count = optional_column("summary.file_count"),
    stringsAsFactors = FALSE
  )
}

fetch_gdc_project_metadata <- function(project_ids, cache_path, refresh = FALSE) {
  json_path <- download_gdc_projects_json(cache_path, refresh = refresh)
  hits <- read_gdc_project_hits(json_path)
  projects <- standardize_gdc_projects(hits)
  projects <- projects[projects$project_id %in% project_ids, , drop = FALSE]
  projects <- projects[match(project_ids, projects$project_id), , drop = FALSE]

  missing <- project_ids[!project_ids %in% projects$project_id]
  if (length(missing) > 0) {
    abort(sprintf("GDC response is missing configured projects: %s", paste(missing, collapse = ", ")))
  }

  if (anyDuplicated(projects$project_id) > 0) {
    abort("GDC response contains duplicate configured project IDs")
  }

  projects$project_code <- sub("^TCGA-", "", projects$project_id)
  projects <- projects[, c(
    "project_id",
    "project_code",
    "disease_type",
    "primary_site",
    "name",
    "program",
    "dbgap_accession_number",
    "released",
    "state",
    "case_count",
    "file_count"
  )]

  rownames(projects) <- NULL
  projects
}

write_project_metadata <- function(projects, output_path) {
  write_tsv(projects, output_path)
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  list(
    refresh = "--refresh" %in% args
  )
}

main <- function() {
  args <- parse_args()
  root <- find_project_root()
  configured <- load_tcga_project_config(root)
  log_info(sprintf("Loaded %d configured TCGA projects", nrow(configured)))

  cache_path <- project_path("data", "raw", "gdc", "projects.json", root = root)
  output_path <- project_path("data", "interim", "tcga_projects.tsv", root = root)
  log_info(sprintf("Raw metadata cache: %s", cache_path))
  log_info(sprintf("Planned output: %s", output_path))

  projects <- fetch_gdc_project_metadata(
    configured$project_id,
    cache_path = cache_path,
    refresh = args$refresh
  )
  write_project_metadata(projects, output_path)
  log_info(sprintf("Retrieved expected TCGA project count: %d", nrow(projects)))
}

if (identical(environment(), globalenv())) {
  main()
}
