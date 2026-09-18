timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

log_info <- function(message) {
  cat(sprintf("[%s] INFO  %s\n", timestamp(), message))
}

log_warn <- function(message) {
  cat(sprintf("[%s] WARN  %s\n", timestamp(), message))
}

abort <- function(message) {
  stop(sprintf("[%s] ERROR %s", timestamp(), message), call. = FALSE)
}

this_script <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
}

find_project_root <- function(start = NULL) {
  current <- if (is.null(start)) dirname(this_script()) else start
  current <- normalizePath(current, mustWork = FALSE)

  repeat {
    marker <- file.path(current, "config", "tcga_projects.yaml")
    if (file.exists(marker)) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      abort("Could not locate project root containing config/tcga_projects.yaml")
    }
    current <- parent
  }
}

project_path <- function(..., root = find_project_root()) {
  file.path(root, ...)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

require_file <- function(path, description = "required file") {
  if (!file.exists(path)) {
    abort(sprintf("Missing %s: %s", description, path))
  }
  invisible(path)
}

read_yaml_config <- function(path) {
  require_file(path, "YAML configuration")
  if (!requireNamespace("yaml", quietly = TRUE)) {
    abort("Package 'yaml' is required")
  }
  yaml::read_yaml(path)
}

read_csv_config <- function(path) {
  require_file(path, "CSV configuration")
  if (!requireNamespace("readr", quietly = TRUE)) {
    abort("Package 'readr' is required")
  }
  readr::read_csv(path, show_col_types = FALSE)
}

write_tsv <- function(data, path) {
  ensure_dir(dirname(path))
  if (!requireNamespace("readr", quietly = TRUE)) {
    abort("Package 'readr' is required")
  }
  readr::write_tsv(data, path)
  log_info(sprintf("Wrote %s", path))
  invisible(path)
}

load_tcga_project_config <- function(root = find_project_root()) {
  cfg <- read_yaml_config(file.path(root, "config", "tcga_projects.yaml"))
  projects <- cfg$projects
  if (length(projects) == 0) {
    abort("No TCGA projects configured")
  }
  data.frame(
    code = vapply(projects, `[[`, character(1), "code"),
    project_id = vapply(projects, `[[`, character(1), "project_id"),
    name = vapply(projects, `[[`, character(1), "name"),
    stringsAsFactors = FALSE
  )
}

not_implemented <- function(step, next_action) {
  abort(sprintf(
    "%s is a starter placeholder. Next implementation step: %s",
    step,
    next_action
  ))
}
