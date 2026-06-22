#!/usr/bin/env Rscript

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

rename_map <- c(
  A = "beta_0",
  B = "beta_1",
  alpha = "alpha_0",
  n_r = "n_p",
  n_rt = "n_pa",
  pi_DR = "pi_P",
  piDR = "pi_P",
  pi_DT = "pi_A",
  piDT = "pi_A"
)

rename_exact <- function(x) {
  if (is.null(x)) return(x)
  idx <- x %in% names(rename_map)
  x[idx] <- unname(rename_map[x[idx]])
  x <- sub("^alpha\\[", "alpha_0[", x)
  x <- sub("^n_r\\[", "n_p[", x)
  x <- sub("^n_rt\\[", "n_pa[", x)
  x <- sub("^pi_DR\\[", "pi_P[", x)
  x <- sub("^piDR\\[", "pi_P[", x)
  x <- sub("^pi_DT\\[", "pi_A[", x)
  x <- sub("^piDT\\[", "pi_A[", x)
  x
}

rename_parameter_names <- function(x) {
  changed <- FALSE

  rename_names <- function(obj) {
    before <- names(obj)
    after <- rename_exact(before)
    if (!identical(before, after)) {
      names(obj) <- after
      attr(obj, "changed_parameter_names") <- TRUE
    }
    obj
  }

  rename_cols <- function(obj) {
    before <- colnames(obj)
    after <- rename_exact(before)
    if (!identical(before, after)) {
      colnames(obj) <- after
      attr(obj, "changed_parameter_names") <- TRUE
    }
    obj
  }

  rename_dimnames <- function(obj) {
    before <- dimnames(obj)
    if (is.null(before)) return(obj)
    after <- lapply(before, rename_exact)
    if (!identical(before, after)) {
      dimnames(obj) <- after
      attr(obj, "changed_parameter_names") <- TRUE
    }
    obj
  }

  updated_names <- rename_names(x)
  if (isTRUE(attr(updated_names, "changed_parameter_names"))) changed <- TRUE
  attr(updated_names, "changed_parameter_names") <- NULL
  x <- updated_names

  if (inherits(x, "mcmc.list")) {
    for (i in seq_along(x)) {
      updated <- rename_parameter_names(x[[i]])
      if (isTRUE(attr(updated, "changed_parameter_names"))) changed <- TRUE
      attr(updated, "changed_parameter_names") <- NULL
      x[[i]] <- updated
    }
  } else if (inherits(x, "mcmc") || is.matrix(x) || is.data.frame(x)) {
    updated <- rename_cols(x)
    if (isTRUE(attr(updated, "changed_parameter_names"))) changed <- TRUE
    attr(updated, "changed_parameter_names") <- NULL
    x <- updated
  } else if (is.array(x)) {
    updated <- rename_dimnames(x)
    if (isTRUE(attr(updated, "changed_parameter_names"))) changed <- TRUE
    attr(updated, "changed_parameter_names") <- NULL
    x <- updated
  } else if (is.list(x)) {
    for (i in seq_along(x)) {
      updated <- rename_parameter_names(x[[i]])
      if (isTRUE(attr(updated, "changed_parameter_names"))) changed <- TRUE
      attr(updated, "changed_parameter_names") <- NULL
      x[[i]] <- updated
    }
  }

  attr(x, "changed_parameter_names") <- changed
  x
}

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

files <- list.files(
  LEPTO_ROOT,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

changed_files <- character()
failed_files <- character()

for (file in files) {
  obj <- tryCatch(readRDS(file), error = function(e) e)
  if (inherits(obj, "error")) {
    failed_files <- c(failed_files, file)
    next
  }

  updated <- rename_parameter_names(obj)
  changed <- isTRUE(attr(updated, "changed_parameter_names"))
  attr(updated, "changed_parameter_names") <- NULL

  if (changed) {
    changed_files <- c(changed_files, file)
    if (!dry_run) {
      saveRDS(updated, file)
    }
  }
}

cat(if (dry_run) "Dry run" else "Updated", length(changed_files), "of", length(files), "RDS files\n")
if (length(changed_files) > 0) {
  cat(paste0(" - ", normalizePath(changed_files, mustWork = FALSE), collapse = "\n"), "\n")
}
if (length(failed_files) > 0) {
  cat("Failed to read", length(failed_files), "RDS files\n")
  cat(paste0(" - ", normalizePath(failed_files, mustWork = FALSE), collapse = "\n"), "\n")
}
