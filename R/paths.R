# Repository-relative path helpers.
#
# Run scripts from the repository root, or set LEPTO_BURDEN_ROOT when launching
# from another working directory, for example on an HPC cluster.

find_repo_root <- function(start = getwd()) {
  env_root <- Sys.getenv("LEPTO_BURDEN_ROOT", unset = NA_character_)
  if (!is.na(env_root) && nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = TRUE))
  }

  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "README.md")) &&
        dir.exists(file.path(path, "Data")) &&
        dir.exists(file.path(path, "Simulation"))) {
      return(path)
    }

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find the lepto_burden repository root. Run from the repo root or set LEPTO_BURDEN_ROOT.")
    }
    path <- parent
  }
}

LEPTO_ROOT <- find_repo_root()

repo_root <- function() {
  LEPTO_ROOT
}

repo_path <- function(...) {
  file.path(LEPTO_ROOT, ...)
}

data_path <- function(...) {
  repo_path("Data", ...)
}

analysis_path <- function(...) {
  repo_path("Data_Analysis", ...)
}

simulation_path <- function(...) {
  repo_path("Simulation", ...)
}
