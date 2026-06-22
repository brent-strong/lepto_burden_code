library(stringr)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

file_list <- list.files(
  simulation_path("Simulation_sensitivity"),
  pattern = "\\.rds$",
  full.names = TRUE
)

read_sensitivity_simulation_file <- function(file_list, job, n_h, sd_sens_spec) {
  pattern <- paste0(
    "model_output_nh_", n_h,
    "_SD_sens_spec_", sd_sens_spec,
    "_job_", job,
    "\\.rds$"
  )

  matched_file <- file_list[str_detect(file_list, pattern)]

  if (length(matched_file) == 0) {
    stop("No file found for the given n_h, sensitivity/specificity prior SD, and job number.")
  } else if (length(matched_file) > 1) {
    warning("Multiple files found, using the first match.")
    matched_file <- matched_file[1]
  }

  readRDS(matched_file)
}

# Example:
# model_output <- read_sensitivity_simulation_file(file_list, job = 1, n_h = 0, sd_sens_spec = 0.11)
# Y <- model_output[[1]]$data_list$Y
