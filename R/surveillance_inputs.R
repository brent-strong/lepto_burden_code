load_surveillance_model_inputs <- function(path = analysis_path("surveillance_model_inputs.csv")) {
  `%>%` <- magrittr::`%>%`
  all_regions <- c("Arecibo", "Bayamón", "Caguas", "Fajardo", "Mayagüez", "Metro", "Ponce")

  inputs <- readr::read_csv(
    path,
    col_types = readr::cols(
      health_region = readr::col_character(),
      site = readr::col_character(),
      input_type = readr::col_character(),
      n_t = readr::col_integer(),
      n_p = readr::col_integer(),
      n_pa = readr::col_integer(),
      passive_surveillance_capture = readr::col_integer(),
      col2_igm = readr::col_integer(),
      col3_pcr = readr::col_integer(),
      count = readr::col_integer()
    )
  )

  counts <- inputs %>%
    dplyr::filter(input_type == "counts")

  pattern_counts <- inputs %>%
    dplyr::filter(input_type == "test_pattern")

  all_sites <- counts$site %>%
    unique() %>%
    sort() %>%
    stats::na.omit() %>%
    as.character()

  make_count_matrix <- function(value_col) {
    counts %>%
      dplyr::select(health_region, site, count = dplyr::all_of(value_col)) %>%
      tidyr::complete(health_region = all_regions, site = all_sites, fill = list(count = 0L)) %>%
      tidyr::pivot_wider(names_from = site, values_from = count) %>%
      dplyr::arrange(factor(health_region, levels = all_regions)) %>%
      tibble::column_to_rownames("health_region") %>%
      as.matrix()
  }

  n_t <- make_count_matrix("n_t")
  n_p <- make_count_matrix("n_p")
  n_pa <- make_count_matrix("n_pa")

  max_rows <- max(n_t)
  all_blocks <- list()

  for (hosp in all_sites) {
    for (reg in all_regions) {
      block_data <- pattern_counts %>%
        dplyr::filter(site == hosp, health_region == reg) %>%
        tidyr::uncount(count) %>%
        dplyr::select(passive_surveillance_capture, col2_igm, col3_pcr)

      actual_row_count <- nrow(block_data)

      if (actual_row_count < max_rows) {
        padding_needed <- max_rows - actual_row_count
        padding <- matrix(-999, nrow = padding_needed, ncol = 3)
        colnames(padding) <- colnames(block_data)
        combined_block <- rbind(as.matrix(block_data), padding)
      } else {
        combined_block <- as.matrix(block_data[1:max_rows, ])
      }

      all_blocks[[paste(hosp, reg, sep = "_")]] <- combined_block
    }
  }

  Y <- array(
    unlist(all_blocks),
    dim = c(max_rows, 3, length(all_blocks))
  )

  list(n_t = n_t, n_p = n_p, n_pa = n_pa, Y = Y)
}
