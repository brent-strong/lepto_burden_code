# base hospitals always included when n_h >= 0
base_hospitals <- c("hosp_2512", "hosp_14","hosp_15", "hosp_23")

# function to get hospitals given n_h
get_hospitals <- function(data, hosp_region, n_h = 0) {
  
  # merge to get regions
  hosp_info <- data %>%
    dplyr::select(hospital) %>%
    distinct() %>%
    left_join(hosp_region, by = "hospital")
  
  # get region membership of base hospitals
  base_regions <- hosp_info %>%
    filter(hospital %in% base_hospitals) %>%
    pull(hosp_region) %>%
    unique()
  
  other_regions <- setdiff(unique(hosp_info$hosp_region), base_regions)
  
  # start with base hospitals
  selected <- base_hospitals
  
  if (n_h > 0) {
    # Case 1: regions without base hospitals → sample n_h hospitals each
    # Case 1: regions without base hospitals → sample n_h hospitals each
    add_other <- other_regions %>%
      map(~ {
        hosp_info %>%
          filter(hosp_region == .x, !(hospital %in% base_hospitals)) %>%
          slice_sample(n = min(n_h, nrow(.)))
      }) %>%
      bind_rows()
    
    # Case 2: regions with base hospitals → sample (n_h - 1) hospitals each
    add_base <- NULL
    if (n_h >= 2) {
      add_base <- base_regions %>%
        map(~ {
          hosp_info %>%
            filter(hosp_region == .x, !(hospital %in% base_hospitals)) %>%
            slice_sample(n = min(n_h - 1, nrow(.)))
        }) %>%
        bind_rows()
    }
    
    selected <- c(selected, add_other$hospital, add_base$hospital)
  }
  
  # filter original data
  data %>%
    left_join(hosp_region, by = "hospital") %>%
    filter(hospital %in% selected)
}
