data_generation <- function(lambda, pi, uafi_lambda,
                            sensitivity, specificity, test_receipt_probs,
                            cell_probabilities = NULL) {
  library(tidyverse)

  if (is.null(cell_probabilities)) {
    if (!exists("data_path")) {
      source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))
    }
    sim_df <- read_csv(
      data_path("processed", "simulation_inputs", "cell_probabilities.csv"),
      show_col_types = FALSE
    ) %>%
      as.data.frame()
  } else {
    sim_df <- as.data.frame(cell_probabilities)
  }
  
  # STEP 1: Negative binomial draw for true disease cases
  sim_df$n_d <- rpois(n = nrow(sim_df),
                        lambda = sim_df$population/100000 * lambda)
  
  hospital_matrix <- as.matrix(sim_df[, 2:65])
  hosp_vars <- names(sim_df) %>% keep(~ startsWith(.x, "hosp_"))
  
  sim_df <- sim_df %>%
    mutate(across(all_of(hosp_vars), .fns = list(cases = ~ .), .names = "{.col}_cases")) 
  
  for (i in 1:nrow(sim_df)) {
    sim_df[i, 69:132] <- as.vector(rmultinom(n = 1, size = sim_df$n_d[i], prob = hospital_matrix[i, ]))
  }
  
  hospital_summary <- sim_df %>%
    dplyr::select(health_region, hosp_832_cases:hosp_678_cases) %>%
    pivot_longer(cols = ends_with("_cases"), names_to = "hospital", values_to = "count") %>%
    group_by(health_region, hospital) %>%
    summarise(n_d = sum(count), .groups = "drop")
  
  # Reported cases
  hospital_summary$n_r <- rbinom(n = nrow(hospital_summary),
                                 size = hospital_summary$n_d,
                                 prob = pi)
  
  region_summary <- hospital_summary %>%
    group_by(health_region) %>%
    summarise(y = sum(n_r))
  
  # Negative binomial for total fever cases
  sim_df$n_t <- rpois(n = nrow(sim_df),
                        lambda = sim_df$population/100000 * uafi_lambda)
  
  sim_n_a_df <- sim_df %>%
    mutate(across(all_of(hosp_vars), .fns = list(cases = ~ .), .names = "{.col}_cases")) %>%
    dplyr::select(-n_d)
  
  for (i in 1:nrow(sim_df)) {
    sim_n_a_df[i, 68:131] <- as.vector(rmultinom(n = 1, size = sim_n_a_df$n_t[i], prob = hospital_matrix[i, ]))
  }
  
  hospital_cases <- sim_n_a_df %>%
    dplyr::select(health_region, hosp_832_cases:hosp_678_cases) %>%
    pivot_longer(cols = ends_with("_cases"), names_to = "hospital", values_to = "count") %>%
    group_by(health_region, hospital) %>%
    summarise(n_t_prelim = sum(count), .groups = "drop") %>%
    right_join(hospital_summary, by = c("health_region", "hospital")) %>%
    mutate(n_t = n_t_prelim + n_d) %>%
    dplyr::select(-n_t_prelim)
  
  # Create patient-level data
  t <- length(sensitivity)
  
  patient_list <- hospital_cases %>%
    pmap(function(health_region, hospital, n_d, n_r, n_t) {
      
      num_tests <- length(sensitivity)   # number of different tests
      if (n_t == 0) {
        return(tibble(
          health_region = character(0),
          hospital = character(0),
          disease_status = integer(0),
          passive_surveillance_capture = integer(0),
          !!!set_names(vector("list", num_tests), paste0("test", 1:num_tests))
        ))
      }
      
      # disease status and passive capture (unchanged)
      disease_status <- c(rep(1, n_d), rep(0, n_t - n_d))
      disease_status <- sample(disease_status)
      
      passive_capture <- integer(n_t)
      pos_indices <- which(disease_status == 1)
      neg_indices <- which(disease_status == 0)
      
      passive_capture[pos_indices] <- 0
      if (n_r > 0) {
        sampled_pos <- sample(pos_indices, size = min(n_r, length(pos_indices)))
        passive_capture[sampled_pos] <- 1
      }
      passive_capture[neg_indices] <- -1
      
      df_patients <- tibble(
        health_region = health_region,
        hospital = hospital,
        disease_status = disease_status,
        passive_surveillance_capture = passive_capture
      )
      
      # --- Test 1: generate receipt and results ---
      prob1 <- test_receipt_probs[1]
      sens1 <- sensitivity[1]
      spec1 <- specificity[1]
      
      received1 <- runif(n_t) < prob1        # TRUE if patient received test1
      test1 <- integer(n_t)
      test1[!received1] <- -1               # not received => -1
      
      # For those who received and are positive:
      idx_pos_recv <- which(received1 & disease_status == 1)
      if (length(idx_pos_recv) > 0) {
        test1[idx_pos_recv] <- rbinom(length(idx_pos_recv), 1, sens1)
      }
      # For those who received and are negative (false positive prob = 1 - spec)
      idx_neg_recv <- which(received1 & disease_status == 0)
      if (length(idx_neg_recv) > 0) {
        test1[idx_neg_recv] <- rbinom(length(idx_neg_recv), 1, 1 - spec1)
      }
      
      tests <- vector("list", num_tests)
      tests[[1]] <- test1
      
      # --- Tests 2..num_tests: conditional on having received test1 ---
      if (num_tests > 1) {
        for (i in 2:num_tests) {
          prob_i <- test_receipt_probs[i]
          sens_i <- sensitivity[i]
          spec_i <- specificity[i]
          
          # only patients who received test1 can possibly receive test i
          # then independently they receive i with probability prob_i
          received_i <- received1 & (runif(n_t) < prob_i)
          
          test_i <- integer(n_t)
          test_i[!received_i] <- -1
          
          idx_pos_i <- which(received_i & disease_status == 1)
          if (length(idx_pos_i) > 0) {
            test_i[idx_pos_i] <- rbinom(length(idx_pos_i), 1, sens_i)
          }
          idx_neg_i <- which(received_i & disease_status == 0)
          if (length(idx_neg_i) > 0) {
            test_i[idx_neg_i] <- rbinom(length(idx_neg_i), 1, 1 - spec_i)
          }
          
          tests[[i]] <- test_i
        }
      }
      
      test_df <- as_tibble(set_names(tests, paste0("test", 1:num_tests)))
      bind_cols(df_patients, test_df)
    })
  
  
  data <- bind_rows(patient_list)  %>%
    mutate(hospital = str_remove(hospital, "_cases"))
  
  return(list(hospital_data = data, region_level_data = region_summary,patient_list=patient_list))
}


