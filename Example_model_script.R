# Example capture-recapture model script.
#
# This script generates synthetic hospital- and region-level surveillance data
# from the processed Puerto Rico inputs, builds the marginalized diagnostic
# testing likelihood, fits the NIMBLE capture-recapture model, and produces
# basic MCMC diagnostics. It is not the main real-data analysis script; use
# Data_Analysis/Analysis_code.R for the paper's Puerto Rico data application.

library(tidyverse)
library(nimble)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

# Simulation code

source(data_path("scripts", "Simulation_Code.R"))

#Get hospital presentation probabilities

pres_probs_region<-read.csv(data_path("processed", "simulation_inputs", "health_regions_probs.csv"))[,-1]

#Get population in regions

region_level_data<-read.csv(data_path("processed", "simulation_inputs", "region_population.csv"))
region_level_data$lepto_cases_2023 <- c(42,36,80,4,47,26,47)

#Get data for each 1 km x 1 km cells

cell_data<-read.csv(data_path("processed", "simulation_inputs", "cell_probabilities.csv"))

#Get rainfall data

rainfall<-read.csv(data_path("processed", "simulation_inputs", "rainfall_annual_by_region.csv"))
rain_matrix <- rainfall %>%
  dplyr::select(region, year, total_rain_in) %>%
  tidyr::pivot_wider(names_from = year, values_from = total_rain_in) %>%
  (\(df) {
    region <- df[[1]]                 # keep region separate
    num_mat <- as.matrix(df[,-1])     # numeric part
    num_mat <- scale(num_mat, center = TRUE, scale = FALSE) / 5
    rownames(num_mat) <- region       # optionally set row names
    num_mat                           # return numeric matrix only
  })()


#Get CAR simulation function

source(simulation_path("Simulation_functions", "CAR_counts_sim.R"))

#Get Gaussian process simulation

source(simulation_path("Simulation_functions", "GP_hosp_prob_sim.R"))

#Get hospital sampling code

source(simulation_path("Simulation_functions", "Hospital_sampling.R"))

#Create rates per 100,000

region_level_data <- region_level_data %>% 
  mutate(lepto_per_100k_2023 = lepto_cases_2023*100000/population)
  
mean_log_rate <- mean(log(1/0.33*region_level_data$lepto_per_100k_2023))

#Import in health region shape file and Puerto Rico shapefile.

hospitals <- st_read(data_path("raw", "geospatial", "hospitals", "Hospitals_PR.shp")) %>%
  filter(NAICS_DESC=="GENERAL MEDICAL AND SURGICAL HOSPITALS") %>% dplyr::select(ID,geometry)

health_regions<-st_read(data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp"))

# Spatial join: assign each hospital to a region
hosp_with_region <- st_join(hospitals, health_regions, join = st_within)

# Build non-spatial dataframe
hosp_region <- hosp_with_region %>%
  transmute(
    hospital = paste0("hosp_", ID),  # prefix "hosp_"
    region   # keep region variable
  ) %>%
  rename(hosp_region=region)

#Create logit function

logit <- function(p) {
  if (any(p <= 0 | p >= 1, na.rm = TRUE)) warning("Input probabilities must be in (0, 1)")
  log(p / (1 - p))
}

#Create inverse logit function

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}


data_gen <- function(beta_0,beta_1,r,sigma_Phi,mean_prob,sigma_phi,zeta,
                     multiplier_lower,multiplier_upper,sensitivity,
                     specificity,test_receipt_probs,rain_matrix,
                     n_h,n_years,max_n_d=200){

#Get lambda and pi
  
car_sim_region<-car_sim(sigma_Phi)
Phi<-as.vector(car_sim_region$y)

hosp_probs_sim<-simulate_hospital_probs(hospitals=hosp_region, sigma=sigma_phi, zeta=zeta, mean_prob=mean_prob)
pi <- hosp_probs_sim$probs
pi_df <- data.frame(hospital=colnames(pres_probs_region),pi=pi)

link_df <- tidyr::crossing(region=region_level_data$region,hospital=pi_df$hospital) %>% 
  left_join(pi_df)

#Simulate data

data_list <- list()
latest_year <- 2025
latest_rain_col <- 5  # column for 2025 in rain_matrix

years <- (latest_year - n_years + 1):latest_year

rain_index_indices<-rep(0,length(years))

log_lambda_matrix <- matrix(0, nrow = 7, ncol = length(years))

for (i in seq_along(years)) {
  year <- years[i]
  rain_index <- latest_rain_col - (latest_year - year)
  rain_index_indices[i]<-rain_index# calculates correct column
  
  lambda <- exp(beta_0 + beta_1 * as.numeric(rain_matrix[, rain_index]) + Phi)
  log_lambda_matrix[,i]<-log(lambda)
  
  lambda_df <- data.frame(health_region = region_level_data$region,
                          lambda = lambda)
  
  lambda_full <- left_join(cell_data, lambda_df, by = "health_region")$lambda
  
  data_list[[as.character(year)]] <- data_generation(
    lambda = lambda_full,
    pi = link_df$pi,
    uafi_lambda = 600,
    sensitivity = sensitivity,
    specificity = specificity,
    test_receipt_probs = test_receipt_probs,
    cell_probabilities = cell_data
  )
}

#Get numbers

# Step 0: pre-select hospitals for years < 2025
selected_hospitals <- get_hospitals(data_list[["2024"]]$hospital_data, hosp_region, n_h = n_h)$hospital

#Get IDS of hospitals

surveilled_indices <- which(hosp_region$hospital %in% unique(selected_hospitals))
not_surveilled_indices <- which(!(hosp_region$hospital %in% unique(selected_hospitals)))

hospital_levels <- hosp_region$hospital[surveilled_indices] 

# Prepare empty lists to store matrices per year
n_p_list <- list()
n_pa_list <- list()
n_t_list <- list()

years <- (latest_year - n_years + 1):(latest_year-1)

y <- matrix(
  0,                                # fill with 0
  nrow = nrow(health_regions),      # number of rows
  ncol = length(years)              # number of columns
)

for (year in years) {
  
  # Get hospital_data for this year
  hospital_data <- data_list[[as.character(year)]]$hospital_data
  hospital_data <- hospital_data %>% filter(hospital %in% unique(selected_hospitals)) 
  
  # Process surveillance data
  surveillance <- hospital_data %>%
    left_join(hosp_region, by = "hospital") %>%
    mutate(hospital = factor(hospital, levels = hospital_levels)) %>%
    arrange(hospital) %>%
    mutate(passive_surveillance_capture_mod = ifelse(passive_surveillance_capture != -1,
                                                     passive_surveillance_capture, 0)) %>%
    mutate(tested_true_passive = ifelse((test1 == -1 & test2 == -1) | passive_surveillance_capture_mod == 0, 0, 1)) %>%
    mutate(tested = ifelse(test1 == -1 & test2 == -1, 0, 1)) %>%
    mutate(tested_true = ifelse((test1 == -1 & test2 == -1) | disease_status == 0, 0, 1))
  
  
  #Subtract from region level case totals
  
  passive_cases_surveilled_hospitals<-surveillance %>%
    group_by(health_region) %>%
    summarize(y=sum(passive_surveillance_capture_mod))
  
  region_cases <- data_list[[as.character(year)]]$region_level_data
  y[,which(years==year)] <- region_cases$y - passive_cases_surveilled_hospitals$y
  
  # Ensure all combinations of region × hospital are present
  all_combinations <- expand.grid(
    health_region = unique(hosp_region$hosp_region),
    hospital = hospital_levels
  )
  
  # Summarize n_p
  n_p_mat <- surveillance %>%
    group_by(health_region, hospital) %>%
    summarize(n_p = sum(passive_surveillance_capture_mod), .groups = "drop") %>%
    right_join(all_combinations, by = c("health_region", "hospital")) %>%
    mutate(n_p = replace_na(n_p, 0)) %>%
    pivot_wider(names_from = hospital, values_from = n_p, values_fill = 0) %>%
    dplyr::select(all_of(hospital_levels))
  
  # Summarize n_pa
  n_pa_mat <- surveillance %>%
    group_by(health_region, hospital) %>%
    summarize(n_pa = sum(tested_true_passive), .groups = "drop") %>%
    right_join(all_combinations, by = c("health_region", "hospital")) %>%
    mutate(n_pa = replace_na(n_pa, 0)) %>%
    pivot_wider(names_from = hospital, values_from = n_pa, values_fill = 0) %>%
    dplyr::select(all_of(hospital_levels))
  
  # Summarize n_t
  n_t_mat <- surveillance %>%
    group_by(health_region, hospital) %>%
    summarize(n_t = sum(tested), .groups = "drop") %>%
    right_join(all_combinations, by = c("health_region", "hospital")) %>%
    mutate(n_t = replace_na(n_t, 0)) %>%
    pivot_wider(names_from = hospital, values_from = n_t, values_fill = 0) %>%
    dplyr::select(all_of(hospital_levels))
  
  # Store in lists
  n_p_list[[as.character(year)]] <- n_p_mat
  n_pa_list[[as.character(year)]] <- n_pa_mat
  n_t_list[[as.character(year)]] <- n_t_mat
}

#Create arrays
n_p_array <- array(unlist(n_p_list), dim = c(nrow(n_p_list[["2024"]]), ncol(n_p_list[["2024"]]), length(n_p_list)))
n_pa_array <- array(unlist(n_pa_list), dim = c(nrow(n_pa_list[["2024"]]), ncol(n_pa_list[["2024"]]), length(n_pa_list)))
n_t_array <- array(unlist(n_t_list), dim = c(nrow(n_t_list[["2024"]]), ncol(n_t_list[["2024"]]), length(n_t_list)))


#Make array of test data 

# Find maximum n_t across all bins
max_n_t <- max(n_t_array)

# Function to reorder and pad matrices with -999 rows
pad_matrix <- function(mat, max_n_t) {
  nrow_current <- nrow(mat)
  ncol_current <- ncol(mat)
  
  if (is.null(nrow_current) || nrow_current == 0) {
    # Empty case: full filler
    return(matrix(-999, nrow = max_n_t, ncol = ncol_current))
  } else {
    # Reorder so rows with first column == 1 come first
    mat <- mat[order(-mat[,1]), , drop = FALSE]
    
    if (nrow_current < max_n_t) {
      filler <- matrix(-999, nrow = max_n_t - nrow_current, ncol = ncol_current)
      return(rbind(mat, filler))
    } else {
      return(mat[1:max_n_t, , drop = FALSE])  # truncate if longer
    }
  }
}


# Collect padded matrices into one big list
all_mats <- list()
regions <- sort(unique(hosp_region$hosp_region))       # 1:S
hospitals <- hospital_levels                # 1:H_t
years <- (latest_year - n_years + 1):(latest_year-1)   # 1:T

all_mats <- vector("list", length = length(regions) * length(hospitals) * length(years))
k <- 1
for (t in seq_along(years)) {
  year <- years[t]
  hospital_data <- data_list[[as.character(year)]]$hospital_data
  hospital_data <- hospital_data %>% filter(hospital %in% hospitals)
  surveillance <- hospital_data %>%
    left_join(hosp_region, by = "hospital") %>%
    mutate(passive_surveillance_capture_mod = ifelse(passive_surveillance_capture != -1,
                                                     passive_surveillance_capture, 0)) %>%
    mutate(tested_true_passive = ifelse((test1 == -1 & test2 == -1) | passive_surveillance_capture_mod == 0, 0, 1)) %>%
    mutate(tested = ifelse(test1 == -1 & test2 == -1, 0, 1)) %>%
    mutate(tested_true = ifelse((test1 == -1 & test2 == -1) | disease_status == 0, 0, 1))
  
  for (h in seq_along(hospitals)) {
    hosp <- hospitals[h]
    for (s in seq_along(regions)) {
      reg <- regions[s]
      df_sub <- surveillance %>%
        filter(hospital == hosp & health_region == reg & tested == 1) %>%
        dplyr::select(passive_surveillance_capture_mod, test1, test2)
      
      mat <- if (nrow(df_sub) == 0) matrix(-999, nrow = max_n_t, ncol = 3) else pad_matrix(as.matrix(df_sub), max_n_t)
      
      all_mats[[k]] <- mat
      k <- k + 1
    }
  }
}

# Finally: stack into a 3D array
test_results_array <- array(
  unlist(all_mats),
  dim = c(
    max_n_t,
    3,  # passive_surveillance_capture_mod, test1, test2
    length(all_mats)  # regions × hospitals × years
  )
)

#Remove columns for processing

pres_probs_region_non_surveilled_hospitals <- pres_probs_region[,not_surveilled_indices]

pres_probs_region_surveilled_hospitals <- pres_probs_region[,surveilled_indices]

#Get 2025 data

true_cases_2025<-data_list[["2025"]]$hospital_data %>% 
  group_by(health_region) %>% 
  summarize(d=sum(disease_status))

data_list <- list(
    y = y,
    Y=test_results_array,
    pop = region_level_data$population / 100000
)

constant_list <- list(S=nrow(region_cases),
                      T=length(years),
                      H_t=length(unique(selected_hospitals)),
                      Not_H_t=ncol(as.matrix(pres_probs_region_non_surveilled_hospitals)),
                      H=length(unique(selected_hospitals))+ncol(as.matrix(pres_probs_region_non_surveilled_hospitals)),
                      prob_pres_ns=as.matrix(pres_probs_region_non_surveilled_hospitals),
                      prob_pres_s=as.matrix(pres_probs_region_surveilled_hospitals),
                      n_t=n_t_array,
                      n_p = n_p_array,
                      n_pa = n_pa_array,
                      L=length(car_sim_region$nb_info$weights),
                      rain=rain_matrix[,rain_index_indices],
                      adj = car_sim_region$nb_info$adj,
                      weights = car_sim_region$nb_info$weights,
                      num = car_sim_region$nb_info$num,
                      distance=hosp_probs_sim$distance,
                      surveilled_indices=surveilled_indices,
                      not_surveilled_indices=not_surveilled_indices,
                      mu_phi_binom=rep(0,64),
                      max_n_d=max_n_d,lgamma_table=lgamma(0:1000+1))
return(list(data_list=data_list,constant_list=constant_list,d_2025=true_cases_2025$d,
            mu_forecast=region_level_data$population/100000*exp(log_lambda_matrix[,length(years)+1]),pi_df=pi_df))
}

sim_results<-data_gen(beta_0=mean_log_rate,beta_1=log(1.05),r=120,mean_prob=0.33,sigma_Phi=0.5,
                      sigma_phi=0.5,zeta=0.06,
    multiplier_lower=4,multiplier_upper=8,sensitivity=c(0.6,0.9),
    specificity=c(0.9,0.9),test_receipt_probs=c(0.5,0.7),
    n_h=0,n_years=3,rain_matrix=rain_matrix)

#Create distribution for vector of test results

#First create necessary nimble functions

#First create necessary nimble functions

logsumexp_vec <- nimbleFunction(
  run = function(vals = double(1)) {
    returnType(double())
    m <- max(vals)
    return(m + log(sum(exp(vals - m))))
  }
)

log_lik_Y <- nimbleFunction(
  run = function(sens = double(1),
                 spec = double(1),
                 n_t = double(0),
                 test_mat = double(2)) {
    returnType(double(1))
    
    
    # Step 1: handle the case of no non-padded rows
    if(n_t == 0) {
      return(rep(0, 1))  # log-likelihood = log(1) = 0
    }
    
    n_total <- n_t
    K <- length(sens)
    
    # Step 2: split passive and active
    passive_idx <- which(test_mat[, 1] == 1)
    active_idx  <- which(test_mat[, 1] == 0)
    
    n_pa <- length(passive_idx)
    nRemain <- length(active_idx)
    
    if(n_pa > 0) {
      Y_passive <- test_mat[passive_idx, 2:(K+1), drop = FALSE]
    }
    if(nRemain > 0) {
      Y_active  <- test_mat[active_idx, 2:(K+1), drop = FALSE]
    }
    
    # Step 3: compute log-term for passive observations
    log_term_passive <- 0
    if (n_pa > 0) {
      for (j in 1:n_pa) {
        for (k in 1:K) {
          y <- Y_passive[j, k]
          if(y==-1){
            log_term_passive <- log_term_passive + 0
          }
          else{
            log_term_passive <- log_term_passive + (y * log(sens[k]) + (1 - y) * log(1 - sens[k]))
          }
        }
      }
    }
    
    if(nRemain == 0) {
      return(rep(log_term_passive,1))
    }
    
    # Step 4: compute log-weights for active observations
    log_w1 <- numeric(nRemain)
    log_w0 <- numeric(nRemain)
    
    for (j in 1:nRemain) {
      lw1 <- 0
      lw0 <- 0
      for (k in 1:K) {
        y <- Y_active[j, k]
        if(y==-1){
          lw1 <- lw1 + 0
          lw0<- lw0 + 0
        }
        
        else{
          lw1 <- lw1 + (y * log(sens[k]) + (1 - y) * log(1 - sens[k]))
          lw0 <- lw0 + (y * log(1 - spec[k]) + (1 - y) * log(spec[k]))
        }
      }
      log_w1[j] <- lw1
      log_w0[j] <- lw0
    }
    
    # Step 5: dynamic programming matrix
    log_xi <- matrix(-Inf, nrow = nRemain + 1, ncol = nRemain + 1)
    log_xi[1, 1] <- 0
    
    for (l in 1:nRemain) {
      for (m in 0:l) {
        nVals <- 0
        vals <- numeric(2)
        
        vals[1] <- log_w0[l] + log_xi[m + 1, l]
        nVals <- nVals + 1
        
        if (m > 0) {
          vals[nVals + 1] <- log_w1[l] + log_xi[m, l]
          nVals <- nVals + 1
        }
        
        if (nVals > 0) {
          log_xi[m + 1, l + 1] <- logsumexp_vec(vals[1:nVals])
        }
      }
    }
    
    # Step 6: log-likelihood vector
    log_likelihood <- numeric(nRemain + 1)
    for (k in 0:nRemain) {
      log_likelihood[k + 1] <- log_xi[k + 1, nRemain + 1] + log_term_passive -
        (lfactorial(nRemain) - lfactorial(k) - lfactorial(nRemain - k))
    }
    
    return(log_likelihood)
  }
)

## ==========================================================
## NIMBLE wrappers for Negative Binomial (size = r, prob = p)
## ==========================================================

## --- Fully compiled Negative Binomial density ---
dMarginalized <- nimbleFunction(
  run = function(x = double(2),
                 n_pa = double(0),   # n^{PA}
                 n_t  = double(0),   # n^{T}
                 n_p  = double(0),   # n^{P}
                 nu   = double(1),
                 kappa = double(1),
                 piP = double(0),
                 piA = double(0),
                 mu_D = double(0),
                 r=double(0),
                 max_nD = double(0, default = 200),
                 lgamma_table=double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    if(n_pa < 0 | n_t < n_pa) {
      if(log) return(-Inf) else return(0.0)
    }
    
    ## likelihood of Y given nDT
    log_lik <- log_lik_Y(sens = nu, spec = kappa, n_t=n_t, test_mat = x)
    total_log_prob <- -Inf
    
    ## precompute logs
    log_piP <- log(piP)
    log1m_piP <- log(1 - piP)
    log_piA <- log(piA)
    log1m_piA <- log(1 - piA)
    
    #Pre-compute gamma values
    
    a1 <- n_pa+1 
    a2 <- n_p - n_pa+1
    lgamma_a1<-lgamma_table[a1]
    lgamma_a2<-lgamma_table[a2]
    
    p <- r/(r + mu_D)
    
    q<-qnbinom(p=1 - 10^(-15), prob=p,size=r)
    
    upper_nDT <- max(min(n_t, q),n_pa)
    
    
    ## sum over nDT
    for(nDT in n_pa:upper_nDT){
      idx <- nDT - n_pa + 1
      log_prob_Y_given_nDT <- log_lik[idx]
      
      nD_min <- nDT + n_p - n_pa
      inner_log_sum <- -Inf
      a3 <- nDT - n_pa+1
      lgamma_a3<-lgamma_table[a3]
      upper_nD <- max(min(max_nD, q),nD_min)
      
      for(nD in nD_min:upper_nD) {
        
        a4 <- nD - (n_p + nDT - n_pa) + 1
        lgamma_a4 <- lgamma_table[a4]
        log_mult_coeff <- lgamma_table[nD + 1] -
          (lgamma_a1 +
             lgamma_a2 +
             lgamma_a3 +
             lgamma_a4)
        
        log_prob_cells <- (a1-1) * (log_piP + log_piA) +
          (a2-1) * (log_piP + log1m_piA) +
          (a3-1) * (log1m_piP + log_piA) +
          (a4-1) * (log1m_piP + log1m_piA)
        
        log_p_nD <- dnbinom(x=nD, prob=p, size=r, log = TRUE)
        
        cur_log <- log_mult_coeff + log_prob_cells +
          log_prob_Y_given_nDT + log_p_nD
        
        ## log-sum-exp accumulate
        if(inner_log_sum == -Inf) inner_log_sum <- cur_log
        else {
          if(cur_log > inner_log_sum)
            inner_log_sum <- cur_log + log1p(exp(inner_log_sum - cur_log))
          else
            inner_log_sum <- inner_log_sum + log1p(exp(cur_log - inner_log_sum))
        }
      }
      
      
      ## accumulate across nDT
      if(inner_log_sum > -Inf) {
        if(total_log_prob == -Inf) total_log_prob <- inner_log_sum
        else {
          if(inner_log_sum > total_log_prob)
            total_log_prob <- inner_log_sum + log1p(exp(total_log_prob - inner_log_sum))
          else
            total_log_prob <- total_log_prob + log1p(exp(inner_log_sum - total_log_prob))
        }
      }
    }
    
    if(log) return(total_log_prob)
    else return(exp(total_log_prob))
  }
)


registerDistributions(list(
  dMarginalized = list(
    BUGSdist = "dMarginalized(n_pa, n_t, n_p, nu, kappa, piP, piA, mu_D,r, max_nD,
    lgamma_table)",
    types = c(
      "value = double(2)",
      "n_pa = double(0)",
      "n_t  = double(0)",
      "n_p  = double(0)",
      "nu   = double(1)",
      "kappa = double(1)",
      "piP = double(0)",
      "piA = double(0)",
      "mu_D = double(0)",
      "r = double(0)",
      "max_nD = double(0)",
      "lgamma_table = double(1)"
    ),
    discrete = TRUE
  )
))

#Write exponential covariance matrix

expcov <- nimbleFunction(
  run = function(distance = double(2), zeta = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(distance)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma * sigma
    for(i in 1:n) {
      for(j in 1:n) {
        result[i, j] <- sigma2 * exp(- distance[i, j] * zeta)
      }
    }
    return(result)
  }
)


#Write model code

m_code <- nimbleCode({
  
  for(h in 1:H){
    phi_binom[h] <- sigma_phi_binom * phi_binom_raw[h]
    pi_P[h] <- expit(alpha_0 + phi_binom[h])
  }
  
  for (h in 1:H_t) {
    pi_P_s[h] <- pi_P[surveilled_indices[h]]
  }
  
  for (h in 1:Not_H_t) {
    pi_P_ns[h] <- pi_P[not_surveilled_indices[h]]
  }
  
  for(s in 1:S){
    Phi_pois[s] <- sigma_Phi_pois * Phi_pois_raw[s] 
  }

for(t in 1:T){
  for(s in 1:S){
    log_lambda[s,t] <- beta_0 + beta_1*rain[s,t] + Phi_pois[s]
    mu[s,t]<-pop[s]*exp(log_lambda[s,t])*sum(prob_pres_ns[s,]*pi_P_ns[])
    p[s,t] <- r / (r + mu[s,t])
    y[s,t] ~ dnegbin(prob=p[s,t],size=r)
  }
}
  #Sampled hospitals
  
for(t in 1:T){
  for(h in 1:H_t){
    for(s in 1:S){
  mu_d[s,h,t]<-pop[s]*exp(log_lambda[s,t])*prob_pres_s[s,h]
  Y[,,(t-1)*(H_t*S)+(h-1)*S+s] ~ dMarginalized(n_pa=n_pa[s,h,t],
                                               n_t=n_t[s,h,t],
                                               n_p=n_p[s,h,t],
                                               nu=nu[1:2],
                                               kappa=kappa[1:2],
                                               piP=pi_P_s[h],
                                               piA=pi_A,
                                               mu_D=mu_d[s,h,t],
                                               r=r,
                                               max_nD=max_n_d,
                                               lgamma_table=lgamma_table[])
    }
  } 
}
  
for(s in 1:S){
log_lambda_forecast[s] <- beta_0 + beta_1*rain[s,T+1] + Phi_pois[s]
mu_forecast[s]<-pop[s]*exp(log_lambda_forecast[s])
}
  
  #Priors
  beta_0 ~ dnorm(0, sd = 5)
  beta_1 ~ dnorm(0, sd=2)
  Phi_pois_raw[1:S] ~ dcar_normal(adj[1:L], weights[1:L], num[1:S], 1, zero_mean = 1)
  sigma_Phi_pois ~ dunif(0,2)
  alpha_0 ~ dnorm(0, sd=1.5)
  Corr[1:H,1:H]<-expcov(distance[1:H,1:H], zeta, 1)
  phi_binom_raw[1:H] ~ dmnorm(mean=mu_phi_binom[],cov=Corr[1:H,1:H])
  log_sigma_phi_binom ~ dnorm(0,1)
  sigma_phi_binom <- exp(log_sigma_phi_binom)
  zeta <- 0.06
  nu[1] ~ dbeta(6000,4000)
  nu[2] ~ dbeta(9000,1000)
  kappa[1] ~ dbeta(9000,1000)
  kappa[2] ~ dbeta(9000,1000)
  pi_A ~ dbeta(50,50)
  r <- exp(log_r)
  log_r ~ dnorm(4.5,1)
})

#Fit and run model

n_p <- sim_results$constant_list$n_p
n_pa <- sim_results$constant_list$n_pa
n_t <- sim_results$constant_list$n_t

# Create function to initialize Phi_pois and sigma_phi_pois

pois_init <- function(upper_pois) {
  sigma_Phi_pois <- runif(1, 0.1, upper_pois)
  Phi_pois_raw <- rnorm(7, 0, 1)
  return(list(sigma_Phi_pois = sigma_Phi_pois, Phi_pois_raw = Phi_pois_raw))
}

# Create function to initialize phi_binom_s, phi_binom_ns, and sigma_phi_binom

binom_init <- function(sd_binom) {
  log_sigma_phi_binom <- rnorm(1,0,sd=sd_binom)
  phi_binom_raw <- rnorm(sim_results$constant_list$H, 0, 1)
  return(list(log_sigma_phi_binom = log_sigma_phi_binom, phi_binom_raw=phi_binom_raw))
}


# Create final inits function

inits_function <- function(upper_pois, sd_binom) {

pois <- pois_init(upper_pois = upper_pois)
binom <- binom_init(sd_binom=sd_binom)
  return(list(
    sigma_Phi_pois = pois$sigma_Phi_pois,
    Phi_pois_raw = pois$Phi_pois_raw,
    log_sigma_phi_binom = binom$log_sigma_phi_binom,
    phi_binom_raw = binom$phi_binom_raw,
    beta_0 = rnorm(1, 0, 1),
    beta_1= rnorm(1, 0, 1),
    alpha_0 = rnorm(1, 0, 1),
    pi_A=rbeta(1,1,1),
    log_r = rnorm(1,4.5,1)
  ))
}

# Initialize and configure model

Rmodel <- nimbleModel(
  m_code,
  constants = sim_results$constant_list,
  data = sim_results$data_list,
  dimensions = list(mu_phi_binom = 64, cov = c(64,64),
                    prob_pres_s = c(7, sim_results$constant_list$H_t),       # or (S, H_t) depending on how you define it
                    prob_pres_ns = c(7,sim_results$constant_list$Not_H_t),
                    pi_P_ns = sim_results$constant_list$Not_H_t)
)



conf <- configureMCMC(Rmodel)
conf$addMonitors(c("mu_forecast","r"))


Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)



library(parallel)
library(coda)  # if you want to keep using mcmc objects

# Function to run a single chain
run_chain <- function(chain_id, nburnin, niter) {
  set.seed(123 + chain_id)  # optional: ensure reproducibility
  inits <- inits_function(upper_pois = 0.5, sd_binom = 0.5)
  
  chain_samples <- runMCMC(
    Cmcmc,
    nburnin = nburnin,
    niter = niter,
    nchains = 1,           # only one chain per mclapply call
    inits = inits,
    samplesAsCodaMCMC = TRUE
  )
  
  return(chain_samples)
}

# Number of chains to run
n_chains <- 4
nburnin <- 40000
niter <- 80000

# Run chains in parallel
samples_list <- mclapply(
  1:n_chains,
  run_chain,
  nburnin = nburnin,
  niter = niter,
  mc.cores = n_chains
)

combined_samples <- do.call(coda::mcmc.list, samples_list)

plot(combined_samples[,c("beta_0","beta_1","sigma_Phi_pois")],density=F)
plot(combined_samples[,c("alpha_0","log_sigma_phi_binom","pi_A")],density=F)
plot(combined_samples[,c("mu_forecast[1]","mu_forecast[2]","mu_forecast[3]","mu_forecast[4]")],density=F)
plot(combined_samples[,c("mu_forecast[5]","mu_forecast[6]","mu_forecast[7]","r")],density=F)
summary(combined_samples[,c("beta_0","beta_1","sigma_Phi_pois",
                            "alpha_0","log_sigma_phi_binom","pi_A",
                            "mu_forecast[1]","mu_forecast[2]","mu_forecast[3]","mu_forecast[4]",
                        "mu_forecast[5]","mu_forecast[6]","mu_forecast[7]","r")])


# Combine chains into an mcmc.list
combined_chains <- mcmc.list(combined_samples[,c("beta_0","beta_1","sigma_Phi_pois",
                                                 "alpha_0","log_sigma_phi_binom","pi_A",
                                                 "mu_forecast[1]","mu_forecast[2]","mu_forecast[3]","mu_forecast[4]",
                                                 "mu_forecast[5]","mu_forecast[6]","mu_forecast[7]","r")])

# Calculate R-hat
gelman.diag(combined_chains)

#Calculate ESS

ess <- effectiveSize(combined_samples)
