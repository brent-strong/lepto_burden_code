library(tidyverse)
library(nimble)
library(parallel)
library(coda)
library(purrr)
library(sf)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))
source(repo_path("R", "surveillance_inputs.R"))
setwd(analysis_path())
example <- read_rds("example_inputs.rds")

hospitals <- st_read(data_path("raw", "geospatial", "hospitals", "Hospitals_PR.shp")) %>%
  filter(NAICS_DESC=="GENERAL MEDICAL AND SURGICAL HOSPITALS")


# Set the path to your folder
folder_path <- "PR_regions/"

# Read the specific shapefile
# Note: Point to the .shp file, but the .dbf, .shx, and .prj must be in the same folder.
regions <- st_read(paste0(folder_path, "pr_health_regions.shp"))

# Custom inverse logit function
inv_logit <- function(x) {
  return(1 / (1 + exp(-x)))
}


#Process incidence data

library(tidyverse)

# 1. Load the data
data <- read_csv("Incidence_data.csv")

# 2. Process and transform the data
y <- data %>%
  # Remove "Imported" cases
  filter(!str_detect(`Health Region`, "Imported")) %>%
  # Standardize region names to ensure 2022 and 2023 match (7 regions total)
  mutate(`Health Region` = recode(`Health Region`, "Metropolitana" = "Metro")) %>%
  # Calculate new total using only Confirmed and Probable cases
  mutate(Case_Sum = `Confirmed Cases` + `Probable Cases`) %>%
  # Keep only the columns needed for the matrix
  dplyr::select(Year, `Health Region`, Case_Sum) %>%
  # Pivot the years into columns
  pivot_wider(names_from = Year, values_from = Case_Sum) %>%
  # Move Health Region names to the matrix row names
  column_to_rownames("Health Region") %>%
  # Convert to a formal matrix object
  as.matrix()

# 3. Display the result
print(y)

#Extract population data

pop <- example$data_list$pop

#Get presentation probabilities

prob_pres <- read.csv("health_regions_probs.csv") %>% dplyr::select(-region)
prob_pres_s <- prob_pres %>% dplyr::select(hosp_2512,hosp_14,hosp_15,hosp_23)


col_names <- c("hosp_2512", "hosp_14", "hosp_15", "hosp_23")
surveilled_indices <- match(col_names, names(prob_pres))
all_indices <- 1:ncol(prob_pres)
not_surveilled_indices <- setdiff(all_indices, surveilled_indices)

surveillance_inputs <- load_surveillance_model_inputs()
n_t <- surveillance_inputs$n_t
n_p <- surveillance_inputs$n_p
n_pa <- surveillance_inputs$n_pa
Y <- surveillance_inputs$Y

# 1. Get dimensions of your original array Y
n_rows  <- dim(Y)[1]
n_depth <- dim(Y)[3]

# 2. Pre-allocate the new array with 2 columns
# Dimensions: [Rows, 2, Depth]
new_array <- array(0, dim = c(n_rows, 3, n_depth))

# 3. Loop through the depth to transform each matrix
for (i in 1:n_depth) {
  # Extract the i-th matrix
  mat <- Y[,,i]
  
  # Column 1: Keep the original ID/Index column
  col1 <- mat[, 1]
  
  # Columns 2-4: The data columns
  data_cols <- mat[, 2:3]
  
  # Logic: Check if any of the three data columns equal 1
  # This creates a vector of 1s and 0s
  combined_val <- as.numeric(rowSums(data_cols == 1) > 0)
  
  # Preserve -999 values
  # If the first data column is -999, set the result to -999
  is_missing <- data_cols[, 1] == -999
  combined_val[is_missing] <- -999
  
  # Assign the 2nd column back to the new array
  new_array[, 1, i] <- col1
  new_array[, 2, i] <- combined_val
  new_array[, 3, i] <- -1
  
  is_missing <- data_cols[, 1] == -999
  new_array[is_missing, 2, i] <- -999
  new_array[is_missing, 3, i] <- -999
}

data_list <- list(y=y,pop=pop,Y=new_array)
constant_list <- example$constant_list
constant_list$t_adj <- c(11.5/12,13/12,18/12,1.5/12)
constant_list$n_p <- n_p
constant_list$n_pa <- n_pa
constant_list$n_t <- n_t
constant_list$prob_pres <- prob_pres
constant_list$prob_pres_s <- prob_pres_s
constant_list$prob_pres_ns <- NULL
constant_list$T <- 2
constant_list$y_2024 <- NULL
constant_list$rain <- NULL
constant_list$surveilled_indices<-surveilled_indices
constant_list$not_surveilled_indices<-not_surveilled_indices
constant_list$max_n_d <- 300


#Create distribution for vector of test results

#First create necessary nimble functions

#First create necessary nimble functions

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
          lw1 <- lw1 + (y * log(sens[k]) + (1 - y) * log(1 - sens[k]+10^(-30)))
          lw0 <- lw0 + (y * log(1 - spec[k]+10^(-30)) + (1 - y) * log(spec[k]))
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
    
    
    q<-qpois(p=1 - 10^(-9), lambda=mu_D)
    
    upper_nDT <- n_t
    
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
        
        log_p_nD <- dpois(x=nD, lambda=mu_D, log = TRUE)
        
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
    BUGSdist = "dMarginalized(n_pa, n_t, n_p, nu, kappa, piP, piA, mu_D, max_nD,
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
    tau_raw[s] ~ dnorm(0,1)
    tau[s] <- sigma_tau * tau_raw[s]
    log_lambda[s] <- A + Phi_pois[s] + tau[s]
    mu[s]<-pop[s]*exp(log_lambda[s])*sum(prob_pres[s,]*pi_P[1:H])
  }
  
  for(t in 1:T){
    for(s in 1:S){
      y[s,t] ~ dpois(mu[s])
    }
  }
  #Sampled hospitals
  
  for(h in 1:H_t){
    for(s in 1:S){
      mu_d[s,h]<-pop[s]*exp(log_lambda[s])*prob_pres_s[s,h]*t_adj[h]
      Y[,,(h-1)*S+s] ~ dMarginalized(n_pa=n_pa[s,h],
                                     n_t=n_t[s,h],
                                     n_p=n_p[s,h],
                                     nu=nu[1:2],
                                     kappa=kappa[1:2],
                                     piP=pi_P_s[h],
                                     piA=pi_A,
                                     mu_D=mu_d[s,h],
                                     max_nD=max_n_d,
                                     lgamma_table=lgamma_table[])
    }
  } 
  
  
  #Forecasts  
  for(s in 1:S){
    mu_forecast[s]<-pop[s]*exp(log_lambda[s])*sum(prob_pres[s,]*(1-pi_P[1:H]))
    d_y_2022_forecast[s] ~ dpois(mu_forecast[s])
    rate_2022[s] <- (d_y_2022_forecast[s] + y[s,1])/pop[s]
    d_y_2023_forecast[s] ~ dpois(mu_forecast[s])
    rate_2023[s] <- (d_y_2023_forecast[s] + y[s,2])/pop[s]
  }
  
  #Forecast rate for entire island
  
  rate_PR_2022 <- sum(d_y_2022_forecast[1:S] + y[,1])/sum(pop[1:S])
  rate_PR_2023 <- sum(d_y_2023_forecast[1:S] + y[,2])/sum(pop[1:S])
  
  #Priors
  A ~ dnorm(3.7, sd = 1.5)
  Phi_pois_raw[1:S] ~ dcar_normal(adj[1:L], weights[1:L], num[1:S], 1, zero_mean = 1)
  log_sigma_Phi_pois ~ dnorm(0, sd=1)
  sigma_Phi_pois <- exp(log_sigma_Phi_pois)
  log_sigma_tau ~ dnorm(0, sd=1)
  sigma_tau <- exp(log_sigma_tau)
  alpha_0 ~ dnorm(-1.4, sd=0.5)
  Corr[1:H,1:H]<-expcov(distance[1:H,1:H], zeta, 1)
  phi_binom_raw[1:H] ~ dmnorm(mean=mu_phi_binom[1:H],cov=Corr[1:H,1:H])
  log_sigma_phi_binom ~ dnorm(0,sd=1)
  sigma_phi_binom <- exp(log_sigma_phi_binom)
  zeta <- 0.06
  pi_A ~ dbeta(5,20)
  nu[1]<-1
  nu[2] <-1
  kappa[1]<-1
  kappa[2]<-1
})

# Create function to initialize Phi_pois and sigma_phi_pois

pois_init <- function(sd_pois) {
  log_sigma_Phi_pois <- rnorm(1, 0, sd=sd_pois)
  Phi_pois_raw <- rnorm(7, 0, 1)
  log_sigma_tau <- rnorm(1, 0, sd=sd_pois)
  tau_raw <- rnorm(7,0,1)
  return(list(log_sigma_Phi_pois = log_sigma_Phi_pois, Phi_pois_raw = Phi_pois_raw,
              log_sigma_tau=log_sigma_tau, tau_raw=tau_raw
  ))
}

# Create function to initialize phi_binom_s, phi_binom_ns, and sigma_phi_binom

binom_init <- function(sd_binom) {
  log_sigma_phi_binom <- rnorm(1,0,sd=sd_binom)
  phi_binom_raw <- rnorm(constant_list$H, 0, 1)
  return(list(log_sigma_phi_binom = log_sigma_phi_binom, phi_binom_raw=phi_binom_raw))
}


# Create final inits function

inits_function <- function(sd_pois, sd_binom) {
  
  pois <- pois_init(sd_pois = sd_pois)
  binom <- binom_init(sd_binom=sd_binom)
  return(list(
    log_sigma_Phi_pois = pois$log_sigma_Phi_pois,
    Phi_pois_raw = pois$Phi_pois_raw,
    log_sigma_tau = pois$log_sigma_tau,
    tau_raw = pois$tau_raw,
    log_sigma_phi_binom = binom$log_sigma_phi_binom,
    phi_binom_raw = binom$phi_binom_raw,
    A = rnorm(1, 0, 1),
    alpha_0 = rnorm(1, 0, 1),
    pi_A=rbeta(1,1,1)
  ))
}

#Create functions to run a single chain


nburnin <- 0
niter <- 500


Rmodel <- nimbleModel(
  m_code,
  constants = constant_list,
  data = data_list,
  dimensions = list(mu_phi_binom = 64,cov = c(64,64),
                    prob_pres_s = c(7, constant_list$H_t),       # or (S, H_t) depending on how you define it
                    pi_P_ns = constant_list$Not_H_t)
)



conf <- configureMCMC(Rmodel)
conf$removeSamplers(c("phi_binom_raw[1:64]"))
conf$addSampler("phi_binom_raw[1:64]", type='ess')
conf$addMonitors(c("rate_2023", "rate_PR_2023","rate_2022", "rate_PR_2022"))

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(123)

library(parallel)

# 1. Define the number of cores (chains)
n_cores <- 3 

# 2. Run in parallel
chain_samples_list <- mclapply(1:n_cores, function(i) {
  
  # It is good practice to provide different inits or seeds per chain
  # Assuming inits_function takes a seed or index
  current_inits <- inits_function(1, 1) 
  
  # Run the MCMC
  # Note: nchains is set to 1 because each core handles one chain
  runMCMC(
    Cmcmc, 
    nburnin = nburnin, 
    niter = niter, 
    nchains = 1, 
    inits = current_inits, 
    samplesAsCodaMCMC = TRUE
  )
}, mc.cores = n_cores)

# 3. Combine into a single mcmc.list (coda package)
library(coda)
chain_samples <- as.mcmc.list(chain_samples_list)

write_rds(chain_samples, file = analysis_path("Archive", "results_sensitivity.rds"))

# Compute R-hat
library(coda)
rhat_results <- gelman.diag(chain_samples, multivariate=F)

# Extract point estimates
rhat_point_estimates <- rhat_results$psrf[, "Point est."]

plot1<-plot(chain_samples[,c("alpha_0","rate_2023[1]","rate_2023[2]","rate_2023[3]")],density=F)
plot2<-plot(chain_samples[,c("rate_2023[4]","rate_2023[5]","rate_2023[6]",
                             "rate_2023[7]")],density=F)



# now summarize, including the new column
summary(
  chain_samples[, c(
    "A","alpha_0",
    "rate_2023[1]","rate_2023[2]","rate_2023[3]",
    "rate_2023[4]","rate_2023[5]","rate_2023[6]",
    "rate_2023[7]","rate_PR_2023","rate_PR_2022",
    "log_sigma_Phi_pois","log_sigma_phi_binom",
    "pi_A","log_sigma_tau"
  )]
)


##Create figures


library(viridis)

# 1. Load the shapefile
regions <- st_read("PR_regions/pr_health_regions.shp")

# 2. Crop the geometry to remove Mona Island
# Mona is located roughly at -67.9. We crop from -67.3 (west of main island) 
# to -65.2 (east of Culebra).
# Note: This assumes your shapefile is in WGS84 (long/lat). 
# If it's in a different projection, run: regions <- st_transform(regions, 4326)
regions <- st_transform(regions, 4326)
regions_cropped <- st_crop(regions, xmin = -67.3, xmax = -65.2, ymin = 17.8, ymax = 18.6)

# 3. Prepare your posterior data 
# (Assuming summ and summ_quant are extracted from your chain_samples summary)
summ <- summary(chain_samples)$statistics
summ_quant <- summary(chain_samples)$quantiles

model_results <- data.frame(
  region_id = 1:7,
  Observed = y[,2]/pop,
  Posterior.Mean = summ[paste0("rate_2023[", 1:7, "]"), "Mean"],
  Lower.PI.Bound = summ_quant[paste0("rate_2023[", 1:7, "]"), "2.5%"],
  Upper.PI.Bound = summ_quant[paste0("rate_2023[", 1:7, "]"), "97.5%"]
)

# 4. Join and Pivot to Long Format
map_data <- cbind(regions_cropped, model_results)

map_long <- map_data %>%
  pivot_longer(
    cols = c(Observed, Posterior.Mean, Lower.PI.Bound, Upper.PI.Bound),
    names_to = "Metric",
    values_to = "Rate"
  ) %>%
  mutate(Metric = case_match(Metric,
                             "Observed"        ~ "Observed",
                             "Posterior.Mean"  ~ "Posterior Mean",
                             "Lower.PI.Bound"  ~ "Lower CI Bound",
                             "Upper.PI.Bound"  ~ "Upper CI Bound"
  )) %>%
  mutate(Metric = factor(Metric, levels = c(
    "Observed", 
    "Posterior Mean", 
    "Lower CI Bound", 
    "Upper CI Bound"
  )))

# 5. Create the 4-panel Figure
p1 <- ggplot(data = map_long) +
  geom_sf(aes(fill = Rate), color = "black", size = 0.1) +
  geom_sf_text(aes(label = region), color = "white", size = 1.5) +
  facet_wrap(~Metric, ncol = 2) +
  scale_fill_viridis_c(
    option = "viridis", 
    name = "Rate",
    limits = c(0, max(map_long$Rate, na.rm = TRUE))
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )


ggsave(
  filename = "PR_Health_Regions_Map_sensitivity.png", 
  plot = p1, 
  width = 10,       # Width in inches
  height = 8,       # Height in inches
  dpi = 300,        # High resolution
  bg = "white"      # Ensures the background isn't transparent
)

library(magick)
img <- image_read("PR_Health_Regions_Map_sensitivity.png") %>%
  image_trim() %>%
  image_border(color = "white", geometry = "50x50") # Adds 50 pixels of padding

image_write(img, "PR_Health_Regions_Map_Cropped_sensitivity.png")
