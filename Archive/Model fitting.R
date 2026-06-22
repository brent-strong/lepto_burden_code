library(tidyverse)
library(nimble)
library(nimbleHMC)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))
source(data_path("scripts", "Simulation_Code.R"))

#Get hospital presentation probabilities

pres_probs_region<-as.data.frame(read_csv(data_path("processed", "simulation_inputs", "health_regions_probs.csv")))

#Get population in regions

region_level_data<-as.data.frame(read_csv(data_path("processed", "simulation_inputs", "region_population.csv")))

region_level_data$lepto_cases_2022 <- c(32,37,53,0,55,31,39)
region_level_data$lepto_cases_2023 <- c(42,36,80,4,47,26,47)

#Get data for each 1 km x 1 km cells

cell_data<-as.data.frame(read_csv(data_path("processed", "simulation_inputs", "cell_probabilities.csv")))

#Get rainfall data

rainfall<-as.data.frame(read_csv(data_path("processed", "simulation_inputs", "rainfall2023.csv")))

#Create rates per 100,000

region_level_data <- region_level_data %>% mutate(lepto_per_100k_2022 = lepto_cases_2022*100000/population,
  lepto_per_100k_2023 = lepto_cases_2023*100000/population) %>%
  mutate(lepto_per_100k_2022=ifelse(region=="Fajardo",0.5,lepto_per_100k_2022))

mean_log_rate <- mean(log(c(3*region_level_data$lepto_per_100k,3*region_level_data$lepto_per_100k_2022)))
sd_log_rate <- sd(log(c(3*region_level_data$lepto_per_100k,3*region_level_data$lepto_per_100k_2022)))

lambda <- exp(rnorm(nrow(region_level_data),mean_log_rate,sd_log_rate))

lambda_df <- data.frame(health_region=region_level_data$region,lambda=lambda)

lambda_full <- left_join(cell_data,lambda_df,by="health_region")$lambda

#Create capture probabilities

#Create logit function

logit <- function(p) {
  if (any(p <= 0 | p >= 1, na.rm = TRUE)) warning("Input probabilities must be in (0, 1)")
  log(p / (1 - p))
}

#Create inverse logit function

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

mean_logit_probability <- logit(0.2)
sd_logit_probability <- 0.5

#Create probabilities

pi <- inv_logit(rnorm(ncol(pres_probs_region)-1,mean_logit_probability,sd_logit_probability))
pi_df <- data.frame(hospital=colnames(pres_probs_region)[-1],pi=pi)

link_df <- tidyr::crossing(region=region_level_data$region,hospital=pi_df$hospital) %>% 
  left_join(pi_df)

#Create data set

data <- data_generation(lambda=lambda_full,rho=lambda_full,pi=link_df$pi,multiplier_lower=4,multiplier_upper=8,
                        sensitivity=c(0.6,0.9),specificity=c(0.9,0.9),test_receipt_probs=c(0.5,0.3))

#Surveil a set of hospitals

# Set the proportion to sample
p <- 0.15625

surveillance <- data$hospital_data %>%
  filter(hospital %in% sample(unique(hospital), size = ceiling(p * n_distinct(hospital)))) %>%
  mutate(passive_surveillance_capture_mod=ifelse(passive_surveillance_capture!=-1,
                                                 passive_surveillance_capture,0))  %>%
  mutate(tested_true_passive=ifelse((test1==-1 & test2==-1) | passive_surveillance_capture_mod==0,0,1)) %>%
  mutate(tested=ifelse(test1==-1 & test2==-1,0,1)) %>%
  mutate(tested_true=ifelse((test1==-1 & test2==-1) | disease_status==0,0,1)) %>%
  group_by(hospital) %>%
  arrange(desc(passive_surveillance_capture), .by_group=T) %>%
  ungroup()

#Create df to split

filtered_surveillance <- surveillance %>% filter(tested==1)

split_df <- split(filtered_surveillance[,4:6], filtered_surveillance$hospital)
mat_conv <- function(matrix){
  matrix <- as.matrix(matrix)
  return(matrix)
}

# Step 2: Convert each subset to a matrix
matrix_test <- lapply(split_df, mat_conv)


#Get vector of cases at surveilled hospitals and total number tested

cases_passive<-surveillance %>% group_by(hospital) %>% 
  summarize(n_r=sum(passive_surveillance_capture_mod),
            n_rt=sum(tested_true_passive))
n_r <- cases_passive$n_r
n_rt <- cases_passive$n_rt

cases_total<-surveillance %>% group_by(hospital) %>% 
  summarize(n_t=sum(tested))
n_t <- cases_total$n_t

cases_true<-surveillance %>% group_by(hospital) %>%
  summarize(n_d=sum(disease_status),
            n_dt=sum(tested_true))
n_d <- cases_true$n_d
n_dt <- cases_true$n_dt

#Create data for the input into model

passive_cases_surveilled_hospitals<-surveillance %>%
  group_by(health_region) %>%
  summarize(y=sum(passive_surveillance_capture_mod))

region_cases <- data$region_level_data
y <- region_cases$y - passive_cases_surveilled_hospitals$y

#Remove columns for processing

pres_probs_region_non_surveilled_hospitals <- pres_probs_region %>%
  dplyr::select(-all_of(surveillance$hospital))

pres_probs_region_surveilled_hospitals <- pres_probs_region %>%
  dplyr::select(all_of(surveillance$hospital))

#Create distribution for vector of test results

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
                 test_mat = double(2)) {
    returnType(double(1))
    
    n_total <- dim(test_mat)[1]
    K <- length(sens)
    
    passive_idx <- which(test_mat[, 1] == 1)
    active_idx  <- which(test_mat[, 1] == 0 | test_mat[, 1] == -1)
    
    n_rt <- length(passive_idx)
    nRemain <- length(active_idx)
    
    if(n_rt > 0) {
      Y_passive <- test_mat[passive_idx, -1, drop = FALSE]
    } else {
      Y_passive <- matrix(0, nrow = 0, ncol = K)
    }
    
    if(nRemain > 0) {
      Y_active  <- test_mat[active_idx, -1, drop = FALSE]
    } else {
      Y_active <- matrix(0, nrow = 0, ncol = K)
    }
    
    log_term_passive <- 0
    if (n_rt > 0) {
      for (j in 1:n_rt) {
        for (k in 1:K) {
          y <- Y_passive[j, k]
          log_term_passive <- log_term_passive + (1 - (y < 0)) * (y * log(sens[k]) + (1 - y) * log(1 - sens[k]))
        }
      }
    }
    
    log_w1 <- numeric(nRemain)
    log_w0 <- numeric(nRemain)
    
    for (j in 1:nRemain) {
      lw1 <- 0
      lw0 <- 0
      for (k in 1:K) {
        y <- Y_active[j, k]
        lw1 <- lw1 + (1 - (y < 0)) * (y * log(sens[k]) + (1 - y) * log(1 - sens[k]))
        lw0 <- lw0 + (1 - (y < 0)) * (y * log(1 - spec[k]) + (1 - y) * log(spec[k]))
      }
      log_w1[j] <- lw1
      log_w0[j] <- lw0
    }
    
    log_xi <- matrix(-1e12, nrow = nRemain + 1, ncol = nRemain + 1)
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
    
    log_likelihood <- numeric(nRemain + 1)
    for (k in 0:nRemain) {
      log_likelihood[k + 1] <- log_xi[k + 1, nRemain + 1] + log_term_passive - (lfactorial(nRemain) - lfactorial(k) - lfactorial(nRemain - k))
    }
    
    return(log_likelihood)
  }
)

dTests <- nimbleFunction(
  run = function(x = double(2),
                 n_dt=double(0),
                 n_rt=double(0),
                 n_t=double(0),
                 nu=double(1),
                 kappa=double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    if (n_dt<n_rt | n_dt>n_t) {
      if (log) return(-Inf)
      else return(0.0)
    }
    else{
    log_lik <- log_lik_Y(sens = nu,
                         spec = kappa,
                         test_mat = x) 
    log_prob <- log_lik[n_d-n_r+1]
    if (log) return(log_prob)
    else return(exp(log_prob))
    }
  },
)

registerDistributions(list(
  dTests = list(
    BUGSdist = "dTests(n_dt,n_rt,n_t,nu,kappa)",
    types = c("value = double(2)","n_dt = double(0)","n_rt = double(0)", "n_t = double(0)",
              "nu = double(1)","kappa = double(1)"),
    discrete = TRUE
  )
))

#Write model code

m_code <- nimbleCode({
  for(i in 1:Not_H_t){
    phi_binom_ns[i] ~ dnorm(0,sd=sigma_phi_binom)
    pi_ns[i] <- expit(alpha + phi_binom_ns[i])
  }
  
  for(r in 1:R){
    y[r] ~ dpois(pop[r]*exp(lambda[r])*sum(prob_pres_ns[r,]*pi_ns[1:Not_H_t]))
    lambda[r] <- A + Phi_pois[r]
    Phi_pois[r] ~ dnorm(0,sd=sigma_Phi_pois)
  }
  
  for(h in 1:H_t){
    n_r[h] ~ dbinom(prob=pi_s[h],size=n_d[h])
    pi_s[h] <- expit(alpha + phi_binom_s[h])
    phi_binom_s[h] ~ dnorm(0,sd=sigma_phi_binom)
    n_d[h] ~ dpois(sum(pop[1:R]*exp(lambda[1:R])*prob_pres_s[,h]))
  }
  
  Y1[,] ~ dTests(log_lik=log_lik_test_1[1:length[1]], n_d=n_d[1], n_r=n_r[1], n_t=n_t[1])
  Y2[,] ~ dTests(log_lik=log_lik_test_2[1:length[2]], n_d=n_d[2], n_r=n_r[2], n_t=n_t[2])
  Y3[,] ~ dTests(log_lik=log_lik_test_3[1:length[3]], n_d=n_d[3], n_r=n_r[3], n_t=n_t[3])
  Y4[,] ~ dTests(log_lik=log_lik_test_4[1:length[4]], n_d=n_d[4], n_r=n_r[4], n_t=n_t[4])
  Y5[,] ~ dTests(log_lik=log_lik_test_5[1:length[5]], n_d=n_d[5], n_r=n_r[5], n_t=n_t[5])
  Y6[,] ~ dTests(log_lik=log_lik_test_6[1:length[6]], n_d=n_d[6], n_r=n_r[6], n_t=n_t[6])
  Y7[,] ~ dTests(log_lik=log_lik_test_7[1:length[7]], n_d=n_d[7], n_r=n_r[7], n_t=n_t[7])
  Y8[,] ~ dTests(log_lik=log_lik_test_8[1:length[8]], n_d=n_d[8], n_r=n_r[8], n_t=n_t[8])
  Y9[,] ~ dTests(log_lik=log_lik_test_9[1:length[9]], n_d=n_d[9], n_r=n_r[9], n_t=n_t[9])
  Y10[,] ~ dTests(log_lik=log_lik_test_10[1:length[10]], n_d=n_d[10], n_r=n_r[10], n_t=n_t[10])
  
  A ~ dnorm(0, sd = 10)
  alpha ~ dnorm(0, sd=10)
  sigma_Phi_pois ~ dunif(0, 10)
  sigma_phi_binom ~ dunif(0,10)
})

#Fit and run model

data_list <- list(y=y,n_r=n_r,Y1 = matrix_test[[1]],
                  Y2 = matrix_test[[2]],
                  Y3 = matrix_test[[3]],
                  Y4 = matrix_test[[4]],
                  Y5 = matrix_test[[5]],
                  Y6 = matrix_test[[6]],
                  Y7 = matrix_test[[7]],
                  Y8 = matrix_test[[8]],
                  Y9 = matrix_test[[9]],
                  Y10 = matrix_test[[10]],
                  pop=region_level_data$population/100000,
                  log_lik_test_1=log_lik_test_1,
                  log_lik_test_2=log_lik_test_2,
                  log_lik_test_3=log_lik_test_3,
                  log_lik_test_4=log_lik_test_4,
                  log_lik_test_5=log_lik_test_5,
                  log_lik_test_6=log_lik_test_6,
                  log_lik_test_7=log_lik_test_7,
                  log_lik_test_8=log_lik_test_8,
                  log_lik_test_9=log_lik_test_9,
                  log_lik_test_10=log_lik_test_10)

constant_list <- list(R=nrow(region_cases),
                      prob_pres_ns=as.matrix(pres_probs_region_non_surveilled_hospitals[,-1]),
                      prob_pres_s=as.matrix(pres_probs_region_surveilled_hospitals),
                      H_t=length(n_r),
                      Not_H_t=ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[,-1])),
                      n_t=n_t,
                      length=n_t-n_r+1)

D_nd_init_function <- function(S, n_t, n_r, i) {
  n_d <- rep(0, length(n_t))
  
  nt <- n_t[i]
  nr <- n_r[i]
  
  # nd must be at least nr and at most nt
  nd <- sample(nr:nt, 1)
  
  mat <- matrix(0, nrow = nt, ncol = S)
  
  for (s in 1:S) {
    vec <- rep(0, nt)
    
    if (nr > 0) {
      vec[1:nr] <- 1
    }
    
    n_extra <- nd - nr
    remaining_indices <- if (nr < nt) (nr + 1):nt else integer(0)
    
    if (n_extra > 0) {
      if (length(remaining_indices) < n_extra) {
        stop(sprintf("Not enough positions to place %d extra 1s for i=%d", n_extra, i))
      }
      
      rand_indices <- sample(remaining_indices, size = n_extra, replace = FALSE)
      vec[rand_indices] <- 1
    }
    
    mat[, s] <- vec
  }
  return(list(D = mat, n_d = nd))
}


# Create function to initialize Phi_pois and sigma_phi_pois

Phi_sigma_pois_init <- function(upper_pois) {
  sigma_Phi_pois <- runif(1, 0, upper_pois)
  Phi_pois <- rnorm(nrow(region_cases), 0, sigma_Phi_pois)
  return(list(sigma_Phi_pois = sigma_Phi_pois, Phi_pois = Phi_pois))
}

# Create function to initialize phi_binom_s, phi_binom_ns, and sigma_phi_binom

phi_sigma_binom_init <- function(n_r, upper_binom) {
  sigma_phi_binom <- runif(1, 0, upper_binom)
  phi_binom_s <- rnorm(length(n_r), 0, sigma_phi_binom)
  phi_binom_ns <- rnorm(ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[, -1])), 0, sigma_phi_binom)
  return(list(sigma_phi_binom = sigma_phi_binom, phi_binom_s = phi_binom_s, phi_binom_ns = phi_binom_ns))
}


# Create final inits function

inits_function <- function(S, n_t, n_r, upper_pois, upper_binom) {
  D_nd1 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 1)
  D_nd2 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 2)
  D_nd3 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 3)
  D_nd4 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 4)
  D_nd5 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 5)
  D_nd6 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 6)
  D_nd7 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 7)
  D_nd8 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 8)
  D_nd9 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 9)
  D_nd10 <- D_nd_init_function(S = S, n_t = n_t, n_r = n_r, i = 10)
  Phi_sigma_pois <- Phi_sigma_pois_init(upper_pois = upper_pois)
  phi_sigma_binom <- phi_sigma_binom_init(n_r = n_r, upper_binom = upper_binom)
  return(list(
    n_d = c(D_nd1$n_d,
            D_nd2$n_d,
            D_nd3$n_d,
            D_nd4$n_d,
            D_nd5$n_d,
            D_nd6$n_d,
            D_nd7$n_d,
            D_nd8$n_d,
            D_nd9$n_d,
            D_nd10$n_d),
    sigma_Phi_pois = Phi_sigma_pois$sigma_Phi_pois,
    Phi_pois = Phi_sigma_pois$Phi_pois,
    sigma_phi_binom = phi_sigma_binom$sigma_phi_binom,
    phi_binom_s = phi_sigma_binom$phi_binom_s,
    phi_binom_ns = phi_sigma_binom$phi_binom_ns,
    A = rnorm(1, 0, 10),
    alpha = rnorm(1, 0, 10)
  ))
}


# Initialize and configure model

Rmodel <- nimbleModel(
  m_code,
  constants = constant_list,
  data = data_list,
  inits = inits_function(S = 2, n_t = n_t, n_r = n_r, upper_pois = 4, upper_binom = 4),
  dimensions = list(pi_ns = ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[, -1]))),
  buildDerivs = TRUE
)

conf <- configureHMC(Rmodel)

## Then add a single joint_RW sampler
conf$addMonitors(c("alpha", "A", "sigma_Phi_pois", "sigma_phi_binom", "n_d"))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

start <- Sys.time()

samples <- runMCMC(Cmcmc,
                   nburnin = 2500,
                   niter = 5000,
                   nchains = 2,
                   inits = inits_function(S = 2, n_t = n_t, n_r = n_r, upper_pois = 4, upper_binom = 4),
                   samplesAsCodaMCMC = TRUE)

end <- Sys.time()
end - start

# Plot trace plots
plot(samples[, c("alpha", "A", "sigma_Phi_pois", "sigma_phi_binom")], density = F)
plot(samples[, c("n_d[1]", "n_d[2]", "n_d[3]", "n_d[4]")], density = F)
plot(samples[, c("n_d[5]", "n_d[6]", "n_d[7]", "n_d[8]")], density = F)
plot(samples[, c("n_d[9]", "n_d[10]")], density = F)

# Create summary
summary(samples[, c("alpha", "A", "sigma_Phi_pois", "sigma_phi_binom")])
summary(samples[, c("n_d[1]", "n_d[2]", "n_d[3]", "n_d[4]",
                    "n_d[5]", "n_d[6]", "n_d[7]", "n_d[8]",
                    "n_d[9]", "n_d[10]")])
