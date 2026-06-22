library(tidyverse)
library(nimble)
library(nimbleHMC)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))
source(data_path("scripts", "Simulation_Code.R"))

#Get hospital presentation probabilities

pres_probs_region<-as.data.frame(read_csv(data_path("processed", "simulation_inputs", "health_regions_probs.csv")))

#Get population in regions

region_level_data<-as.data.frame(read_csv(data_path("processed", "simulation_inputs", "region_population.csv")))

region_level_data$lepto_cases_2023 <- c(42,36,80,4,47,26,47)

#Get data for each 1 km x 1 km cells

cell_data<-as.data.frame(read_csv(data_path("processed", "simulation_inputs", "cell_probabilities.csv")))


#Create rates per 100,000

region_level_data <- region_level_data %>% mutate(lepto_per_100k = lepto_cases_2023*100000/population)

mean_log_rate <- mean(log(28/12*region_level_data$lepto_per_100k))
sd_log_rate <- sd(log(28/12*region_level_data$lepto_per_100k))

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

mean_logit_probability <- logit(12/28)
sd_logit_probability <- 0.5

#Create probabilities

pi <- inv_logit(rnorm(ncol(pres_probs_region)-1,mean_logit_probability,sd_logit_probability))
pi_df <- data.frame(hospital=colnames(pres_probs_region)[-1],pi=pi)

link_df <- tidyr::crossing(region=region_level_data$region,hospital=pi_df$hospital) %>% 
  left_join(pi_df)

#Create data set

data <- data_generation(lambda=lambda_full,pi=link_df$pi,multiplier_lower=4,multiplier_upper=8,
                        sensitivity=c(0.6,0.9),specificity=c(0.9,0.9),test_receipt_probs=c(1,0.6))

#Surveil a set of hospitals

# Set the proportion to sample
p <- 0.15625

surveillance <- data$hospital_data %>%
  filter(hospital %in% sample(unique(hospital), size = ceiling(p * n_distinct(hospital)))) %>%
  group_by(hospital) %>%
  arrange(desc(passive_surveillance_capture), .by_group=T) %>%
  ungroup()

split_df <- split(surveillance[, 4:6], surveillance$hospital)
mat_conv <- function(matrix){
  matrix <- as.matrix(matrix)
  return(matrix)
}


# Step 2: Convert each subset to a matrix
matrix_test <- lapply(split_df, mat_conv)


#Create a function to calculate likelihood

likelihood_nD_log <- function(sens, spec, test_mat) {
  # sens: vector of sensitivities (nu_k)
  # spec: vector of specificities (kappa_k)
  # test_mat: matrix, first col = passive surveillance indicator,
  #           remaining cols = test results (0/1 or -1 for missing)
  
  n_total <- nrow(test_mat)
  K <- length(sens)
  
  # Split into R (captured passively) and T\R (to handle with recursion)
  passive_idx <- which(test_mat[, 1] == 1)
  active_idx  <- which(test_mat[, 1] == 0 | test_mat[, 1] == -1 )
  
  nR <- length(passive_idx)
  nRemain <- length(active_idx)
  
  Y_passive <- test_mat[passive_idx, -1, drop = FALSE]
  Y_active  <- test_mat[active_idx, -1, drop = FALSE]
  
  ## Contribution from passively observed positives (on log scale)
  log_term_passive <- 0
  if (nR > 0) {
    for (j in 1:nR) {
      for (k in 1:K) {
        y <- Y_passive[j, k]
        if (y == -1) next
        log_term_passive <- log_term_passive + 
          (y * log(sens[k]) + (1 - y) * log(1 - sens[k]))
      }
    }
  }
  
  ## Compute log(w_j^1) and log(w_j^0) for the recursive part
  log_w1 <- numeric(nRemain)
  log_w0 <- numeric(nRemain)
  
  for (j in 1:nRemain) {
    lw1 <- 0
    lw0 <- 0
    for (k in 1:K) {
      y <- Y_active[j, k]
      if (y == -1) next
      lw1 <- lw1 + (y * log(sens[k]) + (1 - y) * log(1 - sens[k]))
      lw0 <- lw0 + (y * log(1 - spec[k]) + (1 - y) * log(spec[k]))
    }
    log_w1[j] <- lw1
    log_w0[j] <- lw0
  }
  
  ## Recursion to compute log xi_{m,l}
  # Initialize with -Inf (log 0)
  log_xi <- matrix(-Inf, nrow = nRemain + 1, ncol = nRemain + 1)
  log_xi[1, 1] <- 0  # log(1) = 0, corresponds to xi_{0,0}
  
  logsumexp <- function(a, b) {
    # helper for log(exp(a)+exp(b)) avoiding underflow
    if (is.infinite(a) && is.infinite(b)) return(-Inf)
    m <- max(a, b)
    return(m + log(exp(a - m) + exp(b - m)))
  }
  
  for (l in 1:nRemain) {
    for (m in 0:l) {
      vals <- c()
      # disease = 0
      if (!is.infinite(log_xi[m + 1, l])) {
        vals <- c(vals, log_w0[l] + log_xi[m + 1, l])
      }
      # disease = 1
      if (m > 0 && !is.infinite(log_xi[m, l])) {
        vals <- c(vals, log_w1[l] + log_xi[m, l])
      }
      if (length(vals) > 0) {
        log_xi[m + 1, l + 1] <- Reduce(logsumexp, vals)
      }
    }
  }
  
  ## Extract log-likelihoods for nD = total # diseased
  log_likelihood <- log_xi[, nRemain + 1] + log_term_passive - lchoose(nRemain, 0:nRemain)
  names(log_likelihood) <- nR:(nR + nRemain)
  
  return(log_likelihood)
}

sens <- c(0.6, 0.9)
spec <- c(0.9, 0.9)

for (i in seq_along(matrix_test)) {
  res <- likelihood_nD_log(sens, spec, matrix_test[[i]])
  assign(paste0("log_lik_test_", i), res)
}


#Get vector of cases at surveilled hospitals and total number tested

cases_passive<-surveillance %>% group_by(hospital) %>% 
  mutate(passive_surveillance_capture_mod=ifelse(passive_surveillance_capture!=-1,
                                                 passive_surveillance_capture,0)) %>%
  summarize(n_c=sum(passive_surveillance_capture_mod))
n_c <- cases_passive$n_c

cases_total<-surveillance %>% group_by(hospital) %>% 
  summarize(n_a=n())
n_a <- cases_total$n_a

#Get vector of true number of cases presenting for comparison

cases_true<-surveillance %>% group_by(hospital) %>% 
  summarize(n_p=sum(disease_status))
n_p <- cases_true$n_p


#Create data for the input into model

passive_cases_surveilled_hospitals<-surveillance %>%
  mutate(passive_surveillance_capture=ifelse(passive_surveillance_capture==-1,0,passive_surveillance_capture)) %>%
  group_by(health_region) %>%
  summarize(y=sum(passive_surveillance_capture))

region_cases <- data$region_level_data
y <- region_cases$y - passive_cases_surveilled_hospitals$y

#Remove columns for processing

#Presentation probabilities for non-surveilled hospitals

pres_probs_region_non_surveilled_hospitals <- pres_probs_region %>%
  dplyr::select(-all_of(surveillance$hospital))

#Presentation probabilities for surveilled hospitals

pres_probs_region_surveilled_hospitals <- pres_probs_region %>%
  dplyr::select(all_of(surveillance$hospital))


#Create distribution for vector of test results


dTests <- nimbleFunction(
  run = function(x = double(2),
                 log_lik = double(1),
                 n_p=double(0),
                 n_c=double(0),
                 n_a=double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    if (n_p<n_c | n_p>n_a) {
      if (log) return(-Inf)
      else return(0.0)
    }
    else{
    log_prob <- log_lik[n_p-n_c+1]
    if (log) return(log_prob)
    else return(exp(log_prob))
    }
  },
)

registerDistributions(list(
  dTests = list(
    BUGSdist = "dTests(log_lik,n_p,n_c,n_a)",
    types = c("value = double(2)","log_lik = double(1)","n_p = double(0)","n_c = double(0)", "n_a = double(0)"),
    discrete = TRUE
  )
))

#Write model code

m_code <- nimbleCode({
  #Cases at region level
  
  for(i in 1:Not_H_t){
    phi_binom_ns[i] ~ dnorm(0,sd=sigma_phi_binom)
    pi_ns[i] <- expit(alpha + phi_binom_ns[i])
  }
  
  
  for(r in 1:R){
    y[r] ~ dpois(pop[r]*exp(lambda[r])*sum(prob_pres_ns[r,]*pi_ns[1:Not_H_t]))
    lambda[r] <- A + Phi_pois[r]
    Phi_pois[r] ~ dnorm(0,sd=sigma_Phi_pois)
  }
  
  #Cases at the hospital level
  
  for(h in 1:H_t){
    n_c[h] ~ dbinom(prob=pi_s[h],size=n_p[h])
    pi_s[h] <- expit(alpha + phi_binom_s[h])
    phi_binom_s[h] ~ dnorm(0,sd=sigma_phi_binom)
    n_p[h] ~ dpois(sum(pop[1:R]*exp(lambda[1:R])*prob_pres_s[,h]))
  }
  
  #Testing data
  
  Y1[,] ~ dTests(log_lik=log_lik_test_1[1:length[1]], n_p=n_p[1], n_c=n_c[1], n_a=n_a[1])
  Y2[,] ~ dTests(log_lik=log_lik_test_2[1:length[2]], n_p=n_p[2], n_c=n_c[2], n_a=n_a[2])
  Y3[,] ~ dTests(log_lik=log_lik_test_3[1:length[3]], n_p=n_p[3], n_c=n_c[3], n_a=n_a[3])
  Y4[,] ~ dTests(log_lik=log_lik_test_4[1:length[4]], n_p=n_p[4], n_c=n_c[4], n_a=n_a[4])
  Y5[,] ~ dTests(log_lik=log_lik_test_5[1:length[5]], n_p=n_p[5], n_c=n_c[5], n_a=n_a[5])
  Y6[,] ~ dTests(log_lik=log_lik_test_6[1:length[6]], n_p=n_p[6], n_c=n_c[6], n_a=n_a[6])
  Y7[,] ~ dTests(log_lik=log_lik_test_7[1:length[7]], n_p=n_p[7], n_c=n_c[7], n_a=n_a[7])
  Y8[,] ~ dTests(log_lik=log_lik_test_8[1:length[8]], n_p=n_p[8], n_c=n_c[8], n_a=n_a[8])
  Y9[,] ~ dTests(log_lik=log_lik_test_9[1:length[9]], n_p=n_p[9], n_c=n_c[9], n_a=n_a[9])
  Y10[,] ~ dTests(log_lik=log_lik_test_10[1:length[10]], n_p=n_p[10], n_c=n_c[10], n_a=n_a[10])
  
  #Priors
  
  A ~ dnorm(0, sd = 10)
  alpha ~ dnorm(0, sd=10)
  sigma_Phi_pois ~ dunif(0, 10)
  sigma_phi_binom ~ dunif(0,10)
  
})


#Fit and run model

data_list <- list(y=y,n_c=n_c,Y1 = matrix_test[[1]],
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
                      H_t=length(n_c),
                      Not_H_t=ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[,-1])),
                      n_a=n_a,
                      length=n_a-n_c+1)

#Create initialization functions for parameters

D_np_init_function <- function(S, n_a, n_c,i) {
  n_p <- rep(0, length(n_a))
  
  na <- n_a[i]
  nc <- n_c[i]
  
  # np must be at least nc and at most na
  np <- sample(nc:na,1)
  
  mat <- matrix(0, nrow = na, ncol = S)
  
  for (s in 1:S) {
    vec <- rep(0, na)
    
    if (nc > 0) {
      vec[1:nc] <- 1
    }
    
    n_extra <- np - nc
    remaining_indices <- if (nc < na) (nc + 1):na else integer(0)
    
    if (n_extra > 0) {
      if (length(remaining_indices) < n_extra) {
        stop(sprintf("Not enough positions to place %d extra 1s for i=%d", n_extra, i))
      }
      
      rand_indices <- sample(remaining_indices, size = n_extra, replace = FALSE)
      vec[rand_indices] <- 1
    }
    
    mat[, s] <- vec
  }
  return(list(D=mat,n_p=np))
}


#Create function to initialize Phi_pois and sigma_phi_pois

Phi_sigma_pois_init<-function(upper_pois){
  sigma_Phi_pois <- runif(1,0,upper_pois)
  Phi_pois <- rnorm(nrow(region_cases),0,sigma_Phi_pois)
  return(list(sigma_Phi_pois=sigma_Phi_pois,Phi_pois=Phi_pois))
}

#Create function to initialize phi_binom_s, phi_binom_ns, and sigma_phi_binom

phi_sigma_binom_init<-function(n_c,upper_binom){
  sigma_phi_binom <- runif(1,0,upper_binom)
  phi_binom_s <- rnorm(length(n_c),0,sigma_phi_binom)
  phi_binom_ns <- rnorm(ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[,-1])),0,sigma_phi_binom)
  return(list(sigma_phi_binom=sigma_phi_binom,phi_binom_s=phi_binom_s,phi_binom_ns=phi_binom_ns))
}


#Create final inits function

inits_function <- function(S,n_a,n_c,upper_pois,upper_binom){
  D_np1 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=1)
  D_np2 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=2)
  D_np3 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=3)
  D_np4 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=4)
  D_np5 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=5)
  D_np6 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=6)
  D_np7 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=7)
  D_np8 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=8)
  D_np9 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=9)
  D_np10 <- D_np_init_function(S=S, n_a=n_a, n_c=n_c, i=10)
  Phi_sigma_pois <- Phi_sigma_pois_init(upper_pois=upper_pois)
  phi_sigma_binom <- phi_sigma_binom_init(n_c=n_c,upper_binom=upper_binom)
  return(list(n_p=c(D_np1$n_p,
                     D_np2$n_p,
                     D_np3$n_p,
                     D_np4$n_p,
                     D_np5$n_p,
                     D_np6$n_p,
                     D_np7$n_p,
                     D_np8$n_p,
                     D_np9$n_p,
                     D_np10$n_p),sigma_Phi_pois=Phi_sigma_pois$sigma_Phi_pois,
              Phi_pois=Phi_sigma_pois$Phi_pois,sigma_phi_binom=phi_sigma_binom$sigma_phi_binom,
              phi_binom_s=phi_sigma_binom$phi_binom_s,phi_binom_ns=phi_sigma_binom$phi_binom_ns,
              A=rnorm(1,0,10),alpha=rnorm(1,0,10)))
}


#Initialize and configure model

Rmodel <- nimbleModel(
  m_code,
  constants = constant_list,
  data = data_list,
  inits=inits_function(S=2,n_a=n_a,n_c=n_c,upper_pois=4,upper_binom=4),
  dimensions = list(pi_ns=ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[,-1]))),
  buildDerivs = TRUE)


conf <- configureHMC(Rmodel)


## Then add a single joint_RW sampler

conf$addMonitors(c("alpha", "A", "sigma_Phi_pois","sigma_phi_binom","n_p"))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel,showCompilerOutput = TRUE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


start <- Sys.time()

samples <- runMCMC(Cmcmc,
                   nburnin = 2500,
                   niter=5000,
                   nchains = 2,
                   inits = inits_function(S=2,n_a=n_a,n_c=n_c,upper_pois=4,upper_binom=4),
                   samplesAsCodaMCMC = TRUE)

end <- Sys.time()
end - start 



#Plot trace plots

plot(samples[,c("alpha", "A", "sigma_Phi_pois","sigma_phi_binom")], density = F)
plot(samples[,c("n_p[1]",'n_p[2]','n_p[3]','n_p[4]')], density = F)
plot(samples[,c("n_p[5]",'n_p[6]','n_p[7]','n_p[8]')], density = F)
plot(samples[,c("n_p[9]",'n_p[10]')], density = F)

#Create summary

summary(samples[,c("alpha", "A", "sigma_Phi_pois","sigma_phi_binom")])
summary(samples[,c("n_p[1]",'n_p[2]','n_p[3]','n_p[4]',"n_p[5]",'n_p[6]','n_p[7]','n_p[8]',
                   'n_p[9]','n_p[10]')])
