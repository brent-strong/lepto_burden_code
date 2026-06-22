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

# Sample x% of unique hospitals and filter the rows
surveillance <- data$hospital_data %>%
  filter(hospital %in% sample(unique(hospital), size = ceiling(p * n_distinct(hospital)))) %>%
  group_by(hospital) %>%
  arrange(desc(passive_surveillance_capture), .by_group=T)

# Step 1: Split the data frame by hospital
split_df <- split(surveillance[, 5:6], surveillance$hospital)
mat_conv <- function(matrix){
  matrix <- as.matrix(matrix)
  colnames(matrix)<-NULL
  return(matrix)
}

# Step 2: Convert each subset to a matrix
matrix_test <- lapply(split_df, mat_conv)


surveillance <- data$hospital_data %>%
  filter(hospital %in% sample(unique(hospital), size = ceiling(p * n_distinct(hospital)))) %>%
  group_by(hospital) %>%
  arrange(desc(passive_surveillance_capture), .by_group=T)

# Step 1: Split the data frame by hospital
split_df <- split(surveillance[, 3:6], surveillance$hospital)
mat_conv <- function(matrix){
  matrix <- as.matrix(matrix)
  colnames(matrix)<-NULL
  return(matrix)
}

# Step 2: Convert each subset to a matrix
matrix_test <- lapply(split_df, as.matrix)

for(i in 1:length(matrix_test)){
  colnames(matrix_test[[i]])<-NULL
}


#Get vector of cases at surveilled hospitals and total number tested

cases_passive<-surveillance %>% group_by(hospital) %>% 
  mutate(passive_surveillance_capture_mod=ifelse(passive_surveillance_capture!=-1,
                                                 passive_surveillance_capture,0)) %>%
  summarize(n_r=sum(passive_surveillance_capture_mod))
n_r<- cases_passive$n_r

cases_total<-surveillance %>% group_by(hospital) %>% 
  summarize(n_t=n())
n_t <- cases_total$n_t

#Get vector of true number of cases presenting for comparison

cases_true<-surveillance %>% group_by(hospital) %>% 
  summarize(n_d=sum(disease_status))
n_d <- cases_true$n_d


#Create data for the input into model

passive_cases_surveilled_hospitals<-surveillance %>% 
  group_by(health_region) %>% 
  filter(passive_surveillance_capture!=-1) %>%
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


#Create distribution for array D

#Define the density function as a nimbleFunction

dFixedSum <- nimbleFunction(
  run = function(x = double(2), n_d = double(0),n_r = double(0), 
                 n_t = double(0),B=double(0),log = integer(0, default = 0)) {
    returnType(double(0))
    flag <- rep(0.0,B)
    for(b in 1:B){
      sum_x <- sum(x[,b])
      if (sum_x != n_d | sum_x < n_r | sum_x > n_t) {
        flag[b]<-1
      } 
    }
    if (sum(flag) != 0) {
      if (log) return(-Inf)
      else return(0.0)
    }
    else{
      log_prob <- B*(lfactorial(n_d-n_r) + lfactorial(n_t-n_r) - lfactorial(n_d-n_r)) 
      if(log) return(log_prob)
      else return(exp(log_prob))
    }
  })


#Register the custom distribution 

registerDistributions(list(
  dFixedSum = list(
    BUGSdist = "dFixedSum(n_d,n_r,n_t,B)",
    types = c("value = double(2)",    # this tells NIMBLE that x is a vector
              "n_d = double(0)",
              "n_r = double(0)",
              "n_t =double(0)",
              "B=double(0)"),
    discrete = TRUE
  )
))

#Create distribution for vector of test results


dTests <- nimbleFunction(
  run = function(x = double(2), nu = double(1), kappa = double(1),
                 status = double(2), B = double(0), K = double(0),
                 n_t = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    log_prob_sample <- rep(0.0, B)
    for (b in 1:B) {
      log_prob_individual <- rep(0.0, n_t)
      for (i in 1:n_t) {
        log_prob_individual_test <- rep(0.0,K)
        for(k in 1:K){
          if (x[i,k] == -1) {
            log_prob_individual_test[k] <- 0.0
          } else {
            if (status[i, b] == 1) {
              if (x[i,k] == 1)
                log_prob_individual_test[k] <- log(nu[k])
              else
                log_prob_individual_test[k] <- log(1 - nu[k])
            } else {
              if (x[i,k] == 1)
                log_prob_individual_test[k] <- log(1 - kappa[k])
              else
                log_prob_individual_test[k] <- log(kappa[k])
            }
          }
        }
        log_prob_individual[i] <- sum(log_prob_individual_test[1:K])
      }
      log_prob_sample[b] <- sum(log_prob_individual[1:n_t])
    }
    log_prob <- log(1/B) + max(log_prob_sample[1:B]) + log(sum(exp(log_prob_sample[1:B]-max(log_prob_sample[1:B]))))
    if (log) return(log_prob)
    else return(exp(log_prob))
  },
)



registerDistributions(list(
  dTests = list(
    BUGSdist = "dTests(nu,kappa,status,B,K,n_t)",
    types = c("value = double(2)", "nu = double(1)", "kappa = double(1)", "status = double(2)", "B = double(0)",
              "K = double(0)", "n_t = double(0)"),
    discrete = TRUE
  )
))

#Create sampler for array D

joint_RW <- nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    B <- control$B
    n_t <- control$n_t
    n_r <- control$n_r
  },
  
  run = function() {
    # initial model logProb
    log_MH_ratio <- 0.0
    model_lp_proposed <- 0.0
    model_lp_initial <- 0.0
    n_d <- 0.0
    proposal_n_d <- 0.0
    model_lp_initial <- model$getLogProb(calcNodes)
    u <- runif(1, 0, 1)
    n_d <- values(model, target[1])[1]
    D <- values(model, target[2])
    
    if(u>=0.5 & n_d<n_t){
      proposal_n_d <- n_d+1
      proposal_D_vec <- values(model, target[2])
      proposal_D <- matrix(proposal_D_vec, nrow = n_t, ncol = B)
      for(b in 1:B){
        ind <- which(proposal_D[,b]==0)
        j <- ceiling(runif(1,0,length(ind)))
        proposal_D[ind[j],b]<-1
        for(i in (n_r+1):(n_t-1)){
          j <- ceiling(runif(1,i-1,n_t))
          temp_i <- proposal_D[i,b]
          temp_j <- proposal_D[j,b]
          proposal_D[i,b] <- temp_j
          proposal_D[j,b] <- temp_i  
        }
      }
      
      # generate proposal
      # store proposal into model
      values(model, target[1]) <<- c(proposal_n_d)
      values(model, target[2]) <<- c(proposal_D)
      
      # proposal model logProb
      model_lp_proposed <- model$calculate(calcNodes)
      
      # log-Metropolis-Hastings ratio
      log_MH_ratio <- model_lp_proposed + B*(lfactorial(n_t-n_r) - 
                                               lfactorial(proposal_n_d-n_r) - lfactorial(n_t-proposal_n_d)) -
        model_lp_initial - B*(lfactorial(n_t-n_r) - 
                                lfactorial(n_d-n_r) - lfactorial(n_t-n_d))
    }
    
    
    if(u>=0.5 & n_d==n_t){
      proposal_n_d <- n_d
      proposal_D_vec <- values(model, target[2])
      proposal_D <- matrix(proposal_D_vec, nrow = n_t, ncol = B)
      for(b in 1:B){
        for(i in (n_r+1):(n_t-1)){
          j <- ceiling(runif(1,i-1,n_t))
          temp_i <- proposal_D[i,b]
          temp_j <- proposal_D[j,b]
          proposal_D[i,b] <- temp_j
          proposal_D[j,b] <- temp_i  
        }
      }
      
      # generate proposal
      # store proposal into model
      values(model, target[1]) <<- c(proposal_n_d)
      values(model, target[2]) <<- c(proposal_D)
      
      # proposal model logProb
      model_lp_proposed <- model$calculate(calcNodes)
      
      # log-Metropolis-Hastings ratio
      log_MH_ratio <- model_lp_proposed  - model_lp_initial
    }
    
    
    if(u<0.5 & n_d==n_r){
      proposal_n_d <- n_d
      proposal_D_vec <- values(model, target[2])
      proposal_D <- matrix(proposal_D_vec, nrow = n_t, ncol = B)
      for(b in 1:B){
        for(i in (n_r+1):(n_t-1)){
          j <- ceiling(runif(1,i-1,n_t))
          temp_i <- proposal_D[i,b]
          temp_j <- proposal_D[j,b]
          proposal_D[i,b] <- temp_j
          proposal_D[j,b] <- temp_i  
        }
      }
      # generate proposal
      # store proposal into model
      values(model, target[1]) <<- c(proposal_n_d)
      values(model, target[2]) <<- c(proposal_D)
      
      # proposal model logProb
      model_lp_proposed <- model$calculate(calcNodes)
      
      # log-Metropolis-Hastings ratio
      log_MH_ratio <- model_lp_proposed - model_lp_initial
    }
    
    
    if(u<0.5 & n_d>n_r){
      proposal_n_d <- n_d-1
      proposal_D_vec <- values(model, target[2])
      proposal_D <- matrix(proposal_D_vec, nrow = n_t, ncol = B)
      for(b in 1:B){
        ind <- which(proposal_D[,b]==1)
        j <- ceiling(runif(1,0,length(ind)))
        proposal_D[ind[j],b]<-0
        for(i in (n_r+1):(n_t-1)){
          j <- ceiling(runif(1,i-1,n_t))
          temp_i <- proposal_D[i,b]
          temp_j <- proposal_D[j,b]
          proposal_D[i,b] <- temp_j
          proposal_D[j,b] <- temp_i  
        }
      }
      # generate proposal
      # store proposal into model
      values(model, target[1]) <<- c(proposal_n_d)
      values(model, target[2]) <<- c(proposal_D)
      
      # proposal model logProb
      model_lp_proposed <- model$calculate(calcNodes)
      
      # log-Metropolis-Hastings ratio
      log_MH_ratio <- model_lp_proposed + B*(lfactorial(n_t-n_r) - 
                                          lfactorial(proposal_n_d-n_r) - lfactorial(n_t-proposal_n_d)) -
                                          model_lp_initial - B*(lfactorial(n_t-n_r) - 
                                          lfactorial(n_d-n_r) - lfactorial(n_t-n_d))
}
      # Metropolis-Hastings step: determine whether or
      # not to accept the newly proposed value
      u_alpha <- runif(1, 0, 1)
      if(u_alpha < exp(log_MH_ratio)) jump <- TRUE
      else                      jump <- FALSE
      
      # keep the model and mvSaved objects consistent
      if(jump) copy(from = model, to = mvSaved, row = 1, 
                    nodes = calcNodes, logProb = TRUE)
      else     copy(from = mvSaved, to = model, row = 1,
                    nodes = calcNodes, logProb = TRUE)
    },
    
    methods = list(   reset = function () {}   )
)


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
    n_r[h] ~ dbinom(prob=pi_s[h],size=n_d[h])
    pi_s[h] <- expit(alpha + phi_binom_s[h])
    phi_binom_s[h] ~ dnorm(0,sd=sigma_phi_binom)
    n_d[h] ~ dpois(sum(pop[1:R]*exp(lambda[1:R])*prob_pres_s[,h]))
  }
  #Testing data
  
  
  Y1[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D1[1:n_t[1],1:B],B=B,K=K,n_t=n_t[1])
  Y2[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D2[1:n_t[2],1:B],B=B,K=K,n_t=n_t[2])
  Y3[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D3[1:n_t[3],1:B],B=B,K=K,n_t=n_t[3])
  Y4[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D4[1:n_t[4],1:B],B=B,K=K,n_t=n_t[4])
  Y5[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D5[1:n_t[5],1:B],B=B,K=K,n_t=n_t[5])
  Y6[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D6[1:n_t[6],1:B],B=B,K=K,n_t=n_t[6])
  Y7[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D7[1:n_t[7],1:B],B=B,K=K,n_t=n_t[7])
  Y8[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D8[1:n_t[8],1:B],B=B,K=K,n_t=n_t[8])
  Y9[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D9[1:n_t[9],1:B],B=B,K=K,n_t=n_t[9])
  Y10[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D10[1:n_t[10],1:B],B=B,K=K,n_t=n_t[10])
  
  
  #Priors
  A ~ dnorm(0, sd = 10)
  alpha ~ dnorm(0, sd=10)
  sigma_Phi_pois ~ dunif(0, 10)
  sigma_phi_binom ~ dunif(0,10)
  nu[1] ~ dbeta(6000,4000)
  nu[2] ~ dbeta(9000,1000)
  kappa[1] ~ dbeta(9000,1000)
  kappa[2] ~ dbeta(9000,1000)
  
  
  #Distribution on array of disease status
  D1[1:n_t[1],1:B] ~ dFixedSum(n_d=n_d[1], n_r=n_r[1],n_t=n_t[1],B=B)
  D2[1:n_t[2],1:B] ~ dFixedSum(n_d=n_d[2], n_r=n_r[2],n_t=n_t[2],B=B)
  D3[1:n_t[3],1:B] ~ dFixedSum(n_d=n_d[3], n_r=n_r[3],n_t=n_t[3],B=B)
  D4[1:n_t[4],1:B] ~ dFixedSum(n_d=n_d[4], n_r=n_r[4],n_t=n_t[4],B=B)
  D5[1:n_t[5],1:B] ~ dFixedSum(n_d=n_d[5], n_r=n_r[5],n_t=n_t[5],B=B)
  D6[1:n_t[6],1:B] ~ dFixedSum(n_d=n_d[6], n_r=n_r[6],n_t=n_t[6],B=B)
  D7[1:n_t[7],1:B] ~ dFixedSum(n_d=n_d[7], n_r=n_r[7],n_t=n_t[7],B=B)
  D8[1:n_t[8],1:B] ~ dFixedSum(n_d=n_d[8], n_r=n_r[8],n_t=n_t[8],B=B)
  D9[1:n_t[9],1:B] ~ dFixedSum(n_d=n_d[9], n_r=n_r[9],n_t=n_t[9],B=B)
  D10[1:n_t[10],1:B] ~ dFixedSum(n_d=n_d[10], n_r=n_r[10],n_t=n_t[10],B=B)
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
                  pop=region_level_data$population/100000)


constant_list <- list(R=nrow(region_cases),
                      prob_pres_ns=as.matrix(pres_probs_region_non_surveilled_hospitals[,-1]),
                      prob_pres_s=as.matrix(pres_probs_region_surveilled_hospitals),
                      H_t=length(n_r),
                      Not_H_t=ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[,-1])),
                      B=200,
                      K=2,
                      n_t=n_t)


#Create initialization functions for parameters


D_nd_init_function <- function(B, n_t, n_r,i) {
  n_d <- rep(0, length(n_t))
  nt <- n_t[i]
  nr <- n_r[i]
  # nd must be at least nr and at most nt
  nd <- sample(nr:nt,1)
  mat <- matrix(0, nrow = nt, ncol = B)
  for (b in 1:B) {
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
    mat[, b] <- vec
  }
  return(list(D=mat,n_d=nd))
}

#Create function to initialize Phi_pois and sigma_phi_pois


Phi_sigma_pois_init<-function(upper_pois){
  sigma_Phi_pois <- runif(1,0,upper_pois)
  Phi_pois <- rnorm(nrow(region_cases),0,sigma_Phi_pois)
  return(list(sigma_Phi_pois=sigma_Phi_pois,Phi_pois=Phi_pois))
}
  
  
#Create function to initialize phi_binom_s, phi_binom_ns, and sigma_phi_binom
  
  
phi_sigma_binom_init<-function(n_r,upper_binom){
    sigma_phi_binom <- runif(1,0,upper_binom)
    phi_binom_s <- rnorm(length(n_r),0,sigma_phi_binom)
    phi_binom_ns <- rnorm(ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[,-1])),0,sigma_phi_binom)
    return(list(sigma_phi_binom=sigma_phi_binom,phi_binom_s=phi_binom_s,phi_binom_ns=phi_binom_ns))
  }
  
  
D1<-D_nd_init_function(B=200, n_t=n_t[1], n_r=n_r[1], i=1)
  
  
  #Create final inits function
  
  
  inits_function <- function(B,n_t,n_r,upper_pois,upper_binom){
    D_nd1 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=1)
    D_nd2 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=2)
    D_nd3 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=3)
    D_nd4 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=4)
    D_nd5 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=5)
    D_nd6 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=6)
    D_nd7 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=7)
    D_nd8 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=8)
    D_nd9 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=9)
    D_nd10 <- D_nd_init_function(B=B, n_t=n_t, n_r=n_r, i=10)
    Phi_sigma_pois <- Phi_sigma_pois_init(upper_pois=upper_pois)
    phi_sigma_binom <- phi_sigma_binom_init(n_r=n_r,upper_binom=upper_binom)
    return(list(D1 = D_nd1$D,
                D2 = D_nd2$D,
                D3 = D_nd3$D,
                D4 = D_nd4$D,
                D5 = D_nd5$D,
                D6 = D_nd6$D,
                D7 = D_nd7$D,
                D8 = D_nd8$D,
                D9 = D_nd9$D,
                D10 = D_nd10$D
                ,n_d=c(D_nd1$n_d,
                       D_nd2$n_d,
                       D_nd3$n_d,
                       D_nd4$n_d,
                       D_nd5$n_d,
                       D_nd6$n_d,
                       D_nd7$n_d,
                       D_nd8$n_d,
                       D_nd9$n_d,
                       D_nd10$n_d),sigma_Phi_pois=Phi_sigma_pois$sigma_Phi_pois,
                Phi_pois=Phi_sigma_pois$Phi_pois,sigma_phi_binom=phi_sigma_binom$sigma_phi_binom,
                phi_binom_s=phi_sigma_binom$phi_binom_s,phi_binom_ns=phi_sigma_binom$phi_binom_ns,
                A=rnorm(1,0,10),alpha=rnorm(1,0,10),nu=c(0.6,0.9),
                kappa=c(0.9,0.9)))
  }


Rmodel <- nimbleModel(
  m_code,
  constants = constant_list,
  data = data_list,
  dimensions = list(pi_ns=ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[,-1]))),
  buildDerivs = TRUE)


conf <- configureHMC(Rmodel)

conf$replaceSamplers(target = "nu[1]",
                     type = 'RW')
conf$replaceSamplers(target = "nu[2]",
                     type = 'RW')
conf$replaceSamplers(target = "kappa[1]",
                     type = 'RW')
conf$replaceSamplers(target = "kappa[2]",
                     type = 'RW')
conf$replaceSamplers(target = conf$getUnsampledNodes(),
                     type = 'NUTS')


#Assign joint_RW sampler to n_d and Ds


for (i in 1:10) {
  conf$replaceSamplers(target = c(sprintf('n_d[%d]', i),
                                  sprintf('D%d', i)),
                       type = 'joint_RW',
                       control = list(B = 200,
                                      n_t = n_t[i],
                                      n_r= n_r[i]))
}


conf$addMonitors(c("alpha", "A", "sigma_Phi_pois","sigma_phi_binom","n_d","nu","kappa"))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel,showCompilerOutput = T)


samples <- runMCMC(Cmcmc,
                   nburnin = 0,
                   niter=1000,
                   nchains = 3,
                   inits = inits_function(B=200,n_t=n_t,n_r=n_r,upper_pois=4,upper_binom=4),
                   samplesAsCodaMCMC = TRUE)


#Plot trace plots


plot(samples[,c("alpha", "A", "sigma_Phi_pois","sigma_phi_binom")], density = F)
plot(samples[,c("n_d[1]",'n_d[2]','n_d[3]','n_d[4]')], density = F)
plot(samples[,c("n_d[5]",'n_d[6]','n_d[7]','n_d[8]')], density = F)
plot(samples[,c("n_d[9]",'n_d[10]')], density = F)


#Create summary


summary(samples[,c("alpha", "A", "sigma_Phi_pois","sigma_phi_binom")])
summary(samples[,c("n_d[1]",'n_d[2]','n_d[3]','n_d[4]',"n_d[5]",'n_d[6]','n_d[7]','n_d[8]',
                   'n_d[9]','n_d[10]')])
                                             
                                             
