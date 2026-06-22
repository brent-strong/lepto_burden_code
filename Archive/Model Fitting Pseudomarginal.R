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

split_df <- split(surveillance[, 5:6], surveillance$hospital)
mat_conv <- function(matrix){
  matrix <- as.matrix(matrix)
  return(matrix)
}


# Step 2: Convert each subset to a matrix
matrix_test <- lapply(split_df, mat_conv)


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


#Create distribution for array D

#Define the density function as a nimbleFunction

dFixedSum <- nimbleFunction(
  run = function(x = double(2), n_p = double(0),n_c=double(0), 
                 n_a=double(0),S=double(0),log = integer(0, default = 0)) {
    returnType(double(0))
    flag <- rep(0.0,S)
    for(s in 1:S){
      sum_x <- sum(x[,s])
      if (sum_x != n_p | sum_x < n_c | sum_x > n_a) {
        flag[s]<-1
      } 
    }
    if (sum(flag) != 0) {
      if (log) return(-Inf)
      else return(0.0)
    }
    else{
      log_prob <- S*(lfactorial(n_p-n_c) + lfactorial(n_a-n_p) - lfactorial(n_a-n_c)) 
      if(log) return(log_prob)
      else return(exp(log_prob))
    }
  })


#Register the custom distribution 

registerDistributions(list(
  dFixedSum = list(
    BUGSdist = "dFixedSum(n_p,n_c,n_a,S)",
    types = c("value = double(2)",    # this tells NIMBLE that x is a vector
              "n_p = double(0)",
              "n_c = double(0)",
              "n_a=double(0)",
              "S=double(0)"),
    discrete = TRUE
  )
))

#Create distribution for vector of test results


dTests <- nimbleFunction(
  run = function(x = double(2), nu = double(1), kappa = double(1),
                 status = double(2), S = double(0), K = double(0),
                 n_a=double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    log_prob_sample <- rep(0.0, S)
    for (s in 1:S) {
      log_prob_individual <- rep(0.0, n_a)
      for (i in 1:n_a) {
        log_prob_individual_test <- rep(0.0,K)
        for(k in 1:K){
          if (x[i,k] == -1) {
            log_prob_individual_test[k] <- 0.0
          } else {
            if (status[i, s] == 1) {
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
      log_prob_sample[s] <- sum(log_prob_individual[1:n_a])
    }
    log_prob <- log(1/S) + max(log_prob_sample[1:S]) + log(sum(exp(log_prob_sample[1:S]-max(log_prob_sample[1:S]))))
    if (log) return(log_prob)
    else return(exp(log_prob))
  },
)



registerDistributions(list(
  dTests = list(
    BUGSdist = "dTests(nu,kappa,status,S,K,n_a)",
    types = c("value = double(2)", "nu = double(1)", "kappa = double(1)", "status = double(2)", "S = double(0)",
              "K = double(0)", "n_a = double(0)"),
    discrete = TRUE
  )
))

#Create sampler for array D

joint_RW <- nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    S <- control$S
    n_a <- control$n_a
    n_c <- control$n_c
  },
  
  run = function() {
    # initial model logProb
    log_MH_ratio <- 0.0
    model_lp_proposed <- 0.0
    model_lp_initial <- 0.0
    n_p <- 0.0
    proposal_n_p <- 0.0
    model_lp_initial <- model$getLogProb(calcNodes)
    u <- runif(1, 0, 1)
    n_p <- values(model, target[1])[1]
    D <- values(model, target[2])
    
    if(u>=0.5 & n_p<n_a){
      proposal_n_p <- n_p+1
      proposal_D_vec <- values(model, target[2])
      proposal_D <- matrix(proposal_D_vec, nrow = n_a, ncol = S)
      for(s in 1:S){
        ind <- which(proposal_D[,s]==0)
        j <- ceiling(runif(1,0,length(ind)))
        proposal_D[ind[j],s]<-1
        for(i in (n_c+1):(n_a-1)){
          j <- ceiling(runif(1,i-1,n_a))
          temp_i <- proposal_D[i,s]
          temp_j <- proposal_D[j,s]
          proposal_D[i,s] <- temp_j
          proposal_D[j,s] <- temp_i  
        }
      }
      
      # generate proposal
      # store proposal into model
      values(model, target[1]) <<- c(proposal_n_p)
      values(model, target[2]) <<- c(proposal_D)
      
      # proposal model logProb
      model_lp_proposed <- model$calculate(calcNodes)
      
      # log-Metropolis-Hastings ratio
      log_MH_ratio <- model_lp_proposed + S*(lfactorial(n_a-n_c) - 
                                               lfactorial(proposal_n_p-n_c) - lfactorial(n_a-proposal_n_p)) -
        model_lp_initial - S*(lfactorial(n_a-n_c) - 
                                lfactorial(n_p-n_c) - lfactorial(n_a-n_p))
    }
    
    
    if(u>=0.5 & n_p==n_a){
      proposal_n_p <- n_p
      proposal_D_vec <- values(model, target[2])
      proposal_D <- matrix(proposal_D_vec, nrow = n_a, ncol = S)
      
      # generate proposal
      # store proposal into model
      values(model, target[1]) <<- c(proposal_n_p)
      values(model, target[2]) <<- c(proposal_D)
      
      # proposal model logProb
      model_lp_proposed <- model$calculate(calcNodes)
      
      # log-Metropolis-Hastings ratio
      log_MH_ratio <- model_lp_proposed  - model_lp_initial
    }
    
    
    if(u<0.5 & n_p==n_c){
      proposal_n_p <- n_p
      proposal_D_vec <- values(model, target[2])
      proposal_D <- matrix(proposal_D_vec, nrow = n_a, ncol = S)
      
      # generate proposal
      # store proposal into model
      values(model, target[1]) <<- c(proposal_n_p)
      values(model, target[2]) <<- c(proposal_D)
      
      # proposal model logProb
      model_lp_proposed <- model$calculate(calcNodes)
      
      # log-Metropolis-Hastings ratio
      log_MH_ratio <- model_lp_proposed - model_lp_initial
    }
    
    
    if(u<0.5 & n_p>n_c){
      proposal_n_p <- n_p-1
      proposal_D_vec <- values(model, target[2])
      proposal_D <- matrix(proposal_D_vec, nrow = n_a, ncol = S)
      for(s in 1:S){
        ind <- which(proposal_D[,s]==1)
        j <- ceiling(runif(1,0,length(ind)))
        proposal_D[ind[j],s]<-0
        for(i in (n_c+1):(n_a-1)){
          j <- ceiling(runif(1,i-1,n_a))
          temp_i <- proposal_D[i,s]
          temp_j <- proposal_D[j,s]
          proposal_D[i,s] <- temp_j
          proposal_D[j,s] <- temp_i  
        }
      }
      # generate proposal
      # store proposal into model
      values(model, target[1]) <<- c(proposal_n_p)
      values(model, target[2]) <<- c(proposal_D)
      
      # proposal model logProb
      model_lp_proposed <- model$calculate(calcNodes)
      
      # log-Metropolis-Hastings ratio
      log_MH_ratio <- model_lp_proposed + S*(lfactorial(n_a-n_c) - 
                                               lfactorial(proposal_n_p-n_c) - lfactorial(n_a-proposal_n_p)) -
        model_lp_initial - S*(lfactorial(n_a-n_c) - 
                                lfactorial(n_p-n_c) - lfactorial(n_a-n_p))
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
    n_c[h] ~ dbinom(prob=pi_s[h],size=n_p[h])
    pi_s[h] <- expit(alpha + phi_binom_s[h])
    phi_binom_s[h] ~ dnorm(0,sd=sigma_phi_binom)
    n_p[h] ~ dpois(sum(pop[1:R]*exp(lambda[1:R])*prob_pres_s[,h]))
  }
  
  #Testing data
  
  Y1[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D1[1:n_a[1],1:S],S=S,K=K,n_a=n_a[1])
  Y2[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D2[1:n_a[2],1:S],S=S,K=K,n_a=n_a[2])
  Y3[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D3[1:n_a[3],1:S],S=S,K=K,n_a=n_a[3])
  Y4[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D4[1:n_a[4],1:S],S=S,K=K,n_a=n_a[4])
  Y5[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D5[1:n_a[5],1:S],S=S,K=K,n_a=n_a[5])
  Y6[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D6[1:n_a[6],1:S],S=S,K=K,n_a=n_a[6])
  Y7[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D7[1:n_a[7],1:S],S=S,K=K,n_a=n_a[7])
  Y8[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D8[1:n_a[8],1:S],S=S,K=K,n_a=n_a[8])
  Y9[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D9[1:n_a[9],1:S],S=S,K=K,n_a=n_a[9])
  Y10[,] ~ dTests(nu=nu[1:K],kappa=kappa[1:K],status=D10[1:n_a[10],1:S],S=S,K=K,n_a=n_a[10])
  
  
  
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
  
  
  D1[1:n_a[1],1:S] ~ dFixedSum(n_p=n_p[1], n_c=n_c[1],n_a=n_a[1],S=S)
  D2[1:n_a[2],1:S] ~ dFixedSum(n_p=n_p[2], n_c=n_c[2],n_a=n_a[2],S=S)
  D3[1:n_a[3],1:S] ~ dFixedSum(n_p=n_p[3], n_c=n_c[3],n_a=n_a[3],S=S)
  D4[1:n_a[4],1:S] ~ dFixedSum(n_p=n_p[4], n_c=n_c[4],n_a=n_a[4],S=S)
  D5[1:n_a[5],1:S] ~ dFixedSum(n_p=n_p[5], n_c=n_c[5],n_a=n_a[5],S=S)
  D6[1:n_a[6],1:S] ~ dFixedSum(n_p=n_p[6], n_c=n_c[6],n_a=n_a[6],S=S)
  D7[1:n_a[7],1:S] ~ dFixedSum(n_p=n_p[7], n_c=n_c[7],n_a=n_a[7],S=S)
  D8[1:n_a[8],1:S] ~ dFixedSum(n_p=n_p[8], n_c=n_c[8],n_a=n_a[8],S=S)
  D9[1:n_a[9],1:S] ~ dFixedSum(n_p=n_p[9], n_c=n_c[9],n_a=n_a[9],S=S)
  D10[1:n_a[10],1:S] ~ dFixedSum(n_p=n_p[10], n_c=n_c[10],n_a=n_a[10],S=S)
  
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
                  pop=region_level_data$population/100000)

constant_list <- list(R=nrow(region_cases),
                      prob_pres_ns=as.matrix(pres_probs_region_non_surveilled_hospitals[,-1]),
                      prob_pres_s=as.matrix(pres_probs_region_surveilled_hospitals),
                      H_t=length(n_c),
                      Not_H_t=ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[,-1])),
                      S=5000,
                      K=2,
                      n_a=n_a)

#Create initialization functions for parameters

D_np_init_function <- function(S, n_a, n_c,i) {
  n_p <- rep(0, length(n_a))
  
  na <- n_a[i]
  nc <- n_c[i]
  
  # np must be at least nc and at most na
  np <- nc
  
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

D1<-D_np_init_function(S=200, n_a=n_a[1], n_c=n_c[1], i=1)

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
  return(list(D1 = D_np1$D,
              D2 = D_np2$D,
              D3 = D_np3$D,
              D4 = D_np4$D,
              D5 = D_np5$D,
              D6 = D_np6$D,
              D7 = D_np7$D,
              D8 = D_np8$D,
              D9 = D_np9$D,
              D10 = D_np10$D
              ,n_p=c(D_np1$n_p,
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
              A=rnorm(1,0,10),alpha=rnorm(1,0,10),nu=c(0.6,0.9),
              kappa=c(0.9,0.9)))
}


#Initialize and configure model

Rmodel <- nimbleModel(
  m_code,
  constants = constant_list,
  data = data_list,
  dimensions = list(pi_ns=ncol(as.matrix(pres_probs_region_non_surveilled_hospitals[,-1]))),
  buildDerivs = TRUE)

## Assume you already have n_a as a vector, length 10 (one per Di)
D_nodes <- sprintf("D%d[1:%d, 1:5000]", seq_along(n_a), n_a)

allStoch <- Rmodel$getNodeNames(stochOnly = TRUE)

manualNodes <- c(sprintf("n_p[%d]", seq_along(n_a)), D_nodes,"nu[1]","nu[2]","kappa[1]","kappa[2]")

autoNodes <- setdiff(allStoch, manualNodes)

conf <- configureHMC(Rmodel, nodes = autoNodes)

for (i in 1:10) {
  conf$addSampler(target = c(sprintf('n_p[%d]', i),
                              sprintf('D%d', i)),
                   type = 'joint_RW',
                   control = list(S = 5000,
                                  n_a = n_a[i],
                                  n_c = n_c[i]))
}

conf$addSampler(target = "nu[1]",
                     type = 'RW')
conf$addSampler(target = "nu[2]",
                     type = 'RW')
conf$addSampler(target = "kappa[1]",
                     type = 'RW')
conf$addSampler(target = "kappa[2]",
                     type = 'RW')



## Then add a single joint_RW sampler

conf$addMonitors(c("alpha", "A", "sigma_Phi_pois","sigma_phi_binom","n_p","nu","kappa"))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


start <- Sys.time()

samples <- runMCMC(Cmcmc,
                   nburnin = 2500,
                   niter=5000,
                   nchains = 2,
                   inits = inits_function(S=5000,n_a=n_a,n_c=n_c,upper_pois=4,upper_binom=4),
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
