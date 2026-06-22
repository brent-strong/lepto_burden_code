library(tidyverse)
library(nimble)
library(nimbleHMC)

sim_results <- readRDS("sim_results.rds")

m_code <- nimbleCode({
  
  for(h in 1:H){
    pi_DR[h] <- expit(alpha)
  }
  
  for (h in 1:H_t) {
    pi_DR_s[h] <- pi_DR[surveilled_indices[h]]
  }
  
  for (h in 1:Not_H_t) {
    pi_DR_ns[h] <- pi_DR[not_surveilled_indices[h]]
  }
  
  for(t in 1:T){
    for(s in 1:S){
      log_lambda[s,t] <- A
      y[s,t] ~ dpois(pop[s]*exp(log_lambda[s,t])*sum(prob_pres_ns[s,]*pi_DR_ns[]))
    }
  }
  
  #Priors
  A ~ dnorm(0,10)
  alpha <- logit(0.33)
  rho <- 0.06
  nu[1] <- 0.6
  nu[2] <- 0.9
  kappa[1] <- 0.9
  kappa[2] <- 0.9
  pi_DT <-0.75
})

#Fit and run model

n_r <- sim_results$constant_list$n_r
n_rt <- sim_results$constant_list$n_rt
n_t <- sim_results$constant_list$n_t

# Create function to initialize Phi_pois and sigma_phi_pois



# Create final inits function

inits_function <- function() {
  
  return(list(
    A=2
  ))
}

# Initialize and configure model

Rmodel <- nimbleModel(
  m_code,
  constants = sim_results$constant_list,
  data = sim_results$data_list,
  inits=inits_function(),
  dimensions = list(
    prob_pres_s = c(7, sim_results$constant_list$H_t),       # or (S, H_t) depending on how you define it
    prob_pres_ns = c(7,sim_results$constant_list$Not_H_t),
    pi_DR_ns = sim_results$constant_list$Not_H_t)
)



conf <- configureMCMC(Rmodel)


Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

chain_samples <- runMCMC(
  Cmcmc,
  nburnin = 5000,
  niter = 10000,
  nchains = 1,
  inits=inits_function(),
  samplesAsCodaMCMC = TRUE
)

plot(chain_samples[,c("A")],density=F)
summary(chain_samples[,c("A")],density=F)
