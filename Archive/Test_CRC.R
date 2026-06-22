library(tidyverse)
library(nimble)
library(nimbleHMC)

sim_results <- readRDS("sim_results.rds")

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
    
    n_rt <- length(passive_idx)
    nRemain <- length(active_idx)
    
    if(n_rt > 0) {
      Y_passive <- test_mat[passive_idx, -1, drop = FALSE]
    }
    if(nRemain > 0) {
      Y_active  <- test_mat[active_idx, -1, drop = FALSE]
    }
    
    # Step 3: compute log-term for passive observations
    log_term_passive <- 0
    if (n_rt > 0) {
      for (j in 1:n_rt) {
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


dMarginalized <- nimbleFunction(
  run = function(x = double(2),
                 n_rt = double(0),   # n^{RT}
                 n_t  = double(0),   # n^{T}
                 n_r  = double(0),   # n^{R}
                 nu   = double(1),
                 kappa = double(1),
                 piDR = double(0),
                 piDT = double(0),
                 lambda_D = double(0),
                 max_nD = double(0, default = 200),
                 lgamma_table=double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    if(n_rt < 0 | n_t < n_rt) {
      if(log) return(-Inf) else return(0.0)
    }
    
    ## likelihood of Y given nDT
    log_lik <- log_lik_Y(sens = nu, spec = kappa, n_t=n_t, test_mat = x)
    total_log_prob <- -Inf
    
    ## precompute logs
    log_piDR <- log(piDR)
    log1m_piDR <- log(1 - piDR)
    log_piDT <- log(piDT)
    log1m_piDT <- log(1 - piDT)
    
    #Pre-compute gamma values
    
    a1 <- n_rt+1 
    a2 <- n_r - n_rt+1
    lgamma_a1<-lgamma_table[a1]
    lgamma_a2<-lgamma_table[a2]
    
    upper_nDT <- max(min(n_t, qpois(1 - 10^(-15), lambda_D)),n_rt)
    
    
    ## sum over nDT
    for(nDT in n_rt:upper_nDT){
      idx <- nDT - n_rt + 1
      log_prob_Y_given_nDT <- log_lik[idx]
      
      nD_min <- nDT + n_r - n_rt
      inner_log_sum <- -Inf
      a3 <- nDT - n_rt+1
      lgamma_a3<-lgamma_table[a3]
      upper_nD <- max(min(max_nD, qpois(1 - 10^(-15), lambda_D)),nD_min)
      
      for(nD in nD_min:upper_nD) {
        
        a4 <- nD - (n_r + nDT - n_rt) + 1
        lgamma_a4 <- lgamma_table[a4]
        log_mult_coeff <- lgamma(nD + 1.0) -
          (lgamma_a1 +
             lgamma_a2 +
             lgamma_a3 +
             lgamma_a4)
        
        log_prob_cells <- a1 * (log_piDR + log_piDT) +
          a2 * (log_piDR + log1m_piDT) +
          a3 * (log1m_piDR + log_piDT) +
          a4 * (log1m_piDR + log1m_piDT)
        
        log_p_nD <- dpois(nD, lambda_D, log = TRUE)
        
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
    BUGSdist = "dMarginalized(n_rt, n_t, n_r, nu, kappa, piDR, piDT, lambda_D, max_nD,
    lgamma_table)",
    types = c(
      "value = double(2)",
      "n_rt = double(0)",
      "n_t  = double(0)",
      "n_r  = double(0)",
      "nu   = double(1)",
      "kappa = double(1)",
      "piDR = double(0)",
      "piDT = double(0)",
      "lambda_D = double(0)",
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
  #Sampled hospitals
  
  for(t in 1:T){
    for(h in 1:H_t){
      for(s in 1:S){
        lambda_d[s,h,t]<-pop[s]*exp(log_lambda[s,t])*prob_pres_s[s,h]
        Y[,,(t-1)*(H_t*S)+(h-1)*S+s] ~ dMarginalized(n_rt=n_rt[s,h,t],
                                                     n_t=n_t[s,h,t],
                                                     n_r=n_r[s,h,t],
                                                     nu=nu[1:2],
                                                     kappa=kappa[1:2],
                                                     piDR=pi_DR_s[h],
                                                     piDT=pi_DT,
                                                     lambda_D=lambda_d[s,h,t],
                                                     max_nD=max_n_d,
                                                     lgamma_table=lgamma_table[])
      }
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
  nburnin = 0,
  niter = 100,
  nchains = 1,
  samplesAsCodaMCMC = TRUE
)

plot(chain_samples[,c("A")],density=F)



#Test of dmarginalized


dMarginalized(sim_results$data_list$Y[,,1],n_rt=0,n_t=5,n_r=0,nu=c(0.6,0.9),
                kappa=c(0.9,0.9),piDR=0.33,piDT=0.75,lambda_D=0.006157309*exp(2)*4.83998,
                lgamma_table=lgamma(0:1000+1),log=T) 

    


