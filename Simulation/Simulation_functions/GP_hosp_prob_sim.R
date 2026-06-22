library(tidyverse)
library(sf)
library(MASS)   # for mvrnorm
library(scales) # for logistic function

#--- Function to simulate hospital probabilities ---
simulate_hospital_probs <- function(hospitals, sigma, zeta, mean_prob) {
  
  hospitals <- st_transform(hospitals, 32161)  # NAD83 / Puerto Rico and Virgin Is.
  coords <- st_coordinates(hospitals)
  n <- nrow(coords)
  
  # pairwise Euclidean distances (in meters → km)
  dists <- as.matrix(dist(coords)) / 1000
  
  # covariance matrix
  Sigma <- sigma^2 * exp(-zeta * dists)
  
  # mean vector at logit(mean_prob)
  mu <- rep(qlogis(mean_prob), n)
  
  # sample one realization from MVN
  latent <- MASS::mvrnorm(1, mu = mu, Sigma = Sigma)
  
  # transform to probabilities
  probs <- plogis(latent)
  
  return(list(probs=probs,distance=dists))
}
