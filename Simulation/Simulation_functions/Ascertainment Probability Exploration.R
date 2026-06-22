library(tidyverse)
library(sf)
library(MASS)   # for mvrnorm
library(scales) # for logistic function

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

#--- Load shapefiles ---
health_regions <- st_read(data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp"))

hospitals <- st_read(data_path("raw", "geospatial", "hospitals", "Hospitals_PR.shp")) %>%
  filter(NAICS_DESC == "GENERAL MEDICAL AND SURGICAL HOSPITALS")

# Ensure CRS matches and get projected coordinates (for Euclidean distances in km)
hospitals <- st_transform(hospitals, 32161)  # NAD83 / Puerto Rico and Virgin Is.
hospital_coords <- st_coordinates(hospitals)

#--- Function to simulate hospital probabilities ---
simulate_hospital_probs <- function(coords, sigma = 1, zeta = 0.1, mean_prob = 0.3) {
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
  
  return(probs)
}

# Example: simulate with sigma=1, zeta=0.05

hospitals$prob <- simulate_hospital_probs(hospital_coords, sigma = 0.5, zeta = 0.03)

#--- Plot ---
ggplot() +
  geom_sf(data = health_regions, fill = "grey95", color = "white") +
  geom_sf(data = hospitals, aes(color = prob), size = 3) +
  scale_color_viridis_c(option = "plasma", name = "Sampled\nprobability") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "transparent"))
