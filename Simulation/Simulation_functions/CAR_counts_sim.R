library(mvtnorm)
library(spdep)      
library(sf)

if (!exists("data_path")) {
  source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))
}

car_sim <- function(tau2){

health_regions <- st_read(data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp"))

nb <- poly2nb(health_regions)
nbw <- nb2listw(nb, style = 'B') # for calculating Moran's I later
nb_info <- nb2WB(nb)

# vector with the number of neighbors for each area
num <- nb_info$num

# W is binary adjacency matrix
W <- nb2mat(nb, style = 'B')

# D is diagonal matrix with elements m_i
D <- diag(num)

n_regions <- nrow(health_regions)

# number of simulations from the CAR model to perform
n_sim <- 1

# Compute precision matrix
Q <- D - W

# eigenvalue decomposition
eigen_Q <- eigen(Q)
V <- eigen_Q$vectors[,1:(n_regions - 1)]
Lambda_pow <- diag(eigen_Q$values[1:(n_regions-1)]^(-1/2))

# independent normal random variables 
u <- matrix(rnorm((n_regions-1) * n_sim,
                  mean = 0, 
                  sd = sqrt(tau2)), nrow = n_regions-1, ncol = n_sim)

y <- V %*% Lambda_pow %*% u
return(list(y=y,nb_info=nb_info))
}
