library(tidyverse)
library(nimble)
library(parallel)
library(coda)
library(readxl)
library(purrr)
library(sf)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

setwd(analysis_path())
example <- read_rds("example_inputs.rds")

hospitals <- st_read(data_path("raw", "geospatial", "hospitals", "Hospitals_PR.shp")) %>%
  filter(NAICS_DESC=="GENERAL MEDICAL AND SURGICAL HOSPITALS")


# Set the path to your folder
folder_path <- data_path("raw", "geospatial", "health_regions")

# Read the specific shapefile
# Note: Point to the .shp file, but the .dbf, .shx, and .prj must be in the same folder.
regions <- st_read(file.path(folder_path, "pr_health_regions.shp"))

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

chain_samples <- readRDS(analysis_path("Output", "results_ratio1_N10000.rds"))
draws <- as.data.frame(as.matrix(chain_samples))
poisson_intercept <- if ("beta_0" %in% names(draws)) "beta_0" else "A"

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
    poisson_intercept,"alpha_0","nu[1]","nu[2]","kappa[1]","kappa[2]",
    "rate_2023[1]","rate_2023[2]","rate_2023[3]",
    "rate_2023[4]","rate_2023[5]","rate_2023[6]",
    "rate_2023[7]","rate_PR_2023","rate_PR_2022",
    "log_sigma_Phi_pois","log_sigma_phi_binom",
    "pi_A","log_sigma_tau"
  )]
)


##Create figures

# Combine chains from mcmc.list
draws <- as.matrix(chain_samples)
draws <- as.data.frame(draws)

# Inverse logit function
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

# Indices of interest
idx <- c(15, 23, 24, 32)

# ---- 1) Summarize inverse-logit(alpha_0) ----
alpha_0_prob <- inv_logit(draws$alpha_0)

alpha_summary <- data.frame(
  parameter = "alpha_0 (inv_logit)",
  mean      = mean(alpha_0_prob),
  sd        = sd(alpha_0_prob),
  q2.5      = quantile(alpha_0_prob, 0.025),
  median    = quantile(alpha_0_prob, 0.5),
  q97.5     = quantile(alpha_0_prob, 0.975)
)

# ---- 2) Summarize inverse-logit(alpha_0 + phi_i * sigma) ----
results <- do.call(rbind, lapply(idx, function(i) {
  
  linear_pred <- draws$alpha_0 +
    draws[[paste0("phi_binom_raw[", i, "]")]] *
    exp(draws$log_sigma_phi_binom)
  
  prob <- inv_logit(linear_pred)
  
  data.frame(
    parameter = paste0("index_", i),
    mean      = mean(prob),
    sd        = sd(prob),
    q2.5      = quantile(prob, 0.025),
    median    = quantile(prob, 0.5),
    q97.5     = quantile(prob, 0.975)
  )
}))

# Combine everything
summary_df <- rbind(alpha_summary, results)



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

p1 <- ggplot(data = map_long) +
  geom_sf(aes(fill = Rate), color = "black", size = 0.1) +
  geom_sf_text(aes(label = region), color = "white", size = 1.5) +
  facet_wrap(~Metric, ncol = 2) +
  scale_fill_gradient(
    low = "#FDE0DD",    # light soft pink
    high = "#A50F15",   # deep muted red
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

image_write(img, "PR_Health_Regions_Map_sensitivity_Cropped.png")
