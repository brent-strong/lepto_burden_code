# --- 1. Load Libraries ---
library(terra)
library(sf)
library(gdistance)
library(tidyverse)
library(magick)
library(exactextractr)
library(ggplot2)
library(patchwork)
library(latex2exp)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

# --- 2. Data Loading & Preprocessing ---

# File Paths
friction_path <- data_path("raw", "geospatial", "friction_surface", "Explorer__2020_motorized_friction_surface_latest_.67.9_17.88_.65.2_18.51_2025_04_11.tiff")
hospitals_path <- data_path("raw", "geospatial", "hospitals", "Hospitals_PR.shp")
regions_path <- data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp")
pop_path <- data_path("raw", "geospatial", "population", "dasymetric_population_pr_2020", "Dasymetric_Population_PR_2020_V1.tif")
output_path <- data_path("figures", "hospital_presentation_probability.png")

# Load friction surface and define target CRS
friction_pr <- rast(friction_path)
target_crs <- crs(friction_pr)

# Load and Filter Hospitals
hospitals <- st_read(hospitals_path) %>%
  filter(NAICS_DESC == "GENERAL MEDICAL AND SURGICAL HOSPITALS")

# Impute missing bed values
hospitals_range <- hospitals %>% 
  filter(BEDS != -999) %>% 
  summarize(min = min(BEDS), max = max(BEDS)) 

set.seed(123)
hospitals <- hospitals %>% 
  mutate(BEDS_imputed = ifelse(BEDS == -999,
                               sample(hospitals_range$min:hospitals_range$max, n(), replace = TRUE), 
                               BEDS))

# Project Hospitals and Define Plotting Groups
hospitals <- st_transform(hospitals, target_crs)
hospitals_all <- hospitals %>%
  mutate(plot_group = case_when(
    ID == "14"   ~ "14",
    ID == "15"   ~ "15",
    ID == "2512" ~ "2512",
    ID == "23"   ~ "23",
    TRUE         ~ "Other"
  ))

hospital_coords <- st_coordinates(hospitals_all)
hospital_target <- hospitals_all %>% filter(ID == "14") 

# Load and project Health Regions
health_regions_shp <- st_read(regions_path) %>%
  st_transform(target_crs)
st_agr(health_regions_shp) <- "constant"

# Load Population and resample
pop_raw <- rast(pop_path)
population_projected <- project(pop_raw, target_crs, method = "near")
population_resampled <- resample(population_projected, friction_pr, method = "sum")
crs(population_resampled) <- target_crs

# --- 3. Travel Time Surface (gdistance) ---

friction_r <- raster::raster(friction_pr)
tr <- transition(1 / friction_r, transitionFunction = mean, directions = 8)
tr <- geoCorrection(tr, type = "c")

travel_time_list <- list()
for (i in 1:nrow(hospital_coords)) {
  hospital_pt <- hospital_coords[i, , drop = FALSE]
  cost_surface <- accCost(tr, hospital_pt)
  r <- rast(cost_surface)
  crs(r) <- target_crs
  travel_time_list[[i]] <- r
}
travel_time_rasters <- rast(travel_time_list)

# --- 4. Probability Function ---

calculate_scenario_probs <- function(tt_rasters, beds, use_beds, use_dist) {
  weighted_list <- list()
  for (i in 1:nlyr(tt_rasters)) {
    d <- tt_rasters[[i]]
    d[d == 0] <- 1e-6
    w_cap <- if(use_beds) beds[i] else 1
    if(use_dist) {
      weighted_list[[i]] <- w_cap / (d^2)
    } else {
      weighted_list[[i]] <- (d * 0 + 1) * w_cap
    }
  }
  w_stack <- rast(weighted_list)
  sum_r <- app(w_stack, sum, na.rm = TRUE)
  probs <- w_stack / sum_r
  probs[is.na(probs)] <- 0
  crs(probs) <- target_crs
  return(probs)
}

# --- 5. Generate Scenario Data ---

scenarios <- list(
  A = list(beds = FALSE, dist = FALSE),
  B = list(beds = TRUE,  dist = FALSE),
  C = list(beds = FALSE, dist = TRUE),
  D = list(beds = TRUE,  dist = TRUE)
)

bbox_no_mona <- st_bbox(c(xmin = -67.3, ymin = 17.8, xmax = -65.2, ymax = 18.6), crs = target_crs)
hosp_idx <- which(hospitals_all$ID == "14")

results_list <- list()
for (lab in names(scenarios)) {
  scen <- scenarios[[lab]]
  p_stack <- calculate_scenario_probs(travel_time_rasters, hospitals_all$BEDS_imputed, scen$beds, scen$dist)
  
  target_layer <- p_stack[[hosp_idx]]
  
  reg_sum_val <- exact_extract(target_layer * population_resampled, health_regions_shp, 'sum', progress = FALSE)
  reg_pop_total <- exact_extract(population_resampled, health_regions_shp, 'sum', progress = FALSE)
  
  res_df <- health_regions_shp
  res_df$prob <- reg_sum_val / reg_pop_total
  results_list[[lab]] <- st_crop(res_df, bbox_no_mona)
}

# Combine for Faceting
region_names <- c("Arecibo", "Bayamón", "Caguas", "Fajardo", "Mayagüez", "Metro", "Ponce")
map_long <- bind_rows(
  lapply(names(results_list), function(lab) {
    df <- results_list[[lab]]
    df$Metric <- lab
    df$region_name <- region_names # Use names for labeling
    return(df)
  })
)

global_max <- max(map_long$prob, na.rm = TRUE)

# --- 6. PLOT 1: Faceted Comparison ---

facet_labels <- c(
  A = "gamma==0 * ',' ~ delta==0",
  B = "gamma==0 * ',' ~ delta==1",
  C = "gamma==2 * ',' ~ delta==0",
  D = "gamma==2 * ',' ~ delta==1"
)

p_faceted <- ggplot(data = map_long) +
  geom_sf(aes(fill = prob), color = "black", size = 0.1) +
  geom_sf_text(aes(label = region_name), color = "black", size = 2, check_overlap = TRUE) +
  # Hospital 14 as a red triangle on every facet
  geom_sf(data = hospital_target, color = "red", fill = "red", shape = 17, size = 3) + 
  facet_wrap(~Metric, ncol = 2, labeller = as_labeller(facet_labels, default = label_parsed)) + 
  scale_fill_gradient(low = "#DEEBF7", high = "#08519C", name = "Weight", limits = c(0, global_max)) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, hjust = 0.5, margin = margin(b = 5)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 15))

print(p_faceted)

# --- 7. PLOT 2: Single Scenario (C) with All Hospitals ---

# 1. Prepare hospital data with a specific plotting order
# We set "Other" first so it's drawn at the bottom, and "14" last so it's on top.
hospitals_ordered <- hospitals_all %>%
  mutate(plot_group = factor(plot_group, 
                             levels = c("Other", "15", "2512", "23", "14"))) %>%
  arrange(plot_group)

# 2. Filter data for Scenario C
data_c <- map_long %>% filter(Metric == "C")
title_c <- parse(text = facet_labels["C"])

# 3. Plot p_single
p_single <- ggplot() +
  # Background regions
  geom_sf(data = data_c, aes(fill = prob), color = "black", size = 0.1) +
  geom_sf_text(data = data_c, aes(label = region_name), 
               color = "black", size = 3, check_overlap = TRUE) +
  
  # Hospitals plotted in the order defined above
  geom_sf(data = hospitals_ordered, aes(color = plot_group), shape = 17, size = 3) +
  
  # Color scale with legend suppressed
  scale_color_manual(
    values = c(
      "Other" = "black",  # Generic hospitals
      "15"    = "purple", # Specific yellow group
      "2512"  = "purple", 
      "23"    = "purple", 
      "14"    = "red"     # Target hospital on top
    ),
    guide = "none"        # Suppress the hospital legend
  ) +
  
  # Probability fill scale
  scale_fill_gradient(
    low = "#DEEBF7", 
    high = "#08519C", 
    name = "Weight", 
    limits = c(0, global_max)
  ) +
  
  labs(title = title_c) +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", barwidth = 15)
  )

print(p_single)
# --- 8. Save & Post-process ---

ggsave(filename = output_path, plot = p_faceted, width = 10, height = 8, dpi = 300, bg = "white")

img <- image_read(output_path) %>%
  image_trim() %>%
  image_border(color = "white", geometry = "50x50")

image_write(img, output_path)
