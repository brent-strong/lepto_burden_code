
# --- 1. Load Libraries ---
library(terra)
library(sf)
library(gdistance)
library(tidyverse)
library(magick)
library(exactextractr)
library(ggplot2)
library(patchwork)
library(latex2exp) # For LaTeX math labels

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

# --- 2. Data Loading & Preprocessing ---

# Load friction surface and define target CRS
friction_pr <- rast(data_path("raw", "geospatial", "friction_surface", "Explorer__2020_motorized_friction_surface_latest_.67.9_17.88_.65.2_18.51_2025_04_11.tiff"))
target_crs <- crs(friction_pr)

# Load Hospitals
hospitals <- st_read(data_path("raw", "geospatial", "hospitals", "Hospitals_PR.shp")) %>%
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

# Align all spatial layers to the same CRS
hospitals <- st_transform(hospitals, target_crs)
hospital_coords <- st_coordinates(hospitals)
hospital_target <- hospitals %>% filter(ID == "14") 

# Load and project Health Regions
health_regions_shp <- st_read(data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp")) %>%
  st_transform(target_crs)
st_agr(health_regions_shp) <- "constant" # Silence spatial constant warning

# Load Population and resample to match friction surface
population_raster <- rast(data_path("raw", "geospatial", "population", "dasymetric_population_pr_2020", "Dasymetric_Population_PR_2020_V1.tif"))
population_resampled <- project(population_raster, friction_pr, method = "near") %>%
  resample(friction_pr, method = "sum")
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

# --- 5. Generate Scenario Data and Labels ---

scenarios <- list(
  A = list(beds = FALSE, dist = FALSE, label = r"($\gamma=0, \; \delta=0$)"),
  B = list(beds = TRUE,  dist = FALSE, label = r"($\gamma=0, \; \delta=1$)"),
  C = list(beds = FALSE, dist = TRUE,  label = r"($\gamma=2, \; \delta=0$)"),
  D = list(beds = TRUE,  dist = TRUE,  label = r"($\gamma=2, \; \delta=1$)")
)

bbox_no_mona <- st_bbox(c(xmin = -67.3, ymin = 17.8, xmax = -65.2, ymax = 18.6), crs = target_crs)
hosp_idx <- which(hospitals$ID == "14")

results_list <- list()
for (lab in names(scenarios)) {
  scen <- scenarios[[lab]]
  p_stack <- calculate_scenario_probs(travel_time_rasters, hospitals$BEDS_imputed, scen$beds, scen$dist)
  
  target_layer <- p_stack[[hosp_idx]]
  
  # Calculate pop-weighted average
  reg_sum_val <- exact_extract(target_layer * population_resampled, health_regions_shp, 'sum', progress = FALSE)
  reg_pop_total <- exact_extract(population_resampled, health_regions_shp, 'sum', progress = FALSE)
  
  res_df <- health_regions_shp
  res_df$prob <- reg_sum_val / reg_pop_total
  results_list[[lab]] <- st_crop(res_df, bbox_no_mona)
}

# Calculate global max for a shared color scale
global_max <- max(sapply(results_list, function(x) max(x$prob, na.rm = TRUE)))

# --- 1. Prepare Data for Faceting ---

# Combine results into a long-format data frame
map_long <- bind_rows(
  lapply(names(results_list), function(lab) {
    df <- results_list[[lab]] %>%
      st_set_geometry(NULL) %>%         
      mutate(
        Metric = lab,                  # Use the simple name (e.g., "A")
        region = health_regions_shp$NAME 
      )
    df$geometry <- results_list[[lab]]$geometry
    st_as_sf(df)
  })
)
# Find global max probability for consistent fill scale
global_max <- max(map_long$prob, na.rm = TRUE)


# --- 1. Prepare Data for Faceting ---

map_long <- bind_rows(
  lapply(names(results_list), function(lab) {
    df <- results_list[[lab]] %>%
      st_set_geometry(NULL) %>%          # Keep attributes for plotting
      mutate(
        Metric = lab,                     # Store simple label for faceting
        region = health_regions_shp$NAME # Add region names for labeling
      )
    df$geometry <- results_list[[lab]]$geometry
    st_as_sf(df)
  })
)

# Define the 7 regions in the order they appear in your data
region_names <- c("Arecibo", "Bayamón", "Caguas", "Fajardo", "Mayagüez", "Metro", "Ponce")

# Repeat them 4 times and add to map_long
map_long$region <- rep(region_names, times = 4)

# Find global max probability for consistent fill scale
global_max <- max(map_long$prob, na.rm = TRUE)

# 1. Update labels with a comma (using * "," ~ to join them)
facet_labels <- c(
  A = "gamma==0 * ',' ~ delta==0",
  B = "gamma==0 * ',' ~ delta==1",
  C = "gamma==2 * ',' ~ delta==0",
  D = "gamma==2 * ',' ~ delta==1"
)

# 2. Plot with explicit centering
p_faceted <- ggplot(data = map_long) +
  geom_sf(aes(fill = prob), color = "black", size = 0.1) +
  geom_sf_text(aes(label = region), 
               color = "black", 
               size = 2, 
               # check_overlap prevents names from cluttering small facets
               check_overlap = TRUE) +
  facet_wrap(~Metric, 
             ncol = 2, 
             labeller = as_labeller(facet_labels, default = label_parsed)) + 
  scale_fill_gradient(
    low = "#DEEBF7",
    high = "#08519C",
    name = "Average Probability",
    limits = c(0, global_max)
  ) +
  theme_void() +
  theme(
    # hjust = 0.5 ensures the text is centered in the facet box
    strip.text = element_text(size = 12, hjust = 0.5, margin = margin(b = 5)),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 15)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

print(p_faceted)

# --- 3. Save & Add Padding Using magick ---

ggsave(
  filename = data_path("figures", "hospital_presentation_probability.png"),
  plot = p_faceted,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

img <- image_read(data_path("figures", "hospital_presentation_probability.png")) %>%
  image_trim() %>%
  image_border(color = "white", geometry = "50x50")

image_write(img, data_path("figures", "hospital_presentation_probability.png"))

#Create alternative plot

# --- 6. Plot Specific Scenario (Bottom Left: Scenario C) ---

# 1. Filter data for Scenario C only
# In a 2x2 grid (ncol=2), Scenario C is the bottom-left
data_c <- map_long %>% filter(Metric == "C")

# 2. Extract the label for Scenario C to use as a title
title_c <- parse(text = facet_labels["C"])

p_single <- ggplot() +
  # Draw the regions
  geom_sf(data = data_c, aes(fill = prob), color = "black", size = 0.1) +
  # Add region labels
  geom_sf_text(data = data_c, aes(label = region), 
               color = "black", size = 3, check_overlap = TRUE) +
  # Add the Hospital (ID 14) as a red triangle
  # shape = 17 is a filled triangle
  geom_sf(data = hospital_target, color = "red", fill = "red", 
          shape = 17, size = 4) +
  # Use the global max for consistency, or remove limits to auto-scale
  scale_fill_gradient(
    low = "#DEEBF7",
    high = "#08519C",
    name = "Average Probability of Presentation",
    limits = c(0, global_max)
  ) +
  labs(title = title_c) +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold", margin = margin(b = 10)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 15))

print(p_single)
