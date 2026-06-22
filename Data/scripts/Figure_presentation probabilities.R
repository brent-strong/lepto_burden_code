# --- 1. Load Libraries ---
library(terra)
library(sf)
library(gdistance)
library(tidyverse)
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

# --- 6. Plotting ---

plot_list <- list()
for (lab in names(scenarios)) {
  scen <- scenarios[[lab]]
  
  p <- ggplot() +
    geom_sf(data = results_list[[lab]], aes(fill = prob), color = "white", linewidth = 0.2) +
    geom_sf(data = hospital_target, color = "red", size = 2.5, shape = 17) +
    scale_fill_viridis_c(option = "D", 
                         name = NULL, 
                         limits = c(0, global_max)) +
    labs(title = TeX(scen$label)) +
    theme_void() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 15, hjust = 0.05, vjust = -1),
      plot.margin = margin(0, 0, 0, 0, "pt")
    ) + 
    coord_sf(
      xlim = c(bbox_no_mona["xmin"], bbox_no_mona["xmax"]),
      ylim = c(bbox_no_mona["ymin"], bbox_no_mona["ymax"]),
      expand = FALSE,
      default_crs = NULL
    )
  
  
  plot_list[[lab]] <- p
}

# --- 7. Final Assembly ---

final_plot <- (plot_list$A + plot_list$B) / (plot_list$C + plot_list$D) +
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    panel.spacing = unit(0, "pt"),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )

print(final_plot)

ggsave(
  filename = data_path("figures", "hospital_presentation_probability.png"),
  plot = final_plot,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white")
  
