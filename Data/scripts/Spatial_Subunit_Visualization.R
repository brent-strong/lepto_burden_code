# --- 1. Load Libraries ---
library(terra)
library(sf)
library(tidyverse)
library(ggplot2)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

# --- 2. Data Loading & Preprocessing ---
# File Paths
friction_path <- data_path("raw", "geospatial", "friction_surface", "Explorer__2020_motorized_friction_surface_latest_.67.9_17.88_.65.2_18.51_2025_04_11.tiff")
pop_path <- data_path("raw", "geospatial", "population", "dasymetric_population_pr_2020", "Dasymetric_Population_PR_2020_V1.tif")
regions_path <- data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp")
output_path <- data_path("figures", "spatial_subunits.png")

# Load friction surface to get the target grid and CRS
friction_pr <- rast(friction_path)
target_crs <- crs(friction_pr)

# Load Population
pop_raw <- rast(pop_path)

# Project and resample population to match the friction surface
population_projected <- project(pop_raw, target_crs, method = "near")
population_resampled <- resample(population_projected, friction_pr, method = "sum")
crs(population_resampled) <- target_crs

# Load Health Regions and assign names
health_regions_shp <- st_read(regions_path) %>%
  st_transform(target_crs)
health_regions_shp$region_name <- c("Arecibo", "Bayamón", "Caguas", "Fajardo", "Mayagüez", "Metro", "Ponce")

# --- 3. Crop, MASK, and Format for ggplot2 ---
# Define bounding box (excluding Mona Island)
bbox_sf <- st_bbox(c(xmin = -67.3, ymin = 17.8, xmax = -65.2, ymax = 18.6), crs = target_crs)
health_regions_cropped <- st_crop(health_regions_shp, bbox_sf)

bbox_ext <- ext(-67.3, -65.2, 17.8, 18.6)
pop_cropped <- crop(population_resampled, bbox_ext)

# NEW: Mask the raster to the health regions polygon
# This converts any pixels over the ocean/outside the borders to NA
pop_masked <- mask(pop_cropped, vect(health_regions_cropped))

# Convert the raster to a dataframe for ggplot2
# na.rm = TRUE automatically drops all the NA pixels outside the borders
pop_df <- as.data.frame(pop_masked, xy = TRUE, na.rm = TRUE)

# Dynamically rename the 3rd column to "population"
colnames(pop_df)[3] <- "population"

# (Optional) If you still want to filter out absolute zeros inside PR:
# pop_df <- pop_df %>% filter(population > 0)

# --- 4. Plotting ---
p_pop <- ggplot() +
  geom_tile(data = pop_df, aes(x = x, y = y, fill = population)) +
  
  # Overlay Health Region Boundaries
  geom_sf(data = health_regions_cropped, fill = NA, color = "black", linewidth = 0.5) +
  
  # Add Health Region Labels
  geom_sf_text(data = health_regions_cropped, aes(label = region_name), 
               color = "black", size = 3.5, fontface = "bold", check_overlap = TRUE) +
  
  # Linear color scale
  scale_fill_gradient(
    low = "#E5F5E0",   
    high = "#00441B",  
    name = "Population"
  ) +
  
  coord_sf(expand = FALSE) + 
  theme_void() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold", margin = margin(b = 10)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA) 
  ) +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 15))

print(p_pop)

# --- 5. Save Output ---
ggsave(filename = output_path, plot = p_pop, width = 10, height = 6, dpi = 300, bg = "white")
