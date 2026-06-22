library(malariaAtlas)
library(terra)
library(sf)
library(gdistance)
library(tidyverse)
library(exactextractr)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

# Load friction surface (already cropped to Puerto Rico)
friction_pr <- rast(data_path("raw", "geospatial", "friction_surface", "Explorer__2020_motorized_friction_surface_latest_.67.9_17.88_.65.2_18.51_2025_04_11.tiff"))

# Load and prepare hospital points
hospitals <- st_read(data_path("raw", "geospatial", "hospitals", "Hospitals_PR.shp")) %>%
  filter(NAICS_DESC=="GENERAL MEDICAL AND SURGICAL HOSPITALS")

#Replace missing values for beds for simulation

hospitals_range <- hospitals %>% 
  filter(BEDS!=-999) %>% 
  summarize(min=min(BEDS),max=max(BEDS)) 

set.seed(123)
hospitals <- hospitals %>% 
  mutate(BEDS_imputed=ifelse(BEDS==-999,
                             sample((hospitals_range$min:hospitals_range$max)),BEDS))

#Get coordinates 

hospitals <- st_transform(hospitals, crs(friction_pr))
hospital_coords <- st_coordinates(hospitals)

# Convert to gdistance-compatible raster
friction_r <- raster::raster(friction_pr)

# Create transition object from friction (1/time = speed = conductance)
tr <- transition(1 / friction_r, transitionFunction = mean, directions = 8)
tr <- geoCorrection(tr, type = "c")

# Create an empty list to hold each hospital's raster
travel_time_stack <- list()

# Loop over hospitals and calculate travel time surfaces
for (i in 1:nrow(hospital_coords)) {
  hospital_pt <- hospital_coords[i, , drop = FALSE]
  
  # Cost surface for hospital i
  cost_surface <- accCost(tr, hospital_pt)
  
  # Convert to terra raster and add to list
  travel_time_stack[[i]] <- rast(cost_surface)
  names(travel_time_stack)[i] <- hospitals$hospital_name[i]  # or ID field
}

# Combine into a raster stack
travel_time_rasters <- rast(travel_time_stack)

#Import in regions

health_regions_shp<-st_read(data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp"))

# Convert health_regions to SpatVector
health_regions <- vect(health_regions_shp)

# Step 1: Calculate unnormalized weights W_{ij}
weighted_rasters <- list()

for (i in 1:nlyr(travel_time_rasters)) {
  beds <- hospitals$BEDS_imputed[i]
  distance_raster <- travel_time_rasters[[i]]
  
  # Replace zero distances with a tiny nonzero value to avoid division by zero
  distance_raster[distance_raster == 0] <- 1e-6
  
  # Weighted accessibility surface
  weighted_rasters[[i]] <- beds / (distance_raster^2)
}

# Step 2: Stack and sum to get total denominator S_i
weighted_stack <- rast(weighted_rasters)
sum_raster <- app(weighted_stack, sum, na.rm = TRUE)

# Step 3: Normalize each hospital's surface to get P_{ij}
probability_rasters <- list()

for (i in 1:nlyr(weighted_stack)) {
  probability_rasters[[i]] <- weighted_stack[[i]] / sum_raster
}

# Convert list to a raster stack
probability_stack <- rast(probability_rasters)

# Assign CRS from friction surface to avoid projection errors
crs(probability_stack) <- crs(friction_pr)

# Replace NA values (unreachable areas) with 0 to reflect zero probability
probability_stack[is.na(probability_stack)] <- 0

# Step 4: Reproject regions to raster CRS
health_regions_projected <- project(health_regions, crs(probability_stack))

#Import in and get population raster

# Load your population raster
population_raster <- rast(data_path("raw", "geospatial", "population", "dasymetric_population_pr_2020", "Dasymetric_Population_PR_2020_V1.tif"))
population_raster_projected <- project(population_raster, crs(probability_stack), method = "near")
population_resampled <- resample(population_raster_projected, friction_pr, method = "sum")

#Create file of population for each health region

# Make sure CRS matches
health_regions_pop <- st_transform(health_regions_shp, crs(population_raster))

# Calculate total population per region using exactextractr
# This assumes that population_raster pixel values are population counts
pop_summary <- round(exact_extract(population_raster, health_regions_pop, 'sum'))

# Add results to data frame
health_regions_pop$population <- pop_summary

pop_df <- health_regions_pop %>%
  st_drop_geometry() 

write.csv(pop_df,data_path("processed", "simulation_inputs", "region_population.csv"),
          row.names=F)

region_probs_weighted <- list()

for (i in 1:nlyr(probability_stack)) {
  weighted_product <- probability_stack[[i]] * population_resampled
  
  # Sum over each region
  sum_weighted <- terra::extract(weighted_product, health_regions_projected, fun = sum, na.rm = TRUE)
  total_population <- terra::extract(population_resampled, health_regions_projected, fun = sum, na.rm = TRUE)
  
  # Weighted average = sum(probability * population) / sum(population)
  region_avg <- sum_weighted[,2] / total_population[,2]
  
  region_probs_weighted[[ paste0("hosp_", hospitals$ID[i]) ]] <- region_avg
}

# Combine to a dataframe
health_regions_probs_weighted <- cbind(health_regions_shp, as.data.frame(region_probs_weighted))

#Create a file for export

health_regions_probs <- as.data.frame(region_probs_weighted)
health_regions_probs$region<-health_regions_shp$region 
health_regions_probs <- health_regions_probs %>% dplyr::select(region,everything())
write.csv(health_regions_probs,data_path("processed", "simulation_inputs", "health_regions_probs.csv"),
          row.names=F)

#Create a plot for one hospital

library(ggplot2)

# Ensure health_regions_with_probs has the region probability data
# Filter the relevant hospital's data from the probabilities
hospital_id <- paste0("hosp_",23)
hospital_probabilities_weighted <- health_regions_probs_weighted[, hospital_id, drop = FALSE]
hospital_probabilities_weighted$region <- health_regions_shp$region  # Add region names

# Extract the coordinates for hospital ID 53
hospital_location <- hospitals[hospitals$ID == "23", ]

#Crop out of Isla Mona

bbox_no_mona <- st_bbox(c(xmin = -67.3, ymin = 17.8, xmax = -65.2, ymax = 18.6), 
                       crs = st_crs(hospital_probabilities_weighted))
weighted_graph <- st_crop(hospital_probabilities_weighted, bbox_no_mona)


# Plot
 ggplot() +
  # Plot the health regions with hospital probabilities
  geom_sf(data = weighted_graph, aes(fill = hosp_23), color = "white") +
  scale_fill_viridis_c(option = "D", 
                       name = "Presentation Probability", 
                       limits = c(0, 0.09),
                       breaks = seq(0, 0.09, by = 0.03)) +
  theme_classic() +
  # Add labels for each health region (adjust size and column name as needed)
  geom_sf_text(data = weighted_graph, aes(label = region), size = 3, color = "white") +
  # Add a point for hospital ID 53
  geom_sf(data = hospital_location, color = "red", size = 3, shape = 17) +
  theme(legend.position = "bottom") +
  labs(
    x = "Longitude", 
    y = "Latitude", 
  )

#Create data frame for use in the simulation

library(terra)

# 1. Ensure inputs are SpatRaster and SpatVector

health_regions_vect <- vect(health_regions_shp)  # convert to SpatVector if needed

# 2. Combine all rasters into one stack
combined_stack <- c(probability_stack, population_resampled)
names(combined_stack) <- c(paste0("hosp_",hospitals$ID), "population")

# 3. Extract cell coordinates and raster values
raster_df <- as.data.frame(combined_stack, xy = F, cells = TRUE, na.rm = FALSE) %>% 
  filter(!is.na(population) & population!=0)

# 4. Create a raster with health region IDs by rasterizing the vector shapefile
health_region_rast <- rasterize(health_regions_vect, population_resampled, field = "region")  # or appropriate column
health_region_df <- as.data.frame(health_region_rast, xy = F, cells = TRUE)
names(health_region_df)[which(names(health_region_df) == names(health_region_rast))] <- "health_region"

# 5. Merge health region info with raster data
full_df <- merge(raster_df, health_region_df[, c("cell","health_region")], by = "cell")

write.csv(full_df,data_path("processed", "simulation_inputs", "cell_probabilities.csv"),
          row.names=F)

