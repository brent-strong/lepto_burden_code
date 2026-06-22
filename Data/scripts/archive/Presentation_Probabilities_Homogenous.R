library(malariaAtlas)
library(terra)
library(sf)
library(gdistance)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

PR_regions <- st_read(data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp"))
raster <- getRaster(dataset_id=c("Explorer__2020_motorized_friction_surface"),
                    file_path=data_path("raw", "geospatial", "friction_surface"),
                    shp=PR_regions)

#Calculate travel times 

library(terra)
library(sf)
library(gdistance)

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

# Step 5: Average probability for each region and hospital
region_probs <- list()

for (i in 1:nlyr(probability_stack)) {
  region_means <- extract(probability_stack[[i]], health_regions_projected, fun = mean, na.rm = TRUE)
  
  # Check if region_means has more than 1 column before trying to access it
  region_probs[[ paste0("hosp_",hospitals$ID[i]) ]] <- region_means[,2]
}

# Combine all hospital region averages into a dataframe
health_regions_probs <- cbind(health_regions_shp,as.data.frame(region_probs))

#Create a plot for one hospital

library(ggplot2)

# Ensure health_regions_with_probs has the region probability data
# Filter the relevant hospital's data from the probabilities
hospital_id <- paste0("hosp_",23)
hospital_probabilities <- health_regions_probs[, hospital_id, drop = FALSE]
hospital_probabilities$region <- health_regions_shp$region  # Add region names

# Ensure hospital coordinates are in the same CRS as the regions
hospital_coords_sf <- st_as_sf(hospitals, coords = c("longitude", "latitude"), crs = crs(hospital_probabilities))

# Extract the coordinates for hospital ID 53
hospital_location <- hospital_coords_sf[hospital_coords_sf$ID == "23", ]

#Crop out of Isla Mona

bbox_no_mona <- st_bbox(c(xmin = -67.3, ymin = 17.8, xmax = -65.2, ymax = 18.6), 
                        crs = st_crs(hospital_probabilities))
graph <- st_crop(hospital_probabilities, bbox_no_mona)


# Plot
plot_unweighted <- ggplot() +
  # Plot the health regions with hospital probabilities
  geom_sf(data = graph, aes(fill = hosp_23), color = "white") +
  scale_fill_viridis_c(option = "D", 
                       name = "Presentation Probability", 
                       limits = c(0, 0.09),
                       breaks = seq(0, 0.09, by = 0.03)) +
  theme_classic() +
  # Add labels for each health region (adjust size and column name as needed)
  geom_sf_text(data = graph, aes(label = region), size = 3, color = "white") +
  # Add a point for hospital ID 53
  geom_sf(data = hospital_location, color = "red", size = 3, shape = 17) +
  theme(legend.position = "bottom") +
  labs(
    x = "Longitude", 
    y = "Latitude", 
  )


ggsave(filename = data_path("figures", "presentation_probabilities_unweighted.png"),
       plot = plot_unweighted, width = 10, height = 5.625, dpi = 300)




