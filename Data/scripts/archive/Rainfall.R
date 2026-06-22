# Required libraries
library(terra)
library(sf)
library(dplyr)
library(lubridate)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

# Path to your region shapefile
health_region <- st_read(data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp")) %>%
  st_transform(crs = 4326)  # adjust as needed

# Target date for year-to-date: Dec 31, 2023
target_date <- as.Date("2023-10-1")
yr  <- format(target_date, "%Y")
mo  <- format(target_date, "%m")
dy  <- format(target_date, "%d")
ymd_str <- format(target_date, "%Y%m%d")

# Specify period and format
period <- "wytd"
fmt    <- "tif"

# Construct the URL
url <- sprintf(
  "https://water.noaa.gov/resources/downloads/precip/stageIV/%s/%s/%s/nws_precip_%s_%s_pr.%s",
  yr, mo, dy, period, ymd_str, fmt
)

# Download (to temp), read as raster
dest <- tempfile(fileext = paste0(".", fmt))
tryCatch({
  download.file(url, dest, mode = "wb", quiet = TRUE)
  qpe_raster <- rast(dest)
}, error = function(e) {
  stop("Failed to download or read raster:\n", url)
})

# Crop to your health regions
# Transform shapefile to match raster
health_region_proj <- st_transform(health_region, crs(qpe_raster))

# Crop
qpe_crop <- crop(qpe_raster, vect(health_region_proj))

# Perform zonal mean
mean_vals <- extract(qpe_crop, vect(health_region), fun = mean, na.rm = TRUE)[,3]

rainfall <- tibble(region=health_region$region,year_to_date_in=mean_vals)

# Optional: convert to total annual mm (1 in = 25.4 mm)
rainfall$year_to_date_mm <- mean_vals * 25.4

#Save
write.csv(rainfall, data_path("processed", "simulation_inputs", "rainfall2023.csv"),
          row.names=F)
