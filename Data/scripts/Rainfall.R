# Required libraries
library(terra)
library(sf)
library(dplyr)
library(lubridate)
library(purrr)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

# Path to your region shapefile
health_region <- st_read(data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp")) %>%
  st_transform(crs = 4326)

# Define date range
date_seq <- seq(as.Date("2020-01-01"), as.Date("2024-12-31"), by = "1 day")

# Storage lists
region_daily_list <- list()
grid_daily_list   <- list()

# Loop through each date

for(d in 1:length(date_seq)){
  
  yr <- format(date_seq[d], "%Y")
  mo <- format(date_seq[d], "%m")
  dy <- format(date_seq[d], "%d")
  ymd_str <- format(date_seq[d], "%Y%m%d")
  period <- "1day"
  fmt    <- "tif"
  
  url <- sprintf(
    "https://water.noaa.gov/resources/downloads/precip/stageIV/%s/%s/%s/nws_precip_%s_%s_pr.%s",
    yr, mo, dy, period, ymd_str, fmt
  )
  
  dest <- tempfile(fileext = paste0(".", fmt))
  tryCatch({
    download.file(url, dest, mode = "wb", quiet = TRUE)
    qpe_raster <- rast(dest)
  }, error = function(e) {
    stop("Failed to download or read raster:\n", url)
  })
  
  health_region_proj <- st_transform(health_region, crs(qpe_raster))
  qpe_crop <- crop(qpe_raster, vect(health_region_proj))
  
  mean_vals <- extract(qpe_crop, vect(health_region_proj), fun = mean, na.rm = TRUE)
  
  # Coerce to numeric to avoid prettyNum error
  rain_vals_num <- as.numeric(mean_vals[[2]])  
  
  df_region <- tibble(
    date = date_seq[d],
    region = health_region$region,
    rainfall_in = rain_vals_num,
    rainfall_mm = rain_vals_num * 25.4
  )
  
  region_daily_list[[length(region_daily_list)+1]] <- df_region
  
  # --- (2) Daily grid-cell rainfall ---
  df_grid <- as.data.frame(qpe_crop, xy = TRUE, cells = TRUE)
  
  # Dynamically find the raster value column (last one)
  
  df_grid <- df_grid %>%
    rename(cell_id = cell) %>% 
    select(cell_id,x,y) %>%
  mutate(rainfall_in=as.numeric(df_grid[,4])) %>%
  mutate(rainfall_mm = rainfall_in * 25.4)
  
df_grid$date <- date_seq[d]
  
  grid_daily_list[[length(grid_daily_list)+1]] <- df_grid
}

# --- Combine all daily region data ---
region_daily <- bind_rows(region_daily_list) %>% mutate(rainfall_in=ifelse(rainfall_in<0,0,rainfall_in),
                                                      rainfall_mm=ifelse(rainfall_mm<0,0,rainfall_mm))


# Annual totals by region
region_annual <- region_daily %>%
  group_by(region, year = year(date)) %>%
  summarise(total_rain_in = sum(rainfall_in, na.rm = TRUE),
            total_rain_mm = sum(rainfall_mm, na.rm = TRUE),
            .groups = "drop")

# --- Combine all daily grid data ---
grid_daily <- bind_rows(grid_daily_list) %>% mutate(rainfall_in=ifelse(rainfall_in<0,0,rainfall_in),
                                                    rainfall_mm=ifelse(rainfall_in<0,0,rainfall_mm))

# Save outputs
write.csv(region_daily, data_path("processed", "simulation_inputs", "rainfall_daily_by_region.csv"), row.names = FALSE)
write.csv(region_annual, data_path("processed", "simulation_inputs", "rainfall_annual_by_region.csv"), row.names = FALSE)
write.csv(grid_daily, data_path("processed", "simulation_inputs", "rainfall_daily_by_grid.csv"), row.names = FALSE)
