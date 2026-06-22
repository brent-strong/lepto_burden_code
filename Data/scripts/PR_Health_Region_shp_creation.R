# Load required packages
library(tigris)
library(tidyverse)
library(sf)
library(ggplot2)
options(tigris_use_cache = TRUE)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

# Step 1: Download municipality shapefile for Puerto Rico
pr_munis <- counties(state = "PR", cb = TRUE, year = 2020) %>%
  st_transform(crs = 4326)  # Projected CRS (Puerto Rico State Plane)

# Step 2: Lookup table for health regions
region_lookup <- tibble::tribble(
  ~NAME,         ~region,
  "Arecibo",     "Arecibo",
  "Barceloneta", "Arecibo",
  "Camuy",       "Arecibo",
  "Ciales",      "Arecibo",
  "Dorado",      "Arecibo",
  "Florida",     "Arecibo",
  "Hatillo",     "Arecibo",
  "Lares",     "Arecibo",
  "Manatí",      "Arecibo",
  "Morovis",     "Arecibo",
  "Quebradillas", "Arecibo",
  "Utuado", "Arecibo",
  "Vega Alta",   "Arecibo",
  "Vega Baja",   "Arecibo",
  
  "Barranquitas","Bayamón",
  "Cataño",      "Bayamón",
  "Comerío",     "Bayamón",
  "Corozal",     "Bayamón",
  "Naranjito",   "Bayamón",
  "Orocovis",    "Bayamón",
  "Toa Alta",    "Bayamón",
  "Toa Baja",    "Bayamón",
  "Bayamón","Bayamón",
  
  "Aguas Buenas","Caguas",
  "Aibonito",    "Caguas",
  "Caguas",      "Caguas",
  "Cayey",       "Caguas",
  "Cidra",       "Caguas",
  "Gurabo",      "Caguas",
  "Humacao",     "Caguas",
  "Juncos",      "Caguas",
  "Las Piedras", "Caguas",
  "Maunabo",     "Caguas",
  "Naguabo",     "Caguas",
  "San Lorenzo", "Caguas",
  "Yabucoa",     "Caguas",
  
  "Ceiba",       "Fajardo",
  "Fajardo",     "Fajardo",
  "Luquillo",    "Fajardo",
  "Río Grande",  "Fajardo",
  "Vieques",     "Fajardo",
  "Culebra",     "Fajardo",
  
  "Aguada",      "Mayagüez",
  "Aguadilla",   "Mayagüez",
  "Añasco",      "Mayagüez",
  "Cabo Rojo",   "Mayagüez",
  "Hormigueros", "Mayagüez",
  "Isabela",     "Mayagüez",
  "Las Marías",  "Mayagüez",
  "Lajas",  "Mayagüez",
  "Maricao",     "Mayagüez",
  "Mayagüez",    "Mayagüez",
  "Moca",        "Mayagüez",
  "Rincón",      "Mayagüez",
  "San Germán",  "Mayagüez",
  "San Sebastián", "Mayagüez",
  "Sabana Grande", "Mayagüez",
  
  "Guaynabo",    "Metro",
  "San Juan",    "Metro",
  "Carolina",    "Metro",
  "Trujillo Alto","Metro",
  "Loíza",       "Metro",
  "Canóvanas",   "Metro",
  
  "Adjuntas",    "Ponce",
  "Arroyo",      "Ponce",
  "Coamo",       "Ponce",
  "Guánica",     "Ponce",
  "Guayama",     "Ponce",
  "Guayanilla",  "Ponce",
  "Jayuya",      "Ponce",
  "Juana Díaz",  "Ponce",
  "Peñuelas",    "Ponce",
  "Ponce",       "Ponce",
  "Salinas",     "Ponce",
  "Santa Isabel","Ponce",
  "Villalba",    "Ponce",
  "Yauco",       "Ponce",
  "Patillas",    "Ponce",
)

# Step 3: Join shapefile with region info
pr_munis_regions <- pr_munis %>%
  left_join(region_lookup, by = "NAME") %>%
  filter(!is.na(region))  # Keep only those in a health region

# Step 4: Dissolve municipalities by health region
health_regions <- pr_munis_regions %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_as_sf()

bbox_no_mona <- st_bbox(c(xmin = -67.3, ymin = 17.8, xmax = -65.2, ymax = 18.6), 
                        crs = st_crs(health_regions))
regions_graph <- st_crop(health_regions, bbox_no_mona)

#Plot
plot_regions <- ggplot(regions_graph) +
  geom_sf(aes(fill = region), color = "white", size = 0.3) +
  scale_fill_viridis_d(name = "Health Region", option = "D") +
  geom_sf_text(data = regions_graph, aes(label = region), size = 3, color = "white") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) + 
  labs(
    x = "Longitude", 
    y = "Latitude"
  )

# Step 5: Save the result to a shapefile
st_write(health_regions, data_path("raw", "geospatial", "health_regions", "pr_health_regions.shp"))


ggsave(filename = data_path("figures", "health_regions.png"), 
       plot = plot_regions)


#Create plot with hospitals 

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

hospital_coords_sf <- st_as_sf(hospitals, coords = c("longitude", "latitude"), crs = crs(health_regions))

plot_regions_hospitals<-ggplot(regions_graph) +
  geom_sf(aes(fill = region), color = "white", size = 0.3) +
  scale_fill_viridis_d(name = "Health Region", option = "D") +
  geom_sf_text(data = regions_graph, aes(label = region), size = 3, color = "white") +
  geom_sf(data = hospital_coords_sf, color = "red", size = 2, shape = 17) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) + labs(
    x = "Longitude", 
    y = "Latitude", 
  )

ggsave(filename = data_path("figures", "health_regions_hospitals.png"), 
       plot = plot_regions_hospitals)

#Plot with only the four study hospitals

hospital_coords_sf <- st_as_sf(hospitals, coords = c("longitude", "latitude"), crs = crs(health_regions))
hospital4 <- hospital_coords_sf[hospital_coords_sf$ID %in% c("15","2512","14","23"), ]

plot_regions_hospitals_four<-ggplot(regions_graph) +
  geom_sf(aes(fill = region), color = "white", size = 0.3) +
  scale_fill_viridis_d(name = "Health Region", option = "D") +
  geom_sf_text(data = regions_graph, aes(label = region), size = 3, color = "white") +
  geom_sf(data = hospital4, aes(size = BEDS),color = "red", shape = 17) +
  scale_size_continuous(range = c(2, 4), name = "Number of Beds") +
  coord_sf(expand = FALSE) + 
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) + labs(
    x = "Longitude", 
    y = "Latitude", 
  )

ggsave(filename = data_path("figures", "health_regions_hospitals_four.png"), 
       plot = plot_regions_hospitals_four)


##New plot with all hospitals plus 4 study hospitals differentiated

library(tidyverse)
library(sf)

# 1. Prepare the data with a grouping column
study_ids <- c("15", "2512", "14", "23")

hospital_coords_sf <- hospital_coords_sf %>%
  mutate(plot_group = ifelse(ID %in% c("15", "2512", "14", "23"), ID, "Other")) %>%
  mutate(plot_group = factor(plot_group, levels = c("Other", "15", "2512", "14", "23"))) %>%
  arrange(desc(plot_group == "Other"))

plot_rate_hospitals <- ggplot(regions_graph_rates) +
  geom_sf(aes(fill = observed_rate), color = "white", size = 0.3, alpha = 0.9) +
  scale_fill_gradient(
    name = "Observed rate",
    low = "#fee5d9",
    high = "#a50f15"
  ) +
  
  geom_sf_text(
    data = regions_graph_rates,
    aes(label = region),
    size = 2.5,
    color = "black",
    fontface = "bold",
    check_overlap = TRUE
  ) +
  
  geom_sf(
    data = hospital_coords_sf %>% filter(plot_group == "Other"),
    aes(color = plot_group),
    shape = 17,
    size = 3,
    show.legend = FALSE
  ) +
  
  geom_sf(
    data = hospital_coords_sf %>% filter(plot_group != "Other"),
    aes(color = plot_group),
    shape = 17,
    size = 3,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = c(
      "Other" = "gray40",
      "15"    = "cyan",
      "2512"  = "yellow",
      "14"    = "red",
      "23"    = "purple"
    )
  ) +
  
  coord_sf(datum = NA) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  ) +
  labs(x = NULL, y = NULL)

library(magick)

plot_rate_hospitals_cropped <- plot_rate_hospitals +
  coord_sf(datum = NA, expand = FALSE) +
  theme(
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave(
  filename = data_path("figures", "vibrant_hospitals_map_uncropped.png"),
  plot = plot_rate_hospitals_cropped,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

image_read(data_path("figures", "vibrant_hospitals_map_uncropped.png")) %>%
  image_trim(fuzz = 2) %>%
  image_border(color = "white", geometry = "30x30") %>%
  image_write(data_path("figures", "vibrant_hospitals_map.png"))
