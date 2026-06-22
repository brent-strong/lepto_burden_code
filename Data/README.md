# Data Directory

This folder is organized by file role rather than by the script that created it.

```text
Data/
├── raw/geospatial/
│   ├── friction_surface/          # travel-time friction rasters
│   ├── health_regions/            # Puerto Rico health-region shapefile
│   ├── hospitals/                 # Puerto Rico hospital shapefile
│   └── population/                # dasymetric population raster
├── processed/simulation_inputs/   # CSV inputs used by model and simulation scripts
├── figures/                       # generated setup maps and presentation-probability figures
└── scripts/                       # data-generation and setup scripts
```

The canonical model inputs are in `Data/processed/simulation_inputs/`. Scripts should use `R/paths.R` helpers rather than hard-coded local paths.
