# lepto_burden

Code and supporting files for the paper **"A Bayesian Spatiotemporal Model to Estimate Disease Burden Using Hospital-Based Active Surveillance"**.

The repository contains the data setup scripts, simulation study code, and real-data analysis code for estimating hospital-presenting leptospirosis burden in Puerto Rico while accounting for under-capture, spatial misalignment, and imperfect diagnostic testing.

## Repository Layout

```text
.
├── Data/
│   ├── raw/geospatial/              # source shapefiles, rasters, and friction surfaces
│   ├── processed/simulation_inputs/ # CSV inputs consumed by model and simulation scripts
│   ├── figures/                     # maps and data-setup figures
│   └── scripts/                     # data preparation and figure-generation scripts
├── Data_Analysis/                   # real-data analysis scripts and manuscript figures
├── Simulation/                      # simulation study scripts, helper functions, and processors
├── R/paths.R                        # repo-relative path helpers
├── Background/                      # background papers and notes
└── Archive/                         # older exploratory scripts retained for reference
```

`Data/scripts/archive/` and `Data_Analysis/Archive/` are historical references and may retain older paths or notation. Use the non-archive scripts for the reproducible workflow.

## Setup

Run scripts from the repository root:

```r
source("R/paths.R")
```

On a cluster or from another working directory, set:

```bash
export LEPTO_BURDEN_ROOT=/path/to/lepto_burden
```

Main R package dependencies include `tidyverse`, `nimble`, `coda`, `sf`, `terra`, `gdistance`, `exactextractr`, `spdep`, `mvtnorm`, `readxl`, `magick`, `patchwork`, `latex2exp`, and `tigris`.

## Workflow

1. **Prepare geospatial inputs**
   - `Data/scripts/PR_Health_Region_shp_creation.R` builds the health-region shapefile and maps.
   - `Data/scripts/Presentation_Probabilities_Inhomogenous.R` creates:
     - `Data/processed/simulation_inputs/region_population.csv`
     - `Data/processed/simulation_inputs/health_regions_probs.csv`
     - `Data/processed/simulation_inputs/cell_probabilities.csv`
   - `Data/scripts/Rainfall.R` regenerates rainfall summaries.

2. **Run the real-data model**
   - `Data_Analysis/Analysis_code.R` fits the main data-analysis model.
   - `Data_Analysis/surveillance_model_inputs.csv` contains the aggregated active-surveillance inputs used by the model. It replaces the private patient-level surveillance workbook and includes only health-region/site counts plus diagnostic-result pattern counts.
   - `Data_Analysis/Results_script.R` summarizes posterior draws and produces the main forest plot.
   - `Data_Analysis/Results_script_sensitivity.R` and `Data_Analysis/Results_script_sensitivity2.R` summarize sensitivity analyses.

3. **Run simulations**
   - `Simulation/Simulation_script.R` runs the main simulation study from Section 4, varying the baseline passive capture probability (`mean_prob`) and the prior standard deviation for `alpha_0` (`sd_alpha_0`). Outputs are written to `Simulation/Simulation_output/`.
   - `Simulation/Simulation_script_sensitivity.R` runs the added-hospitals/test-prior sensitivity simulation, varying the number of active-surveillance hospitals and the prior strength for sensitivity/specificity. Outputs are written to `Simulation/Simulation_sensitivity/`.
   - `Simulation/Simulation_processing.Rmd` processes the main simulation outputs.
   - `Simulation/Simulation_processing_sensitivity.Rmd` processes the added-hospitals/test-prior sensitivity outputs.

4. **Main CRC sandbox**
   - `Example_model_script.R` contains an example capture-recapture model-fitting sandbox using the same data inputs and manuscript-aligned notation.

Generated output directories such as `Simulation/Simulation_output`, `Simulation/Simulation_sensitivity`, `Data_Analysis/Output`, and `Data_Analysis/Sensitivity_Analysis` are ignored by Git.

## Notation Crosswalk

| Paper notation | Code name | Meaning |
| --- | --- | --- |
| `beta_0` | `beta_0` | Poisson log-rate intercept |
| `beta_1` | `beta_1` | Rainfall coefficient in simulation models |
| `Phi_s` | `Phi_pois`, `Phi_pois_raw` | spatial ICAR effect for disease rates |
| `tau_s` | `tau`, `tau_raw` | unstructured region-level heterogeneity |
| `alpha_0` | `alpha_0` | baseline logit passive case-capture probability |
| `phi_h` | `phi_binom`, `phi_binom_raw` | hospital-level capture-probability effect |
| `pi^P_h` | `pi_P` / `piP` | passive surveillance case-capture probability |
| `pi^A` | `pi_A` / `piA` | active surveillance capture/testing probability among diseased patients |
| `pi^C_{h,s}` | `prob_pres`, `prob_pres_s`, `prob_pres_ns` | hospital presentation probabilities |
| `nu` | `nu` | diagnostic test sensitivity |
| `kappa` | `kappa` | diagnostic test specificity |
| `n^P` | `n_p` | passively captured disease cases |
| `n^{PA}` | `n_pa` | cases captured by both passive and active surveillance |
| `n^T` | `n_t` | patients tested under active surveillance |
| `Y` | `Y` | padded active-surveillance test-result array |

Older prototype scripts under `Archive/` may still contain scratch notation. Maintained scripts and saved `.rds` outputs use the names in this crosswalk.
