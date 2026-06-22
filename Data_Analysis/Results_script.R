# 1. Load Libraries (Consolidated)
library(tidyverse)
library(nimble)
library(coda)
library(sf)
library(readxl)
library(viridis)

source(file.path(Sys.getenv("LEPTO_BURDEN_ROOT", "."), "R", "paths.R"))

# 2. Load Data & MCMC Results
setwd(analysis_path())
example       <- read_rds("example_inputs.rds")
chain_samples <- readRDS("Output/results_ratio1_N10.rds")
draws <- as.data.frame(as.matrix(chain_samples))
poisson_intercept <- if ("beta_0" %in% names(draws)) "beta_0" else "A"
summary(
  chain_samples[, c(
    poisson_intercept,"alpha_0","nu[1]","nu[2]","kappa[1]","kappa[2]",
    "rate_2023[1]","rate_2023[2]","rate_2023[3]",
    "rate_2023[4]","rate_2023[5]","rate_2023[6]",
    "rate_2023[7]","rate_PR_2023","rate_PR_2022",
    "log_sigma_Phi_pois","log_sigma_phi_binom",
    "pi_A","log_sigma_tau"
  )]
)

# Inverse logit function
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

# Indices of interest
idx <- c(15, 23, 24, 32)

# ---- 1) Summarize inverse-logit(alpha_0) ----
alpha_0_prob <- inv_logit(draws$alpha_0)

alpha_summary <- data.frame(
  parameter = "alpha_0 (inv_logit)",
  mean      = mean(alpha_0_prob),
  sd        = sd(alpha_0_prob),
  q2.5      = quantile(alpha_0_prob, 0.025),
  median    = quantile(alpha_0_prob, 0.5),
  q97.5     = quantile(alpha_0_prob, 0.975)
)

# ---- 2) Summarize inverse-logit(alpha_0 + phi_i * sigma) ----
results <- do.call(rbind, lapply(idx, function(i) {
  
  linear_pred <- draws$alpha_0 +
    draws[[paste0("phi_binom_raw[", i, "]")]] *
    exp(draws$log_sigma_phi_binom)
  
  prob <- inv_logit(linear_pred)
  
  data.frame(
    parameter = paste0("index_", i),
    mean      = mean(prob),
    sd        = sd(prob),
    q2.5      = quantile(prob, 0.025),
    median    = quantile(prob, 0.5),
    q97.5     = quantile(prob, 0.975)
  )
}))

# Combine everything
summary_df <- rbind(alpha_summary, results)



data          <- read_csv("Incidence_data.csv")
y <- data %>%
  # Remove "Imported" cases
  filter(!str_detect(`Health Region`, "Imported")) %>%
  # Standardize region names to ensure 2022 and 2023 match (7 regions total)
  mutate(`Health Region` = recode(`Health Region`, "Metropolitana" = "Metro")) %>%
  # Calculate new total using only Confirmed and Probable cases
  mutate(Case_Sum = `Confirmed Cases` + `Probable Cases`) %>%
  # Keep only the columns needed for the matrix
  dplyr::select(Year, `Health Region`, Case_Sum) %>%
  # Pivot the years into columns
  pivot_wider(names_from = Year, values_from = Case_Sum) %>%
  # Move Health Region names to the matrix row names
  column_to_rownames("Health Region") %>%
  # Convert to a formal matrix object
  as.matrix()
pop           <- example$data_list$pop # Ensure this is a vector of 7 values

# 2. Extract MCMC Summaries (The FIX for summ_quant)
mcmc_summ_list <- summary(chain_samples)
summ           <- mcmc_summ_list$statistics
summ_quant     <- mcmc_summ_list$quantiles

# 1. Define Colors and Names
color_2023 <- "#238B45" # Green
color_2022 <- "#2171B5" # Blue
region_names <- c("Arecibo", "Bayamón", "Caguas", "Fajardo", "Mayagüez", "Metro", "Ponce")

# 2. Build the Regional Data Frame
# We create 2023 and 2022 dataframes then bind them
df_2023 <- data.frame(
  Region   = region_names,
  Year     = "2023",
  Observed = y[, "2023"] / pop,
  Median   = summ_quant[paste0("rate_2023[", 1:7, "]"), "50%"],
  Lower    = summ_quant[paste0("rate_2023[", 1:7, "]"), "2.5%"],
  Upper    = summ_quant[paste0("rate_2023[", 1:7, "]"), "97.5%"]
)

df_2022 <- data.frame(
  Region   = region_names,
  Year     = "2022",
  Observed = y[, "2022"] / pop,
  # Note: Assuming regional 2022 rates follow the same naming convention
  Median   = summ_quant[paste0("rate_2022[", 1:7, "]"), "50%"], 
  Lower    = summ_quant[paste0("rate_2022[", 1:7, "]"), "2.5%"],
  Upper    = summ_quant[paste0("rate_2022[", 1:7, "]"), "97.5%"]
)

# 3. Build Puerto Rico Data Frame
df_island <- data.frame(
  Region   = "Puerto Rico",
  Year     = c("2023", "2022"),
  Observed = c(sum(y[, "2023"]) / sum(pop), sum(y[, "2022"]) / sum(pop)),
  Median   = c(summ_quant["rate_PR_2023", "50%"], summ_quant["rate_PR_2022", "50%"]),
  Lower    = c(summ_quant["rate_PR_2023", "2.5%"], summ_quant["rate_PR_2022", "2.5%"]),
  Upper    = c(summ_quant["rate_PR_2023", "97.5%"], summ_quant["rate_PR_2022", "97.5%"])
)

# --- 1. Data Preparation for "Broken Lines" ---
# Define how much space the median number needs (the "break" in the line)
# Adjust 'gap_size' if the numbers are too close or too far from the lines
gap_size <- max(plot_data$Upper, na.rm = TRUE) * 0.03 
overlap_limit <- max(plot_data$Upper, na.rm = TRUE) * 0.05

# --- 1. Update Factor Levels and Logic ---
region_names_sorted <- sort(region_names) # A to Z
target_order <- c("Puerto Rico", region_names_sorted)

plot_data_final <- plot_data %>%
  mutate(
    # Set the order: Puerto Rico first, then alphabetical
    Region = factor(Region, levels = target_order),
    Year = factor(Year, levels = c("2023", "2022")),
    
    # "Broken line" logic
    seg1_end = Median - gap_size,
    seg2_start = Median + gap_size,
    
    # Label logic: hide lower predictive label if it crashes into the triangle
    low_label_clean = ifelse(abs(Lower - Observed) < overlap_limit, "", round(Lower)),
    up_label_clean = round(Upper)
  )

# --- 2. Final Plotting Code ---
forest_p <- ggplot(plot_data_final, aes(x = Year, y = Median, color = Year)) +
  facet_grid(Region ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # Segmented Lines (The "Break" for the number)
  geom_segment(aes(x = Year, xend = Year, y = Lower, yend = seg1_end), 
               linewidth = 0.8) +
  geom_segment(aes(x = Year, xend = Year, y = seg2_start, yend = Upper), 
               linewidth = 0.8) +
  
  # Large End Caps
  geom_segment(aes(x = as.numeric(Year) - 0.25, xend = as.numeric(Year) + 0.25, 
                   y = Lower, yend = Lower), linewidth = 0.8) +
  geom_segment(aes(x = as.numeric(Year) - 0.25, xend = as.numeric(Year) + 0.25, 
                   y = Upper, yend = Upper), linewidth = 0.8) +
  
  # Posterior Median Number (Point Estimate)
  geom_text(aes(label = med_lab), 
            fontface = "bold", size = 4.5, show.legend = FALSE) + 
  
  # Observed Point (Triangle)
  geom_point(aes(y = Observed), color = "black", size = 2.5, shape = 17) +
  
  # Large Labels: Observed, Lower, and Upper
  geom_text(aes(y = Observed, label = obs_lab), 
            hjust = 2.2, vjust = 0.5, color = "black", size = 3.5) +
  
  geom_text(aes(y = Lower, label = low_label_clean), 
            hjust = 1.4, vjust = 0.5, size = 3.5, alpha = 0.9) +
  
  geom_text(aes(y = Upper, label = up_label_clean), 
            hjust = -0.5, vjust = 0.5, size = 3.5, alpha = 0.9) +
  
  # Orientation and Scales
  coord_flip() +
  scale_color_manual(values = c("2023" = "#238B45", "2022" = "#2171B5")) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) + 
  
  labs(x = NULL, y = "Rate per 100,000") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 11, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.6, "lines")
  )

print(forest_p)

ggsave("PR_ForestPlot_Final.png", plot = forest_p, width = 10, height = 6, dpi = 300, bg = "white")
