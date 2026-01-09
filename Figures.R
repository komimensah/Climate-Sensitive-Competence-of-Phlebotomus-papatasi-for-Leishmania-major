

library(ggplot2)
setwd('/Users/kagboka/Desktop/Lesh0.0/')
# === 1. Load rasters ===
suitability_rast <- rast('/Users/kagboka/Desktop/Lesh0.0/composite_additive_10_90.tif')      # Your suitability score raster
healthcare_rast  <- rast('/Users/kagboka/Desktop/Lesh0.0/distance_to_health.tif')            # Your healthcare access raster

# === Convert healthcare raster from meters to kilometers ===
healthcare_rast <- healthcare_rast / 1000
# === 2. Align healthcare raster to suitability raster ===
healthcare_rast <- resample(healthcare_rast, suitability_rast, method = "bilinear")

# === 3. Stack and convert to dataframe ===
r_stack <- c(suitability_rast, healthcare_rast)
names(r_stack) <- c("Suitability", "Healthcare")

df <- as.data.frame(r_stack, xy = TRUE, na.rm = TRUE)

# === 4. Sample 10,000 points ===
set.seed(123)  # for reproducibility
df_sample <- df[sample(nrow(df), min(10000, nrow(df))), ]

# === 5. Compute correlation ===
r_val <- cor(df_sample$Suitability, df_sample$Healthcare, use = "complete.obs")

library(ggplot2)

# Create the scatter plot with legend labels
p <- ggplot(df_sample, aes(x = Suitability, y = Healthcare)) +
  geom_point(aes(color = "Suitability Points"), alpha = 0.3, size = 0.7) +
  geom_smooth(aes(color = "Linear Trend"), method = "lm", se = FALSE, linewidth = 1.3) +
  scale_color_manual(name = "Legend", values = c("Suitability Points" = "steelblue", "Linear Trend" = "darkred")) +
  labs(
    title = "Suitability vs Health Care Access",
    subtitle = "Sample of 10,000 spatial points",
    x = "Suitability Index (λₘₐₓ)",
    y = "Distance to Health Care (km)"
  ) +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1.3,
           label = paste0("r = ", round(r_val, 3)),
           size = 5.2, color = "black") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )
print(p)
# Save high-resolution plot
ggsave("suitability_vs_healthcare.png", plot = p, width = 8, height = 6, dpi = 600)
#####################
library(terra)


# === 1. Load all rasters from folder ===
folder_path <- "/Users/kagboka/Desktop/Lesh0.0/moi/"
raster_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)

# Extract base names for labeling
raster_names <- tools::file_path_sans_ext(basename(raster_files))

# Load and stack rasters with names
r_stack <- rast(raster_files)
names(r_stack) <- raster_names  # Assign meaningful layer names

# === 2. Load coordinate CSV ===
coords_df <- read.csv("Final data revised_2.csv")  # Should include 'lon' and 'lat'

# === 3. Create SpatVector in WGS84 (common GPS) ===
points_vect <- vect(coords_df, geom = c("lon", "lat"), crs = "EPSG:4326")

# === 4. Reproject to raster CRS
points_vect_proj <- project(points_vect, crs(r_stack))

# === 5. Extract raster values at point locations
extracted_vals <- extract(r_stack, points_vect_proj)

# === 6. Merge with original coordinates
results <- cbind(coords_df, extracted_vals[,-1])  # Drop ID column
# === 7. Load required libraries for plotting ===
# Install if needed
# install.packages(c("ggplot2","ggpubr","dplyr","tidyr"))


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# 1. Prepare long format
long <- results %>%
  select(delta_SSP245, delta_SSP585) %>%
  pivot_longer(cols = everything(),
               names_to = "Scenario",
               values_to = "Delta") %>%
  mutate(Scenario = factor(case_when(
    Scenario == "delta_SSP245" ~ "SSP2-4.5",
    Scenario == "delta_SSP585" ~ "SSP5-8.5"
  ), levels = c("SSP2-4.5", "SSP5-8.5")))

# 2. Plot
p <- ggplot(long, aes(x = Scenario, y = Delta, fill = Scenario)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               fill = "white", color = "black", size = 0.8) +  # <-- bold box border
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, fill = "black") +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.1, color = "grey30") +
  scale_fill_manual(values = c("SSP2-4.5" = "#4C72B0", "SSP5-8.5" = "#C44E52")) +
  labs(title = "",
       subtitle = "",
       y = "Δ Suitability",
       x = NULL,
       fill = "Scenario",
       caption = "") +
  theme_pubr() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

# 3. Print
print(p)

# 4. Save
ggsave("suitability_violin_boxplot.png", plot = p, width = 7.5, height = 5.5, dpi = 300)
ggsave("suitability_violin_boxplot.pdf", plot = p, width = 7.5, height = 5.5)
print(p)
# Ensure required variables exist and have no NAs
paired_data <- na.omit(results[, c("delta_SSP245", "delta_SSP585")])

# Run Wilcoxon signed-rank test (non-parametric paired test)
wilcox_test <- wilcox.test(paired_data$delta_SSP585, 
                           paired_data$delta_SSP245, 
                           paired = TRUE, 
                           alternative = "two.sided")

# Print results
print(wilcox_test)
summary (results)
write.csv(results, "climate change analysis.csv", row.names = FALSE)


###########################Figure 6-7

##########################################
library(terra)
library(boot)

setwd('/Users/kagboka/Desktop/socioeeconomic_modelling_review/')

# Load rasters
suitability_rast <- rast('/Users/kagboka/Desktop/socioeeconomic_modelling_review/Composite risk_index_present_nor/composite_mean/composite_additive_10_MEAN.tif')
healthcare_rast  <- rast('distance_to_health.tif') / 1000  # convert to km
population_rast  <- rast('Population_Maghreb.tif')

# Align rasters
healthcare_rast <- resample(healthcare_rast, suitability_rast, method="bilinear")
population_rast <- project(population_rast, suitability_rast, method="bilinear")
population_rast <- resample(population_rast, suitability_rast, method="bilinear")

# Create combined data frame
r_stack <- c(suitability_rast, healthcare_rast, population_rast)
names(r_stack) <- c("Suitability", "Healthcare", "Population")
df_all <- as.data.frame(r_stack, xy=TRUE, na.rm=TRUE)

# Define thresholds to evaluate
thresholds <- c(1, 5, 10, 25, 50, 100)
results <- data.frame(
  Threshold=thresholds,
  N=NA,
  Pearson_r=NA,
  Pearson_CI_low=NA,
  Pearson_CI_high=NA,
  Spearman_rho=NA,
  Spearman_p=NA
)

# Function to compute pearson r for bootstrap
boot_pearson <- function(data, indices) {
  d <- data[indices,]
  return(cor(d$Suitability, d$Healthcare, method="pearson"))
}

# Loop through thresholds
for (i in seq_along(thresholds)) {
  thr <- thresholds[i]
  
  df_thresh <- df_all[df_all$Population >= thr, ]
  results$N[i] <- nrow(df_thresh)
  
  if (nrow(df_thresh) < 30) {
    # Too few for reliable correlation
    next
  }
  
  # Pearson with bootstrap confidence intervals
  set.seed(123)
  boot_out <- boot(df_thresh, statistic=boot_pearson, R=500)
  ci <- boot.ci(boot_out, type="perc")
  
  pearson_r <- cor(df_thresh$Suitability, df_thresh$Healthcare, method="pearson")
  results$Pearson_r[i] <- pearson_r
  results$Pearson_CI_low[i]  <- ci$perc[4]
  results$Pearson_CI_high[i] <- ci$perc[5]
  
  # Spearman
  spearman_test <- cor.test(df_thresh$Suitability, df_thresh$Healthcare,
                            method="spearman", exact=FALSE)
  
  results$Spearman_rho[i] <- spearman_test$estimate
  results$Spearman_p[i]   <- spearman_test$p.value
}

print(results)
library(ggplot2)

# Prepare data for plotting
plot_df <- results
plot_df$Threshold <- factor(plot_df$Threshold, levels = thresholds)

# ===============================
# High-standard sensitivity plot
# ===============================
p_sens <- ggplot(plot_df, aes(x = Threshold)) +
  
  # Pearson r with CI
  geom_point(aes(y = Pearson_r), size = 3.2, color = "black") +
  geom_errorbar(
    aes(ymin = Pearson_CI_low, ymax = Pearson_CI_high),
    width = 0.15,
    linewidth = 0.9
  ) +
  geom_line(aes(y = Pearson_r, group = 1),
            linewidth = 1.1, color = "black") +
  
  # Spearman rho (scaled for visual comparability)
  geom_point(
    aes(y = Spearman_rho),
    size = 3,
    shape = 21,
    fill = "white",
    color = "darkred"
  ) +
  geom_line(
    aes(y = Spearman_rho, group = 1),
    linewidth = 1,
    linetype = "dashed",
    color = "darkred"
  ) +
  
  labs(
    x = "Population density threshold (persons per km²)",
    y = "Correlation coefficient",
    title = "",
    subtitle = "Pearson r (solid, with 95% CI) and Spearman ρ (dashed) across population thresholds"
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 13),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8)
  )

print(p_sens)

# Save figure (high resolution)
ggsave(
  "Figure_Sensitivity_PopulationThresholds.png",
  plot = p_sens,
  width = 7.5,
  height = 5.5,
  dpi = 600
)
