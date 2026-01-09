# ============================================================
# BIVARIATE RASTER CLASS DEFINITION (Suitability × Healthcare)
# ============================================================
# Each raster cell is assigned an integer value (1–9) encoding
# the joint tercile-based classification of:
#   (i) Climatic suitability
#   (ii) Healthcare access (distance to nearest facility)
#
# Classification meaning:
#
# Raster value | Suitability | Healthcare access | Interpretation
# ---------------------------------------------------------------
# 1            | Low         | Good              | Low suitability, good access
# 2            | Low         | Moderate          | Low suitability, moderate access
# 3            | Low         | Poor              | Low suitability, poor access
# 4            | Medium      | Good              | Moderate suitability, good access
# 5            | Medium      | Moderate          | Moderate suitability, moderate access
# 6            | Medium      | Poor              | Moderate suitability, poor access
# 7            | High        | Good              | High suitability, good access
# 8            | High        | Moderate          | High suitability, moderate access
# 9            | High        | Poor              | High suitability, poor access
#
# Notes:
# - Suitability and healthcare access are discretized into terciles
#   (Low / Medium / High and Good / Moderate / Poor, respectively)
#   within the population-masked domain.
# - Higher raster values indicate increasing co-location of
#   epidemiological suitability and poor healthcare access.
# - Class 9 ("High_Poor") represents priority zones where
#   high disease suitability coincides with limited healthcare access.
# ============================================================
# ============================================================
# FULL PIPELINE (UPDATED):
# Bootstrap sensitivity + COMPARABLE bivariate rasters
# ============================================================

library(terra)
library(boot)
library(ggplot2)

setwd('/Users/kagboka/Desktop/socioeeconomic_modelling_review/')

# ============================================================
# 1. LOAD RASTERS
# ============================================================

suitability_rast <- rast('composite_additive_10_MEAN.tif')
healthcare_rast  <- rast('distance_to_health.tif') / 1000  # meters → km
population_rast  <- rast('Population_Maghreb.tif')

# ============================================================
# 2. ALIGN RASTERS (STRICT)
# ============================================================

healthcare_rast <- resample(healthcare_rast, suitability_rast, method = "bilinear")
population_rast <- project(population_rast, suitability_rast, method = "bilinear")
population_rast <- resample(population_rast, suitability_rast, method = "near")

compareGeom(suitability_rast, healthcare_rast, stopOnError = TRUE)
compareGeom(suitability_rast, population_rast, stopOnError = TRUE)

# ============================================================
# 3. STACK & DATAFRAME (ALL PIXELS)
# ============================================================

r_stack <- c(suitability_rast, healthcare_rast, population_rast)
names(r_stack) <- c("Suitability", "Healthcare", "Population")

df_all <- as.data.frame(r_stack, xy = TRUE, na.rm = TRUE)

# ============================================================
# 4. BOOTSTRAP SENSITIVITY ANALYSIS
# ============================================================

thresholds <- c(1, 5, 10, 25, 50, 100)

results <- data.frame(
  Threshold = thresholds,
  N = NA,
  Pearson_r = NA,
  Pearson_CI_low = NA,
  Pearson_CI_high = NA,
  Spearman_rho = NA,
  Spearman_p = NA
)

boot_pearson <- function(data, indices) {
  d <- data[indices, ]
  cor(d$Suitability, d$Healthcare, method = "pearson")
}

for (i in seq_along(thresholds)) {
  
  thr <- thresholds[i]
  df_thr <- df_all[df_all$Population >= thr, ]
  
  results$N[i] <- nrow(df_thr)
  if (nrow(df_thr) < 50) next
  
  set.seed(123)
  boot_out <- boot(df_thr, statistic = boot_pearson, R = 500)
  ci <- boot.ci(boot_out, type = "perc")
  
  results$Pearson_r[i] <- cor(df_thr$Suitability, df_thr$Healthcare)
  results$Pearson_CI_low[i]  <- ci$perc[4]
  results$Pearson_CI_high[i] <- ci$perc[5]
  
  sp <- cor.test(df_thr$Suitability, df_thr$Healthcare,
                 method = "spearman", exact = FALSE)
  
  results$Spearman_rho[i] <- as.numeric(sp$estimate)
  results$Spearman_p[i]   <- sp$p.value
}

write.csv(results, "Bootstrap_Sensitivity_Table.csv", row.names = FALSE)
print(results)

# ============================================================
# 5. FIXED GLOBAL BREAKS (CRITICAL CHANGE)
# ============================================================
# These breaks DEFINE the meaning of Low / Medium / High
# and are reused for ALL thresholds.

mask_ref <- population_rast >= 1

suit_ref <- mask(suitability_rast, mask_ref)
acc_ref  <- mask(healthcare_rast,  mask_ref)

suit_breaks <- quantile(values(suit_ref),
                        probs = c(0, 1/3, 2/3, 1),
                        na.rm = TRUE)

acc_breaks  <- quantile(values(acc_ref),
                        probs = c(0, 1/3, 2/3, 1),
                        na.rm = TRUE)

# ============================================================
# 6. BIVARIATE COLOR PALETTE + CLASS MEANING
# ============================================================

# Raster values and interpretation:
# 1 = Low suitability, good access
# 2 = Low suitability, moderate access
# 3 = Low suitability, poor access
# 4 = Medium suitability, good access
# 5 = Medium suitability, moderate access
# 6 = Medium suitability, poor access
# 7 = High suitability, good access
# 8 = High suitability, moderate access
# 9 = High suitability, poor access (PRIORITY ZONES)

bivar_classes <- c(
  "Low_Good", "Low_Moderate", "Low_Poor",
  "Medium_Good", "Medium_Moderate", "Medium_Poor",
  "High_Good", "High_Moderate", "High_Poor"
)

bivar_colors <- c(
  "#e8e8e8", "#ace4e4", "#5ac8c8",
  "#dfb0d6", "#a5add3", "#5698b9",
  "#be64ac", "#8c62aa", "#3b4994"
)
names(bivar_colors) <- bivar_classes

# ============================================================
# 7. LOOP: CREATE & SAVE COMPARABLE BIVARIATE RASTERS
# ============================================================

for (thr in thresholds) {
  
  message("Processing population threshold ≥ ", thr)
  
  mask_thr <- population_rast >= thr
  suit_thr <- mask(suitability_rast, mask_thr)
  acc_thr  <- mask(healthcare_rast,  mask_thr)
  
  df_bi <- as.data.frame(c(suit_thr, acc_thr), xy = TRUE, na.rm = TRUE)
  names(df_bi) <- c("x", "y", "Suitability", "Healthcare")
  
  # FIXED classification (NO recomputation of quantiles)
  df_bi$Suit_cat <- cut(
    df_bi$Suitability,
    breaks = suit_breaks,
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE
  )
  
  df_bi$Access_cat <- cut(
    df_bi$Healthcare,
    breaks = acc_breaks,
    labels = c("Good", "Moderate", "Poor"),
    include.lowest = TRUE
  )
  
  df_bi$bivar_class <- paste(df_bi$Suit_cat, df_bi$Access_cat, sep = "_")
  df_bi$bivar_id <- match(df_bi$bivar_class, bivar_classes)
  
  # Rasterize
  bivar_rast <- rast(suit_thr)
  bivar_rast <- rasterize(df_bi[, c("x", "y")],
                          bivar_rast,
                          df_bi$bivar_id)
  
  levels(bivar_rast) <- data.frame(
    value = 1:9,
    class = bivar_classes
  )
  
  writeRaster(
    bivar_rast,
    paste0("Bivariate_Suitability_Healthcare_pop", thr, ".tif"),
    overwrite = TRUE
  )
  
  # ============================================================
  # 8. MAP FROM SAVED RASTER
  # ============================================================
  
  df_plot <- as.data.frame(bivar_rast, xy = TRUE, na.rm = TRUE)
  names(df_plot) <- c("x", "y", "id")
  df_plot$class <- bivar_classes[df_plot$id]
  
  p <- ggplot(df_plot) +
    geom_raster(aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = bivar_colors, guide = "none") +
    coord_equal() +
    theme_void() +
    labs(
      title = "Bivariate suitability × healthcare access",
      subtitle = paste("Population density ≥", thr, "persons/km²")
    )
  
  ggsave(
    paste0("Bivariate_Map_pop", thr, ".png"),
    plot = p,
    width = 9,
    height = 4.5,
    dpi = 600
  )
}

# ============================================================
# DONE – Maps are now COMPARABLE across thresholds
# ============================================================