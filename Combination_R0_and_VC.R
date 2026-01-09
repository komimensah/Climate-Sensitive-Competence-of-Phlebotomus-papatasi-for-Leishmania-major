

# ==============================================================================
# Composite Risk Index – Quarterly Pairing (WITH RODENT BACKGROUND CONSTRAINT)
# ==============================================================================

rm(list = ls()); gc()
library(terra)

setwd("/Users/kagboka/Desktop/socioeeconomic_modelling_review/Composite risk_index/")

lambda_dir <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/RI_quarterly"
vc_dir     <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/VC_quarterly"

lambda_files <- sort(list.files(lambda_dir, "\\.tif$", full.names = TRUE))
vc_files     <- sort(list.files(vc_dir,     "\\.tif$", full.names = TRUE))

stopifnot(length(lambda_files) == 4, length(vc_files) == 4)

quarters <- c("Q1","Q2","Q3","Q4")

weights <- list(
  "10_90" = c(0.10,0.90),
  "20_80" = c(0.20,0.80),
  "30_70" = c(0.30,0.70),
  "40_60" = c(0.40,0.60),
  "50_50" = c(0.50,0.50),
  "60_40" = c(0.60,0.40),
  "70_30" = c(0.70,0.30),
  "80_20" = c(0.80,0.20),
  "90_10" = c(0.90,0.10)
)

# ------------------------------------------------------------------
# NEW: Load rodent suitability BACKGROUND raster
# ------------------------------------------------------------------
rodent_bg_path <- '/Users/kagboka/Desktop/socioeeconomic_modelling_review/final_ensemble.tif'
rodent_bg <- rast(rodent_bg_path)

for (i in 1:4) {
  
  message("▶ Processing ", quarters[i])
  
  lambda_r <- rast(lambda_files[i])
  vc_r     <- rast(vc_files[i])
  
  # ------------------------------------------------------------
  # CRS alignment ONLY (as before)
  # ------------------------------------------------------------
  if (!same.crs(lambda_r, vc_r)) {
    vc_r <- project(vc_r, crs(lambda_r))
  }
  vc_r <- resample(vc_r, lambda_r, method = "bilinear")
  
  # ------------------------------------------------------------
  # NEW: align rodent background to lambda grid
  # ------------------------------------------------------------
  if (!same.crs(rodent_bg, lambda_r)) {
    rodent_bg <- project(rodent_bg, crs(lambda_r))
  }
  rodent_bg_q <- resample(rodent_bg, lambda_r, method = "bilinear")
  
  # ------------------------------------------------------------
  # VALID MASK (UNCHANGED)
  # ------------------------------------------------------------
  valid <- is.finite(lambda_r) & is.finite(vc_r) & is.finite(rodent_bg_q)
  
  # ------------------------------------------------------------
  # Geometric
  # ------------------------------------------------------------
  gm <- ifel(valid, sqrt(lambda_r * vc_r) * rodent_bg_q, NA)
  writeRaster(gm, paste0("composite_geometric_", quarters[i], ".tif"), overwrite=TRUE)
  
  # ------------------------------------------------------------
  # Harmonic
  # ------------------------------------------------------------
  hm <- ifel(valid & (lambda_r + vc_r) > 0,
             (2 * lambda_r * vc_r / (lambda_r + vc_r)) * rodent_bg_q,
             NA)
  writeRaster(hm, paste0("composite_harmonic_", quarters[i], ".tif"), overwrite=TRUE)
  
  # ------------------------------------------------------------
  # Additive (UNCHANGED weights)
  # ------------------------------------------------------------
  for (w in names(weights)) {
    
    a <- weights[[w]][1]
    b <- weights[[w]][2]
    
    am <- ifel(valid, (a * lambda_r + b * vc_r) * rodent_bg_q, NA)
    
    writeRaster(
      am,
      paste0("composite_additive_", w, "_", quarters[i], ".tif"),
      overwrite = TRUE
    )
  }
  
  # ------------------------------------------------------------
  # Multiplicative
  # ------------------------------------------------------------
  mul <- ifel(valid, (lambda_r * vc_r) * rodent_bg_q, NA)
  writeRaster(mul, paste0("composite_multiplicative_", quarters[i], ".tif"), overwrite=TRUE)
}

cat("\n✔ Composite rasters recomputed WITH rodent background constraint\n")

# ------------------------------------------------------------------------------
# MEAN COMPOSITES ACROSS QUARTERS (UNCHANGED)
# ------------------------------------------------------------------------------

mean_dir <- file.path(getwd(), "composite_mean")
dir.create(mean_dir, showWarnings = FALSE, recursive = TRUE)

quarterly_files <- list.files(
  getwd(),
  pattern = "^composite_.*_Q[1-4]\\.tif$",
  full.names = TRUE
)

get_category <- function(f) {
  nm <- basename(f)
  nm <- sub("_Q[1-4]\\.tif$", "", nm)
  nm <- sub("_[0-9]+$", "", nm)
  nm
}

categories <- unique(vapply(quarterly_files, get_category, character(1)))

for (category in categories) {
  
  message("▶ Computing MEAN for ", category)
  
  cat_files <- quarterly_files[
    vapply(quarterly_files, function(f) {
      get_category(f) == category
    }, logical(1))
  ]
  
  if (length(cat_files) != 4) {
    warning("Skipping ", category,
            ": expected 4 quarters, found ", length(cat_files))
    next
  }
  
  r_stack <- rast(cat_files)
  r_mean  <- app(r_stack, mean, na.rm = TRUE)
  
  writeRaster(
    r_mean,
    file.path(mean_dir, paste0(category, "_MEAN.tif")),
    overwrite = TRUE
  )
}

cat("\n✔ Mean composite rasters created in 'composite_mean/'\n")

