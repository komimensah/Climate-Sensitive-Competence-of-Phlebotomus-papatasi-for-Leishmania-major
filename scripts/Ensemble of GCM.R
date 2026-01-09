# ==============================================================================
# Multi-Model Temperature Ensemble: Compute taverage (robust to extents)
# ==============================================================================
library(terra)
library(sf)
library(viridis)

# ------------------------------------------------------------------------------
# 0) DEFINE INPUT RASTER STACK FILES (one per GCM)
# ------------------------------------------------------------------------------
tmin_gcm_stacks <- c(
  '/Users/kagboka/Desktop/socioeeconomic_modelling_review/wc2.1_2.5m_tmin_EC-Earth3-Veg_ssp585_2041-2060.tif',
  '/Users/kagboka/Desktop/socioeeconomic_modelling_review/wc2.1_2.5m_tmin_IPSL-CM6A-LR_ssp585_2021-2040.tif',
  '/Users/kagboka/Desktop/socioeeconomic_modelling_review/wc2.1_2.5m_tmin_UKESM1-0-LL_ssp585_2021-2040.tif'
)

tmax_gcm_stacks <- c(
  '/Users/kagboka/Desktop/socioeeconomic_modelling_review/wc2.1_2.5m_tmax_EC-Earth3-Veg_ssp585_2041-2060.tif',
  '/Users/kagboka/Desktop/socioeeconomic_modelling_review/wc2.1_2.5m_tmax_IPSL-CM6A-LR_ssp585_2021-2040.tif',
  '/Users/kagboka/Desktop/socioeeconomic_modelling_review/wc2.1_2.5m_tmax_UKESM1-0-LL_ssp585_2021-2040.tif'
)

stopifnot(length(tmin_gcm_stacks) == length(tmax_gcm_stacks))

# ------------------------------------------------------------------------------
# 1) Load study area
# ------------------------------------------------------------------------------
study_shp_path <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/Study_area/study area dissolved_B.shp"
study_area <- st_read(study_shp_path, quiet = TRUE)

# ------------------------------------------------------------------------------
# 2) Load monthly stacks
# ------------------------------------------------------------------------------
load_monthly_stack <- function(file_path) {
  r <- rast(file_path)
  stopifnot(nlyr(r) == 12)
  names(r) <- month.abb
  r
}

tmin_stacks <- lapply(tmin_gcm_stacks, load_monthly_stack)
tmax_stacks <- lapply(tmax_gcm_stacks, load_monthly_stack)
tmin_all <- do.call(c, tmin_stacks)
tmax_all <- do.call(c, tmax_stacks)

# ------------------------------------------------------------------------------
# 3) Compute monthly ensemble means
# ------------------------------------------------------------------------------
tmin_ensemble <- lapply(month.abb, function(mon) {
  layers <- tmin_all[[grep(mon, names(tmin_all))]]
  app(layers, mean, na.rm = TRUE)
})
tmin_ensemble <- rast(tmin_ensemble)
names(tmin_ensemble) <- paste0("tmin_", month.abb)

tmax_ensemble <- lapply(month.abb, function(mon) {
  layers <- tmax_all[[grep(mon, names(tmax_all))]]
  app(layers, mean, na.rm = TRUE)
})
tmax_ensemble <- rast(tmax_ensemble)
names(tmax_ensemble) <- paste0("tmax_", month.abb)

# ------------------------------------------------------------------------------
# 4) Align grids
# ------------------------------------------------------------------------------
tmax_aligned <- resample(tmax_ensemble, tmin_ensemble, method = "bilinear")

# ------------------------------------------------------------------------------
# 5) Compute taverage
# ------------------------------------------------------------------------------
tavg_ensemble <- (tmin_ensemble + tmax_aligned) / 2
names(tavg_ensemble) <- paste0("tavg_ens_", month.abb)

# ------------------------------------------------------------------------------
# 6) Reproject study area
# ------------------------------------------------------------------------------
study_area_r <- vect(st_transform(study_area, crs(tavg_ensemble)))

# ------------------------------------------------------------------------------
# 7) Create a raster mask template from study area
# ------------------------------------------------------------------------------
#  — rasterize the study area polygon into a raster that *matches* your tavg grid
template_raster <- terra::rasterize(study_area_r, tavg_ensemble[[1]], field=1)

# ------------------------------------------------------------------------------
# 8) Apply mask using raster template
# ------------------------------------------------------------------------------
tavg_masked <- mask(tavg_ensemble, template_raster)

# ------------------------------------------------------------------------------
# 9) Save outputs (taverage)
# ------------------------------------------------------------------------------
out_dir <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/ens_taverage/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nlyr(tavg_masked))) {
  writeRaster(
    tavg_masked[[i]],
    filename = file.path(out_dir, paste0("tavg_ensemble_", month.abb[i], ".tif")),
    overwrite = TRUE
  )
}

# ------------------------------------------------------------------------------
# 10) Quick plot
# ------------------------------------------------------------------------------
plot(tavg_masked[[1]], main="tavg Ensemble (Jan)", col=viridis(10))
plot(study_area_r, add=TRUE, border="black")

cat("✔ taverage ensemble computed and saved to:\n", out_dir, "\n")
