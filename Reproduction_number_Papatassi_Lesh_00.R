# =============================================================================
# =============================================================================
# RI raster workflow with thresholded probability transform
# =============================================================================

library(terra)
library(sf)
library(viridis)

# -----------------------------------------------------------------------------
# 1) Define Risk Index functions
# -----------------------------------------------------------------------------

logistic <- function(x) 1 / (1 + exp(-x))

D1 <- function(T) exp(-17.30323 + 1.10520 * T - 0.01962 * T^2)
D2 <- function(T) exp(-6.539502 + 0.130549 * T)
D3 <- function(T) exp(-6.97742 + 0.17405 * T)

M1 <- function(T) logistic(10 - 0.917877 * T + 0.019355 * T^2)
M2 <- function(T) logistic(5.4124 - 0.258013 * T + 0.001726 * T^2)
M3 <- function(T) logistic(10 - 0.996416 * T + 0.020652 * T^2)
M4 <- function(T) logistic(-3.4409589 + 0.0525869 * T + 0.0001129 * T^2)

Fec <- function(T) 54.8478 * exp(-((T - 27.6971)^2) / (2 * 2.3114^2))

RI_fun <- function(T) {
  num <- Fec(T) * D1(T) * D2(T) * D3(T)
  den <- (D1(T) + M1(T)) * (D2(T) + M2(T)) * (D3(T) + M3(T)) * M4(T)
  ri  <- num / den
  ri[!is.finite(ri)] <- NA
  ri
}

# -----------------------------------------------------------------------------
# 2) Define probability transform function
# -----------------------------------------------------------------------------

RI_prob <- function(x) {
  # x is a numeric vector from terra::app
  out <- x * NA_real_     # initialize with NA
  out[x <= 1] <- 0        # RI <= 1 → 0
  idx <- which(x > 1)     # RI > 1
  if (length(idx) > 0) {
    out[idx] <- 1 - 1 / x[idx]
  }
  out
}

# -----------------------------------------------------------------------------
# 3) Load monthly temperature rasters
# -----------------------------------------------------------------------------

temp_folder <- '/Users/kagboka/Desktop/socioeeconomic_modelling_review/output_present/ens_taverage_585/'
temp_files  <- sort(list.files(temp_folder, pattern="\\.tif$", full.names=TRUE))
stopifnot(length(temp_files) == 12)

T_stack <- rast(temp_files)
names(T_stack) <- month.abb[1:12]

# -----------------------------------------------------------------------------
# 4) Apply study area mask
# -----------------------------------------------------------------------------

shp <- vect("/Users/kagboka/Desktop/socioeeconomic_modelling_review/Study_area/study area dissolved_B.shp")
shp <- project(shp, crs(T_stack))

T_stack <- crop(T_stack, shp)
T_stack <- mask(T_stack, shp)

print(summary(T_stack))

# -----------------------------------------------------------------------------
# 5) Compute raw RI stack
# -----------------------------------------------------------------------------

RI_stack <- app(T_stack, fun=function(vals) RI_fun(vals))
names(RI_stack) <- paste0("RI_", month.abb[1:12])

# Optional: visualize
for (i in 1:nlyr(RI_stack)) {
  plot(RI_stack[[i]],
       main=paste("Raw RI -", month.abb[i]),
       col = viridis(10))
}

# -----------------------------------------------------------------------------
# 6) Apply probability transform and save
# -----------------------------------------------------------------------------

out_RI_raw_dir  <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/RI_raw/"
out_RI_prob_dir <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/RI_prob/"

dir.create(out_RI_raw_dir,  recursive=TRUE, showWarnings=FALSE)
dir.create(out_RI_prob_dir, recursive=TRUE, showWarnings=FALSE)

RI_prob_stack <- rast(RI_stack)  # empty stack to fill

for (i in seq_len(nlyr(RI_stack))) {
  raw_r <- RI_stack[[i]]
  
  # Apply transform
  prob_r <- terra::app(raw_r, RI_prob)
  
  # Name
  names(prob_r) <- paste0("RIprob_", month.abb[i])
  RI_prob_stack[[i]] <- prob_r
  
  # Save
  writeRaster(raw_r,
              filename = file.path(out_RI_raw_dir, paste0("RI_", month.abb[i], ".tif")),
              overwrite = TRUE)
  writeRaster(prob_r,
              filename = file.path(out_RI_prob_dir, paste0("RIprob_", month.abb[i], ".tif")),
              overwrite = TRUE)
}

# -----------------------------------------------------------------------------
# 7) Plot example
# -----------------------------------------------------------------------------

plot(
  RI_prob_stack[[1]],
  main = "RI Probability (Jan)",
  col  = viridis::viridis(10)
)
plot(shp, add=TRUE, border="black", lwd=1.2)

cat("\n✔ RI probability rasters complete.\n",
    "✔ Raw in:", out_RI_raw_dir, "\n",
    "✔ Probability in:", out_RI_prob_dir, "\n")
# =============================================================================
# 8) Compute aggregated RI probabilities
# =============================================================================

# Make sure RI_prob_stack has named layers like RIprob_Jan … RIprob_Dec
print(names(RI_prob_stack))

# Define seasons
seasons_3m <- list(
  Q1 = c("RIprob_Jan", "RIprob_Feb", "RIprob_Mar"),
  Q2 = c("RIprob_Apr", "RIprob_May", "RIprob_Jun"),
  Q3 = c("RIprob_Jul", "RIprob_Aug", "RIprob_Sep"),
  Q4 = c("RIprob_Oct", "RIprob_Nov", "RIprob_Dec")
)

seasons_4m <- list(
  JFMAMJ = c("RIprob_Jan", "RIprob_Feb", "RIprob_Mar", "RIprob_Apr"),
  MJJASO = c("RIprob_May", "RIprob_Jun", "RIprob_Jul", "RIprob_Aug"),
  SONDJF = c("RIprob_Sep", "RIprob_Oct", "RIprob_Nov", "RIprob_Dec")
)

# Annual all 12
annual_names <- names(RI_prob_stack)

# Output folders
out_Q3_dir  <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/RI_quarterly/"
out_4m_dir  <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/RI_4month/"
out_ann_dir <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/RI_annual/"

dir.create(out_Q3_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(out_4m_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(out_ann_dir, recursive=TRUE, showWarnings=FALSE)

# -----------------------------------------------------------------------------
# 8a) Quarterly (3-month) means
# -----------------------------------------------------------------------------
for (q in names(seasons_3m)) {
  layers <- seasons_3m[[q]]
  q_stack <- RI_prob_stack[[layers]]
  q_mean  <- app(q_stack, mean, na.rm=TRUE)
  names(q_mean) <- paste0("RIprob_", q)
  
  # Save
  writeRaster(
    q_mean,
    filename = file.path(out_Q3_dir, paste0("RIprob_", q, ".tif")),
    overwrite = TRUE
  )
  
  # Plot (optional)
  plot(q_mean, main=paste("RI Probability (Quarter:", q, ")"), col=viridis(10))
}

# -----------------------------------------------------------------------------
# 8b) 4-month seasonal means
# -----------------------------------------------------------------------------
for (s in names(seasons_4m)) {
  layers <- seasons_4m[[s]]
  s_stack <- RI_prob_stack[[layers]]
  s_mean  <- app(s_stack, mean, na.rm=TRUE)
  names(s_mean) <- paste0("RIprob_", s)
  
  # Save
  writeRaster(
    s_mean,
    filename = file.path(out_4m_dir, paste0("RIprob_", s, ".tif")),
    overwrite = TRUE
  )
  
  # Plot (optional)
  plot(s_mean, main=paste("RI Probability (4-mo:", s, ")"), col=viridis(10))
}

# -----------------------------------------------------------------------------
# 8c) Annual mean
# -----------------------------------------------------------------------------
annual_mean <- app(RI_prob_stack, mean, na.rm=TRUE)
names(annual_mean) <- "RIprob_Annual"

writeRaster(
  annual_mean,
  filename = file.path(out_ann_dir, "RIprob_Annual.tif"),
  overwrite = TRUE
)

plot(annual_mean, main="RI Probability (Annual Mean)", col=viridis(10))

cat("\n✔ Quarterly, 4-mo, and Annual RI probabilities complete.\n",
    "✔ Quarterlies in:", out_Q3_dir, "\n",
    "✔ 4-month means in:", out_4m_dir, "\n",
    "✔ Annual mean in:", out_ann_dir, "\n")

