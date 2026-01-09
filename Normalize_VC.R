
library(terra)
library(viridis)

# =============================================================================
# 0) Vectorial Capacity functions (UNCHANGED – mechanistic)
# =============================================================================
EIP <- function(T,
                a     = 7.8024,
                Topt  = 24.7428,
                sigma = 3.0263,
                Tmin  = 17.6,
                Tmax  = 35) {
  
  eip <- a * exp(-((T - Topt)^2) / (2 * sigma^2))
  eip[T < Tmin | T > Tmax] <- Inf
  eip
}

logistic <- function(x) 1 / (1 + exp(-x))

M4 <- function(T) logistic(-3.4409589 + 0.0525869 * T + 0.0001129 * T^2)
p  <- function(T) 1 - M4(T)
a  <- function(T) 0.2

VC_T <- function(T) {
  pT   <- p(T)
  eipT <- EIP(T)
  VC   <- (a(T)^2 * pT^eipT) / (-log(pT))
  VC[!is.finite(VC) | VC < 0] <- 0
  VC
}

# =============================================================================
# 1) Load monthly temperature rasters (12 layers)
# =============================================================================
temp_dir <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/wc2.1_2.5m_tavg/"
shp_path <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/Study_area/study area dissolved_B.shp"

out_raw_dir  <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/VC_raw/"
out_norm_dir <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/VC_norm/"
out_q_dir    <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/VC_quarterly/"
out_4m_dir   <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/VC_4month/"
out_ann_dir  <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/VC_annual/"

dir.create(out_raw_dir,  recursive=TRUE, showWarnings=FALSE)
dir.create(out_norm_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(out_q_dir,    recursive=TRUE, showWarnings=FALSE)
dir.create(out_4m_dir,   recursive=TRUE, showWarnings=FALSE)
dir.create(out_ann_dir,  recursive=TRUE, showWarnings=FALSE)

temp_files <- sort(list.files(temp_dir, "\\.tif$", full.names = TRUE))
stopifnot(length(temp_files) == 12)

T_stack <- rast(temp_files)
names(T_stack) <- month.abb

shp <- vect(shp_path) |> project(crs(T_stack))
T_stack <- mask(crop(T_stack, shp), shp)

# =============================================================================
# 2) Compute monthly VC (RAW, mechanistic)
# =============================================================================
VC_stack <- app(T_stack, VC_T)
names(VC_stack) <- paste0("VC_", month.abb)

# Save raw VC
for (i in 1:nlyr(VC_stack)) {
  writeRaster(
    VC_stack[[i]],
    file.path(out_raw_dir, paste0("VC_", month.abb[i], ".tif")),
    overwrite = TRUE
  )
}

# =============================================================================
# 3) Min–max normalization (CORRECT PLACE)
# =============================================================================
minmax_norm <- function(r) {
  rmin <- global(r, "min", na.rm = TRUE)[1,1]
  rmax <- global(r, "max", na.rm = TRUE)[1,1]
  if (!is.finite(rmin) || !is.finite(rmax) || rmax == rmin) {
    return(r * NA)
  }
  (r - rmin) / (rmax - rmin)
}

VC_norm_stack <- rast(VC_stack)

for (i in 1:nlyr(VC_stack)) {
  VC_norm_stack[[i]] <- minmax_norm(VC_stack[[i]])
  writeRaster(
    VC_norm_stack[[i]],
    file.path(out_norm_dir, paste0("VCnorm_", month.abb[i], ".tif")),
    overwrite = TRUE
  )
}

names(VC_norm_stack) <- paste0("VCnorm_", month.abb)

# =============================================================================
# 4) Aggregate NORMALIZED VC (THIS WAS YOUR BUG)
# =============================================================================

# ---- Quarterly (3-month) ----
vc_q <- list(
  Q1 = c("VCnorm_Jan","VCnorm_Feb","VCnorm_Mar"),
  Q2 = c("VCnorm_Apr","VCnorm_May","VCnorm_Jun"),
  Q3 = c("VCnorm_Jul","VCnorm_Aug","VCnorm_Sep"),
  Q4 = c("VCnorm_Oct","VCnorm_Nov","VCnorm_Dec")
)

for (q in names(vc_q)) {
  r <- app(VC_norm_stack[[vc_q[[q]]]], mean, na.rm=TRUE)
  writeRaster(r, file.path(out_q_dir, paste0("VCnorm_", q, ".tif")), overwrite=TRUE)
}

# ---- 4-month seasons ----
vc_4m <- list(
  JFMA = c("VCnorm_Jan","VCnorm_Feb","VCnorm_Mar","VCnorm_Apr"),
  MJJA = c("VCnorm_May","VCnorm_Jun","VCnorm_Jul","VCnorm_Aug"),
  SOND = c("VCnorm_Sep","VCnorm_Oct","VCnorm_Nov","VCnorm_Dec")
)

for (s in names(vc_4m)) {
  r <- app(VC_norm_stack[[vc_4m[[s]]]], mean, na.rm=TRUE)
  writeRaster(r, file.path(out_4m_dir, paste0("VCnorm_", s, ".tif")), overwrite=TRUE)
}

# ---- Annual ----
VCnorm_annual <- app(VC_norm_stack, mean, na.rm=TRUE)
writeRaster(VCnorm_annual,
            file.path(out_ann_dir, "VCnorm_Annual.tif"),
            overwrite=TRUE)

# =============================================================================
# 5) SANITY CHECK (THIS SHOULD NOW BE TRUE)
# =============================================================================
print(global(VCnorm_annual, c("min","max"), na.rm=TRUE))
plot(VCnorm_annual, col=viridis(10), main="VC Normalized – Annual Mean")
