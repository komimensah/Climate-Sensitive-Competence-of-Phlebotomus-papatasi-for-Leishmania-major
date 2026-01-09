###############################################################################
# FULL END-TO-END PIPELINE (ONE SCRIPT)
# Uncertainty-conserving raster workflow:
#  1) Fit + Monte Carlo propagate thermal uncertainty (RI_prob and VC_norm)
#  2) Spatialize MONTHLY rasters (mean / lwr / upr) from MONTHLY temperature stack
#  3) Aggregate MONTHLY -> QUARTERLY (and optional ANNUAL)
#  4) Composite risk maps with rodent background constraint:
#       - geometric, harmonic, multiplicative
#       - additive for ALL weights (10_90 ... 90_10)
#  5) Optional: mean composites across quarters
#
# NOTE: You must set the paths to your temperature folder, mask, and output dirs.
###############################################################################

rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(stats)
  library(terra)
})

set.seed(123)

###############################################################################
## A) PART 1 — THERMAL UNCERTAINTY (NON-SPATIAL) -> CURVES
###############################################################################

## ---- Utilities ----
logistic <- function(x) 1 / (1 + exp(-x))

RI_prob <- function(x) {
  out <- rep(NA_real_, length(x))
  out[is.finite(x) & x <= 1] <- 0
  idx <- which(is.finite(x) & x > 1)
  if (length(idx) > 0) out[idx] <- 1 - 1 / x[idx]
  out
}

q_safe <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

bootstrap_nls <- function(df, formula, start,
                          lower = NULL, upper = NULL,
                          nboot = 1000) {
  
  pnames <- names(start)
  P <- matrix(NA_real_, nboot, length(start))
  colnames(P) <- pnames
  
  for (i in seq_len(nboot)) {
    df_b <- df[sample.int(nrow(df), replace = TRUE), , drop = FALSE]
    
    fit_b <- try(
      nls(
        formula,
        data = df_b,
        start = start,
        algorithm = if (!is.null(lower)) "port" else "default",
        lower = lower,
        upper = upper,
        control = nls.control(warnOnly = TRUE, maxiter = 600)
      ),
      silent = TRUE
    )
    
    if (!inherits(fit_b, "try-error")) {
      P[i, ] <- coef(fit_b)
    }
  }
  
  P <- P[complete.cases(P), , drop = FALSE]
  P
}

minmax_norm_vec <- function(x, xmin, xmax, clamp = TRUE) {
  if (!is.finite(xmin) || !is.finite(xmax) || xmax == xmin) {
    return(rep(NA_real_, length(x)))
  }
  z <- (x - xmin) / (xmax - xmin)
  if (clamp) z <- pmin(1, pmax(0, z))
  z
}

## ---- Data (your provided points) ----
df_D1 <- data.frame(T = c(25.5, 26.5, 29.5, 23, 28),
                    y = c(0.164203612, 0.179533214, 0.183150183, 0.09596929, 0.1447178))

df_D2 <- data.frame(T = c(25.5, 26.5, 29.5),
                    y = c(0.039323634, 0.047147572, 0.06779661))

df_D3 <- data.frame(T = c(25.5, 26.5, 29.5),
                    y = c(0.057636888, 0.11778563, 0.154798762))

df_M1 <- data.frame(T = c(15,18,20,25,28,32), M = c(0.72,0.40,0.30,0.38,0.28,0.65))
df_M2 <- data.frame(T = c(15,18,20,25,28,32), M = c(0.95,0.85,0.50,0.70,0.30,0.25))
df_M3 <- data.frame(T = c(15,18,20,25,28,32), M = c(0.50,0.10,0.15,0.10,0.25,0.30))
df_M4 <- data.frame(T = c(15,18,20,25,28,32), M = c(0.0525,0.0859,0.1008,0.1295,0.0992,0.1736))

df_F <- data.frame(T = c(25.5, 26.5, 29.5),
                   F = c(42.49, 42.71, 41.30))

df_EIP <- data.frame(T = c(23, 26, 28),
                     EIP = c(7, 6.5, 5.0))

## ---- Bootstrap fits ----
nboot <- 500

boot_D1 <- bootstrap_nls(df_D1, y ~ exp(a + b*T + c*T^2),
                         start = list(a=-5,b=0.1,c=-0.002), nboot=nboot)

boot_D2 <- bootstrap_nls(df_D2, y ~ exp(a + b*T),
                         start = list(a=-3,b=0.1), nboot=nboot)

boot_D3 <- bootstrap_nls(df_D3, y ~ exp(a + b*T),
                         start = list(a=-3,b=0.1), nboot=nboot)

boot_M <- function(df) {
  bootstrap_nls(
    df,
    M ~ 1 / (1 + exp(-(a + b*T + c*T^2))),
    start = list(a=-3,b=0.2,c=-0.01),
    lower = c(-10,-5,-1),
    upper = c(10,5,1),
    nboot = nboot
  )
}

boot_M1 <- boot_M(df_M1)
boot_M2 <- boot_M(df_M2)
boot_M3 <- boot_M(df_M3)
boot_M4 <- boot_M(df_M4)

boot_F <- bootstrap_nls(
  df_F,
  F ~ A * exp(-((T - mu)^2)/(2*sigma^2)),
  start = list(A=45,mu=27,sigma=2),
  lower = c(0,20,0.1),
  upper = c(100,35,10),
  nboot = nboot
)

boot_EIP <- bootstrap_nls(
  df_EIP,
  EIP ~ a * exp(-((T - Topt)^2)/(2*sigma^2)),
  start = list(a=7,Topt=25,sigma=3),
  lower = c(0.1,15,0.2),
  upper = c(60,35,15),
  nboot = nboot
)

stopifnot(nrow(boot_D1)>20, nrow(boot_D2)>20, nrow(boot_D3)>20,
          nrow(boot_M1)>20, nrow(boot_M2)>20, nrow(boot_M3)>20, nrow(boot_M4)>20,
          nrow(boot_F)>20,  nrow(boot_EIP)>20)

## ---- Trait functions ----
D1  <- function(T,p) exp(p[1] + p[2]*T + p[3]*T^2)
D2  <- function(T,p) exp(p[1] + p[2]*T)
D3  <- function(T,p) exp(p[1] + p[2]*T)
M   <- function(T,p) logistic(p[1] + p[2]*T + p[3]*T^2)
Fec <- function(T,p) p[1]*exp(-((T-p[2])^2)/(2*p[3]^2))

EIP_fun <- function(T,p,Tmin=17.6,Tmax=35){
  e <- p[1]*exp(-((T-p[2])^2)/(2*p[3]^2))
  e[T < Tmin | T > Tmax] <- Inf
  e
}

a_fun <- function(T) 0.2

## ---- Monte Carlo propagation over bootstrapped parameters ----
T_seq <- seq(15, 35, by = 0.1)

Bprop <- 1000   # Monte Carlo draws for propagation
RIprob_draws <- matrix(NA_real_, nrow = length(T_seq), ncol = Bprop)
VC_draws     <- matrix(NA_real_, nrow = length(T_seq), ncol = Bprop)

for (i in seq_len(Bprop)) {
  
  pD1  <- boot_D1[sample.int(nrow(boot_D1), 1), ]
  pD2  <- boot_D2[sample.int(nrow(boot_D2), 1), ]
  pD3  <- boot_D3[sample.int(nrow(boot_D3), 1), ]
  pM1  <- boot_M1[sample.int(nrow(boot_M1), 1), ]
  pM2  <- boot_M2[sample.int(nrow(boot_M2), 1), ]
  pM3  <- boot_M3[sample.int(nrow(boot_M3), 1), ]
  pM4  <- boot_M4[sample.int(nrow(boot_M4), 1), ]
  pF   <- boot_F [sample.int(nrow(boot_F ), 1), ]
  pEIP <- boot_EIP[sample.int(nrow(boot_EIP), 1), ]
  
  RI_vals <- sapply(T_seq, function(T) {
    num <- Fec(T,pF) * D1(T,pD1) * D2(T,pD2) * D3(T,pD3)
    den <- (D1(T,pD1) + M(T,pM1)) *
      (D2(T,pD2) + M(T,pM2)) *
      (D3(T,pD3) + M(T,pM3)) *
      M(T,pM4)
    r <- num/den
    if (!is.finite(r)) NA_real_ else r
  })
  RIprob_draws[, i] <- RI_prob(RI_vals)
  
  VC_vals <- sapply(T_seq, function(T) {
    pT <- 1 - M(T,pM4)
    if (!is.finite(pT) || pT <= 0 || pT >= 1) return(0)
    eip <- EIP_fun(T, pEIP)
    if (!is.finite(eip)) return(0)
    vc <- (a_fun(T)^2 * pT^eip) / (-log(pT))
    if (!is.finite(vc) || vc < 0) 0 else vc
  })
  VC_draws[, i] <- VC_vals
}

## ---- Summaries (curves) ----
RI_curve <- data.frame(
  Temperature  = T_seq,
  RI_prob_mean = rowMeans(RIprob_draws, na.rm = TRUE),
  RI_prob_lwr  = apply(RIprob_draws, 1, q_safe, 0.025),
  RI_prob_upr  = apply(RIprob_draws, 1, q_safe, 0.975)
)
RI_curve$RI_prob_width <- RI_curve$RI_prob_upr - RI_curve$RI_prob_lwr

VC_mean <- rowMeans(VC_draws, na.rm = TRUE)
VC_lwr  <- apply(VC_draws, 1, q_safe, 0.025)
VC_upr  <- apply(VC_draws, 1, q_safe, 0.975)

###############################################################################
# VC NORMALIZATION — FIXED (RASTER-AWARE, BIOLOGICALLY CONSISTENT)
###############################################################################

# --- load the SAME temperature stack used for mapping ---
temp_folder <- '/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/ens_taverage_585/'
temp_files  <- sort(list.files(temp_folder, pattern="\\.tif$", full.names=TRUE))
stopifnot(length(temp_files) == 12)
T_stack <- rast(temp_files)

# --- temperature range actually present in rasters ---
t_min_r <- min(global(T_stack, "min", na.rm=TRUE)[,1], na.rm=TRUE)
t_max_r <- max(global(T_stack, "max", na.rm=TRUE)[,1], na.rm=TRUE)

# --- clip to biological model domain ---
Tmin_dom <- min(T_seq)
Tmax_dom <- max(T_seq)

t_min_eff <- max(t_min_r, Tmin_dom)
t_max_eff <- min(t_max_r, Tmax_dom)

# --- indices of temperature support actually used spatially ---
idx_eff <- which(T_seq >= t_min_eff & T_seq <= t_max_eff)

# --- robust VC normalization range (prevents tail domination) ---
VC_min <- as.numeric(quantile(VC_mean[idx_eff], probs = 0.01, na.rm = TRUE))
VC_max <- as.numeric(quantile(VC_mean[idx_eff], probs = 0.99, na.rm = TRUE))

# safety fallback
if (!is.finite(VC_min) || !is.finite(VC_max) || VC_max <= VC_min) {
  VC_min <- min(VC_mean[idx_eff], na.rm = TRUE)
  VC_max <- max(VC_mean[idx_eff], na.rm = TRUE)
}

# --- FINAL VC CURVE ---
VC_curve <- data.frame(
  Temperature  = T_seq,
  VC_norm_mean = minmax_norm_vec(VC_mean, VC_min, VC_max, clamp = TRUE),
  VC_norm_lwr  = minmax_norm_vec(VC_lwr,  VC_min, VC_max, clamp = TRUE),
  VC_norm_upr  = minmax_norm_vec(VC_upr,  VC_min, VC_max, clamp = TRUE)
)
VC_curve$VC_norm_width <- VC_curve$VC_norm_upr - VC_curve$VC_norm_lwr
###############################################################################
## B) PART 2 — BRIDGE: CURVES -> MONTHLY RASTERS (MEAN/LWR/UPR)
###############################################################################

## ---- Safe interpolator for terra::app ----
safe_interp <- function(x, y) {
  f <- approxfun(x, y, rule = 2)
  function(v) {
    out <- rep(NA_real_, length(v))
    ok  <- is.finite(v)
    if (any(ok)) out[ok] <- f(v[ok])
    out
  }
}

f_lam_mean <- safe_interp(RI_curve$Temperature, RI_curve$RI_prob_mean)
f_lam_lwr  <- safe_interp(RI_curve$Temperature, RI_curve$RI_prob_lwr)
f_lam_upr  <- safe_interp(RI_curve$Temperature, RI_curve$RI_prob_upr)

f_vc_mean  <- safe_interp(VC_curve$Temperature, VC_curve$VC_norm_mean)
f_vc_lwr   <- safe_interp(VC_curve$Temperature, VC_curve$VC_norm_lwr)
f_vc_upr   <- safe_interp(VC_curve$Temperature, VC_curve$VC_norm_upr)

## ---- INPUTS (EDIT PATHS) ----
temp_folder <- '/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/ens_taverage_585/'
temp_files  <- sort(list.files(temp_folder, pattern="\\.tif$", full.names=TRUE))
stopifnot(length(temp_files) == 12)

T_stack <- rast(temp_files)
names(T_stack) <- month.abb[1:12]
summary(T_stack)
# -----------------------------------------------------------------------------
# 4) Apply study area mask (SHAPEFILE or RASTER) — SAFE
# -----------------------------------------------------------------------------

# --- provide ONE of these ---
mask_path_vec <- '/Users/kagboka/Desktop/socioeeconomic_modelling_review/Study_area/study area dissolved_B.shp'  # shapefile
mask_path_rst <- NULL  # e.g. "/Users/kagboka/Desktop/socioeeconomic_modelling_review/mask.tif"

if (!is.null(mask_path_vec)) {
  
  # Vector mask (shapefile)
  aoi <- vect(mask_path_vec)
  
  # align CRS
  if (!same.crs(aoi, T_stack)) aoi <- project(aoi, crs(T_stack))
  
  # crop to AOI extent first (faster), then mask by polygon
  T_stack <- crop(T_stack, aoi)
  T_stack <- mask(T_stack, aoi)
  
} else if (!is.null(mask_path_rst)) {
  
  # Raster mask (0/1 or NA mask)
  mask_r <- rast(mask_path_rst)
  if (!same.crs(mask_r, T_stack)) mask_r <- project(mask_r, T_stack)
  mask_r <- resample(mask_r, T_stack[[1]], method = "near")
  
  # If your raster mask is 0/1:
  T_stack <- mask(T_stack, mask_r, maskvalues = 0)
  
  # If your raster mask is NA outside area (common), use instead:
  # T_stack <- mask(T_stack, mask_r)
  
} else {
  message("No mask applied (mask_path_vec and mask_path_rst are both NULL).")
}
## ---- OUTPUT DIRS ----
out_lambda_m <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/lambda_monthly_unc/"
out_vc_m     <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/vc_monthly_unc/"
dir.create(out_lambda_m, recursive=TRUE, showWarnings=FALSE)
dir.create(out_vc_m, recursive=TRUE, showWarnings=FALSE)

## ---- Build monthly λ and VC rasters (mean/lwr/upr) ----
lambda_mean_stack <- rast()
lambda_lwr_stack  <- rast()
lambda_upr_stack  <- rast()

vc_mean_stack <- rast()
vc_lwr_stack  <- rast()
vc_upr_stack  <- rast()

for (m in 1:nlyr(T_stack)) {
  
  mon <- month.abb[m]
  message("▶ Spatializing month: ", mon)
  
  Tm <- T_stack[[m]]
  
  lam_mean <- app(Tm, f_lam_mean); names(lam_mean) <- paste0("lambda_mean_", mon)
  lam_lwr  <- app(Tm, f_lam_lwr ); names(lam_lwr ) <- paste0("lambda_lwr_",  mon)
  lam_upr  <- app(Tm, f_lam_upr ); names(lam_upr ) <- paste0("lambda_upr_",  mon)
  
  vc_mean  <- app(Tm, f_vc_mean ); names(vc_mean ) <- paste0("vc_mean_", mon)
  vc_lwr   <- app(Tm, f_vc_lwr  ); names(vc_lwr  ) <- paste0("vc_lwr_",  mon)
  vc_upr   <- app(Tm, f_vc_upr  ); names(vc_upr  ) <- paste0("vc_upr_",  mon)
  
  writeRaster(lam_mean, file.path(out_lambda_m, paste0("lambda_mean_", mon, ".tif")), overwrite=TRUE)
  writeRaster(lam_lwr,  file.path(out_lambda_m, paste0("lambda_lwr_",  mon, ".tif")), overwrite=TRUE)
  writeRaster(lam_upr,  file.path(out_lambda_m, paste0("lambda_upr_",  mon, ".tif")), overwrite=TRUE)
  
  writeRaster(vc_mean,  file.path(out_vc_m, paste0("vc_mean_", mon, ".tif")), overwrite=TRUE)
  writeRaster(vc_lwr,   file.path(out_vc_m, paste0("vc_lwr_",  mon, ".tif")), overwrite=TRUE)
  writeRaster(vc_upr,   file.path(out_vc_m, paste0("vc_upr_",  mon, ".tif")), overwrite=TRUE)
  
  lambda_mean_stack <- c(lambda_mean_stack, lam_mean)
  lambda_lwr_stack  <- c(lambda_lwr_stack,  lam_lwr)
  lambda_upr_stack  <- c(lambda_upr_stack,  lam_upr)
  
  vc_mean_stack <- c(vc_mean_stack, vc_mean)
  vc_lwr_stack  <- c(vc_lwr_stack,  vc_lwr)
  vc_upr_stack  <- c(vc_upr_stack,  vc_upr)
}

###############################################################################
## C) PART 3 — AGGREGATE MONTHLY -> QUARTERLY (MEAN/LWR/UPR PATHS SEPARATE)
###############################################################################

out_lambda_q <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/RI_quarterly_unc/"
out_vc_q     <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/output/VC_quarterly_unc/"
dir.create(out_lambda_q, recursive=TRUE, showWarnings=FALSE)
dir.create(out_vc_q, recursive=TRUE, showWarnings=FALSE)

quarters <- list(
  Q1 = c("Jan","Feb","Mar"),
  Q2 = c("Apr","May","Jun"),
  Q3 = c("Jul","Aug","Sep"),
  Q4 = c("Oct","Nov","Dec")
)

get_layers <- function(stack, prefix, mons) {
  nms <- paste0(prefix, "_", mons)
  idx <- match(nms, names(stack))
  stack[[idx]]
}

for (q in names(quarters)) {
  
  mons <- quarters[[q]]
  message("▶ Aggregating quarter: ", q)
  
  lam_mean_q <- app(get_layers(lambda_mean_stack, "lambda_mean", mons), mean, na.rm=TRUE)
  lam_lwr_q  <- app(get_layers(lambda_lwr_stack,  "lambda_lwr",  mons), mean, na.rm=TRUE)
  lam_upr_q  <- app(get_layers(lambda_upr_stack,  "lambda_upr",  mons), mean, na.rm=TRUE)
  
  names(lam_mean_q) <- paste0("lambda_mean_", q)
  names(lam_lwr_q)  <- paste0("lambda_lwr_",  q)
  names(lam_upr_q)  <- paste0("lambda_upr_",  q)
  
  writeRaster(lam_mean_q, file.path(out_lambda_q, paste0("lambda_mean_", q, ".tif")), overwrite=TRUE)
  writeRaster(lam_lwr_q,  file.path(out_lambda_q, paste0("lambda_lwr_",  q, ".tif")), overwrite=TRUE)
  writeRaster(lam_upr_q,  file.path(out_lambda_q, paste0("lambda_upr_",  q, ".tif")), overwrite=TRUE)
  
  vc_mean_q <- app(get_layers(vc_mean_stack, "vc_mean", mons), mean, na.rm=TRUE)
  vc_lwr_q  <- app(get_layers(vc_lwr_stack,  "vc_lwr",  mons), mean, na.rm=TRUE)
  vc_upr_q  <- app(get_layers(vc_upr_stack,  "vc_upr",  mons), mean, na.rm=TRUE)
  
  names(vc_mean_q) <- paste0("vc_mean_", q)
  names(vc_lwr_q)  <- paste0("vc_lwr_",  q)
  names(vc_upr_q)  <- paste0("vc_upr_",  q)
  
  writeRaster(vc_mean_q, file.path(out_vc_q, paste0("vc_mean_", q, ".tif")), overwrite=TRUE)
  writeRaster(vc_lwr_q,  file.path(out_vc_q, paste0("vc_lwr_",  q, ".tif")), overwrite=TRUE)
  writeRaster(vc_upr_q,  file.path(out_vc_q, paste0("vc_upr_",  q, ".tif")), overwrite=TRUE)
}

###############################################################################
## D) PART 4 — COMPOSITE RISK (FULL PIPELINE, UNCERTAINTY PRESERVED)
###############################################################################

setwd("/Users/kagboka/Desktop/socioeeconomic_modelling_review/Composite risk_index/")

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

rodent_bg_path <- '/Users/kagboka/Desktop/socioeeconomic_modelling_review/future_predictions/ensemble_prediction_SSP585_2050.tif'
rodent_bg <- rast(rodent_bg_path)

for (q in names(quarters)) {
  
  message("▶ Composite risk ", q)
  
  lam_mean <- rast(file.path(out_lambda_q, paste0("lambda_mean_", q, ".tif")))
  lam_lwr  <- rast(file.path(out_lambda_q, paste0("lambda_lwr_",  q, ".tif")))
  lam_upr  <- rast(file.path(out_lambda_q, paste0("lambda_upr_",  q, ".tif")))
  
  vc_mean <- rast(file.path(out_vc_q, paste0("vc_mean_", q, ".tif")))
  vc_lwr  <- rast(file.path(out_vc_q, paste0("vc_lwr_",  q, ".tif")))
  vc_upr  <- rast(file.path(out_vc_q, paste0("vc_upr_",  q, ".tif")))
  
  if (!same.crs(vc_mean, lam_mean)) {
    vc_mean <- project(vc_mean, lam_mean)
    vc_lwr  <- project(vc_lwr,  lam_mean)
    vc_upr  <- project(vc_upr,  lam_mean)
  }
  
  vc_mean <- resample(vc_mean, lam_mean, method="bilinear")
  vc_lwr  <- resample(vc_lwr,  lam_mean, method="bilinear")
  vc_upr  <- resample(vc_upr,  lam_mean, method="bilinear")
  
  if (!same.crs(rodent_bg, lam_mean)) rodent_bg <- project(rodent_bg, lam_mean)
  rb <- resample(rodent_bg, lam_mean, method="bilinear")
  
  valid_mean <- is.finite(lam_mean) & is.finite(vc_mean) & is.finite(rb)
  valid_lwr  <- is.finite(lam_lwr ) & is.finite(vc_lwr ) & is.finite(rb)
  valid_upr  <- is.finite(lam_upr ) & is.finite(vc_upr ) & is.finite(rb)
  
  ## ---- Geometric ----
  gm_mean <- ifel(valid_mean, sqrt(lam_mean * vc_mean) * rb, NA)
  gm_lwr  <- ifel(valid_lwr,  sqrt(lam_lwr  * vc_lwr ) * rb, NA)
  gm_upr  <- ifel(valid_upr,  sqrt(lam_upr  * vc_upr ) * rb, NA)
  
  writeRaster(gm_mean, paste0("composite_geometric_mean_", q, ".tif"), overwrite=TRUE)
  writeRaster(gm_lwr,  paste0("composite_geometric_lwr_",  q, ".tif"), overwrite=TRUE)
  writeRaster(gm_upr,  paste0("composite_geometric_upr_",  q, ".tif"), overwrite=TRUE)
  
  ## ---- Harmonic ----
  hm_mean <- ifel(valid_mean & (lam_mean + vc_mean) > 0,
                  (2 * lam_mean * vc_mean / (lam_mean + vc_mean)) * rb, NA)
  hm_lwr  <- ifel(valid_lwr & (lam_lwr + vc_lwr) > 0,
                  (2 * lam_lwr * vc_lwr / (lam_lwr + vc_lwr)) * rb, NA)
  hm_upr  <- ifel(valid_upr & (lam_upr + vc_upr) > 0,
                  (2 * lam_upr * vc_upr / (lam_upr + vc_upr)) * rb, NA)
  
  writeRaster(hm_mean, paste0("composite_harmonic_mean_", q, ".tif"), overwrite=TRUE)
  writeRaster(hm_lwr,  paste0("composite_harmonic_lwr_",  q, ".tif"), overwrite=TRUE)
  writeRaster(hm_upr,  paste0("composite_harmonic_upr_",  q, ".tif"), overwrite=TRUE)
  
  ## ---- Multiplicative ----
  mul_mean <- ifel(valid_mean, (lam_mean * vc_mean) * rb, NA)
  mul_lwr  <- ifel(valid_lwr,  (lam_lwr  * vc_lwr ) * rb, NA)
  mul_upr  <- ifel(valid_upr,  (lam_upr  * vc_upr ) * rb, NA)
  
  writeRaster(mul_mean, paste0("composite_multiplicative_mean_", q, ".tif"), overwrite=TRUE)
  writeRaster(mul_lwr,  paste0("composite_multiplicative_lwr_",  q, ".tif"), overwrite=TRUE)
  writeRaster(mul_upr,  paste0("composite_multiplicative_upr_",  q, ".tif"), overwrite=TRUE)
  
  ## ---- Additive (ALL weights) ----
  for (w in names(weights)) {
    
    a <- weights[[w]][1]
    b <- weights[[w]][2]
    
    am_mean <- ifel(valid_mean, (a * lam_mean + b * vc_mean) * rb, NA)
    am_lwr  <- ifel(valid_lwr,  (a * lam_lwr  + b * vc_lwr ) * rb, NA)
    am_upr  <- ifel(valid_upr,  (a * lam_upr  + b * vc_upr ) * rb, NA)
    
    writeRaster(am_mean, paste0("composite_additive_", w, "_mean_", q, ".tif"), overwrite=TRUE)
    writeRaster(am_lwr,  paste0("composite_additive_", w, "_lwr_",  q, ".tif"), overwrite=TRUE)
    writeRaster(am_upr,  paste0("composite_additive_", w, "_upr_",  q, ".tif"), overwrite=TRUE)
  }
}

cat("\n✔ COMPLETE: monthly → quarterly → composite (mean/lwr/upr) saved.\n")

###############################################################################
## E) OPTIONAL — MEAN COMPOSITES ACROSS QUARTERS (MEAN/LWR/UPR SEPARATELY)
###############################################################################

mean_dir <- file.path(getwd(), "composite_mean_unc")
dir.create(mean_dir, showWarnings = FALSE, recursive = TRUE)

all_q_files <- list.files(getwd(), pattern="^composite_.*_(mean|lwr|upr)_Q[1-4]\\.tif$",
                          full.names=TRUE)

get_category <- function(f) {
  nm <- basename(f)
  nm <- sub("_Q[1-4]\\.tif$", "", nm)
  nm
}

cats <- unique(vapply(all_q_files, get_category, character(1)))

for (catg in cats) {
  files <- all_q_files[vapply(all_q_files, function(f) get_category(f) == catg, logical(1))]
  if (length(files) != 4) next
  r <- rast(files)
  r_mean <- app(r, mean, na.rm=TRUE)
  writeRaster(r_mean, file.path(mean_dir, paste0(catg, "_MEAN.tif")), overwrite=TRUE)
}

cat("\n✔ MEAN across quarters written to composite_mean_unc/\n")

