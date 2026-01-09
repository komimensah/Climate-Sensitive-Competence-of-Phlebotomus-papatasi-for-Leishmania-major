# ==============================================================================
# Proper background generation for SDM evaluation (North Africa)
# - avoids overlap by ensuring environmental coverage
# - supports bias file (population) and distance buffer
# ==============================================================================

rm(list=ls()); gc()
library(terra)
library(sf)
library(dplyr)

# ------------------------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------------------------
target_crs <- "EPSG:3395"  # meters

# Inputs
pop_path   <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/Population_Maghreb.tif"
obs_path   <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/Final data DIsease.csv"
obs_lon    <- "lon"
obs_lat    <- "lat"

# Environmental stack used for stratification (IMPORTANT!)
# Use your main predictor stack or a compact set (e.g., tavg, aridity, ndvi proxy, soil texture)
env_stack_path <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/Temperature_mean.tif"

# Study boundary polygon (recommended). If you don't have one, use pop extent.
boundary_path <- NA  # e.g. "/path/to/north_africa_boundary.shp"

# Background settings
n_background <- 1000
buffer_km    <- 100

# Choose background method:
# "env_strat"  = stratified by environmental space (best baseline)
# "env_bias"   = env_strat + bias weights from population (recommended if surveillance bias is strong)
# "target_grp" = use provided target-group points (requires target_points_path)
bg_method <- "env_bias"

# Only needed for target-group background
target_points_path <- NA  # e.g. "/path/to/other_surveillance_points.shp" or CSV

set.seed(123)

km_to_m <- function(km) km * 1000

# ------------------------------------------------------------------------------
# 1) LOAD DATA
# ------------------------------------------------------------------------------
pop_r <- rast(pop_path) |> project(target_crs)

env <- rast(env_stack_path)
env <- project(env, target_crs)
# align env to population grid (or vice versa) to avoid extract problems
env <- resample(env, pop_r, method = "bilinear")

obs_df <- read.csv(obs_path)
obs_sf <- st_as_sf(obs_df, coords = c(obs_lon, obs_lat), crs = 4326)
obs_sf <- st_transform(obs_sf, target_crs)
pts_presence <- vect(obs_sf)

# Boundary
if (!is.na(boundary_path) && file.exists(boundary_path)) {
  boundary <- vect(boundary_path) |> project(target_crs)
} else {
  boundary <- as.polygons(ext(pop_r), crs = crs(pop_r))
}

# Allowed area = boundary minus buffer around presences
pres_buffer <- terra::buffer(
  x     = pts_presence,
  width = km_to_m(buffer_km)
)
allowed_area <- erase(boundary, pres_buffer)

# Mask env and pop to allowed area
env_allowed <- mask(env, allowed_area)
pop_allowed <- mask(log1p(pop_r), allowed_area)
pop_allowed[pop_allowed <= 0] <- NA

# ------------------------------------------------------------------------------
# 2) HELPER: build environmental strata
# ------------------------------------------------------------------------------
# We stratify using the first 2–3 PCA-like axes would be ideal, but we keep it simple:
# choose 2–3 representative layers (first 3 layers by default) and bin them.
pick_layers <- function(env_rast, k = 3) {
  k <- min(k, nlyr(env_rast))
  env_rast[[1:k]]
}

make_env_strata <- function(env_rast, n_bins = 6) {
  E <- pick_layers(env_rast, k = 3)
  
  # create bins per layer using quantiles
  strata <- NULL
  for (j in 1:nlyr(E)) {
    v <- values(E[[j]], mat = FALSE)
    qs <- quantile(v, probs = seq(0,1,length.out = n_bins+1), na.rm = TRUE)
    qs <- unique(qs)
    if (length(qs) < 3) {
      # fallback: fewer bins
      qs <- quantile(v, probs = c(0, .5, 1), na.rm = TRUE)
      qs <- unique(qs)
    }
    strata_j <- classify(E[[j]], rcl = cbind(qs[-length(qs)], qs[-1], seq_len(length(qs)-1)), include.lowest = TRUE)
    if (is.null(strata)) strata <- strata_j else strata <- strata * 100 + strata_j
  }
  strata
}

# ------------------------------------------------------------------------------
# 3) METHOD A: Environmental stratified background (recommended baseline)
# ------------------------------------------------------------------------------
bg_env_strat <- function(env_rast, allowed_poly, n_bg, n_bins = 6, oversample = 15) {
  
  strata <- make_env_strata(env_rast, n_bins = n_bins)
  strata <- mask(strata, allowed_poly)
  
  # Oversample candidates uniformly within allowed area
  cand <- spatSample(allowed_poly, size = n_bg * oversample, method = "random")
  
  # Extract stratum id for each candidate
  sid <- terra::extract(strata, cand)[, 2]
  ok  <- which(is.finite(sid))
  cand <- cand[ok]
  sid  <- sid[ok]
  
  # Allocate equal samples per stratum (as much as possible)
  tab <- table(sid)
  strata_ids <- as.numeric(names(tab))
  per <- ceiling(n_bg / length(strata_ids))
  
  pick <- c()
  for (s in strata_ids) {
    idx <- which(sid == s)
    if (length(idx) == 0) next
    take <- min(per, length(idx))
    pick <- c(pick, sample(idx, take))
    if (length(pick) >= n_bg) break
  }
  
  # If still short, fill randomly
  if (length(pick) < n_bg) {
    remaining <- setdiff(seq_along(sid), pick)
    add <- sample(remaining, size = n_bg - length(pick))
    pick <- c(pick, add)
  }
  
  cand[pick[1:n_bg]]
}

# ------------------------------------------------------------------------------
# 4) METHOD B: Environmental stratified + bias raster (population) weights
#    (best compromise: keeps env coverage but respects sampling bias)
# ------------------------------------------------------------------------------
bg_env_bias <- function(env_rast, bias_rast, allowed_poly, n_bg, n_bins = 6, oversample = 30) {
  
  strata <- make_env_strata(env_rast, n_bins = n_bins)
  strata <- mask(strata, allowed_poly)
  
  # Oversample candidates uniformly within allowed area
  cand <- spatSample(allowed_poly, size = n_bg * oversample, method = "random")
  
  sid <- terra::extract(strata, cand)[, 2]
  bw  <- terra::extract(bias_rast, cand)[, 2]
  
  ok <- which(is.finite(sid) & is.finite(bw) & bw > 0)
  cand <- cand[ok]
  sid  <- sid[ok]
  bw   <- bw[ok]
  
  tab <- table(sid)
  strata_ids <- as.numeric(names(tab))
  per <- ceiling(n_bg / length(strata_ids))
  
  pick <- c()
  for (s in strata_ids) {
    idx <- which(sid == s)
    if (length(idx) == 0) next
    take <- min(per, length(idx))
    # weighted sample within stratum
    pick <- c(pick, sample(idx, size = take, replace = TRUE, prob = bw[idx]))
    if (length(pick) >= n_bg) break
  }
  
  # Fill if short
  if (length(pick) < n_bg) {
    remaining <- setdiff(seq_along(sid), pick)
    add <- sample(remaining, size = n_bg - length(pick), replace = TRUE, prob = bw[remaining])
    pick <- c(pick, add)
  }
  
  cand[pick[1:n_bg]]
}

# ------------------------------------------------------------------------------
# 5) METHOD C: Target-group background (if you have surveillance points)
# ------------------------------------------------------------------------------
bg_target_group <- function(target_points_vect, allowed_poly, n_bg) {
  # keep only points inside allowed area
  tg <- crop(target_points_vect, allowed_poly)
  tg <- mask(tg, allowed_poly)
  
  if (nrow(tg) < n_bg) {
    stop("Not enough target-group points in allowed area for requested n_background.")
  }
  tg[sample(seq_len(nrow(tg)), n_bg)]
}

# ------------------------------------------------------------------------------
# 6) GENERATE BACKGROUND (choose method)
# ------------------------------------------------------------------------------
if (bg_method == "env_strat") {
  bg_points <- bg_env_strat(env_allowed, allowed_area, n_background, n_bins = 6, oversample = 15)
  
} else if (bg_method == "env_bias") {
  bg_points <- bg_env_bias(env_allowed, pop_allowed, allowed_area, n_background, n_bins = 6, oversample = 30)
  
} else if (bg_method == "target_grp") {
  
  if (is.na(target_points_path)) stop("Provide target_points_path for target-group background.")
  if (grepl("\\.csv$", tolower(target_points_path))) {
    tg_df <- read.csv(target_points_path)
    tg_sf <- st_as_sf(tg_df, coords = c(obs_lon, obs_lat), crs = 4326) |> st_transform(target_crs)
    tg    <- vect(tg_sf)
  } else {
    tg <- vect(target_points_path) |> project(target_crs)
  }
  
  bg_points <- bg_target_group(tg, allowed_area, n_background)
  
} else {
  stop("Unknown bg_method. Use 'env_strat', 'env_bias', or 'target_grp'.")
}

# ------------------------------------------------------------------------------
# 7) QUICK DIAGNOSTICS (HIGHLY RECOMMENDED)
# ------------------------------------------------------------------------------
cat("\nBackground generated with method:", bg_method, "\n")
cat("Presences:", nrow(pts_presence), " | Background:", nrow(bg_points), "\n")

plot(pop_allowed, main = "Bias raster (log1p population) in allowed area")
plot(bg_points, add = TRUE, pch = 16, cex = 0.3)
plot(pts_presence, add = TRUE, col = "red", pch = 16)



# ==============================================================================
# NEXT STEPS (AFTER YOUR BACKGROUND CODE)
# - Build presence+background table
# - Spatial block CV
# - Evaluate all composite rasters
# - Metrics: AUC-PR, TSS (thr that maximizes TSS), partial ROC
# - Save fold metrics + summary ranking
# ==============================================================================

# ---- Packages needed for evaluation ----
suppressPackageStartupMessages({
  library(blockCV)
  library(stringr)
  library(pROC)
  library(PRROC)
  library(ecospat)
  library(ggplot2)
  library(zoo)
})

# ---- Paths / settings ----
composite_dir   <- '/Users/kagboka/Desktop/socioeeconomic_modelling_review/Composite risk_index_PRESENT copy/'
out_dir         <- file.path(composite_dir, "evaluation_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_fold_csv    <- file.path(out_dir, "fold_metrics.csv")
out_summary_csv <- file.path(out_dir, "summary_metrics.csv")
out_calib_dir   <- file.path(out_dir, "calibration_plots")
dir.create(out_calib_dir, showWarnings = FALSE, recursive = TRUE)

# Spatial block settings
block_km   <- 100
k_folds    <- 5
block_iter <- 100

# Partial ROC settings
fpr_max <- 0.10

# Rare disease safe thresholds
min_test_points <- 10          # minimum points per fold
min_pos_fold    <- 3           # minimum presences per fold
min_bg_fold     <- 10          # minimum background per fold

km_to_m <- function(km) km * 1000

# ------------------------------------------------------------------------------
# 8) Build Presence + Background table (NO geom() nonsense, robust)
# ------------------------------------------------------------------------------
pres_xy <- terra::crds(pts_presence, df = TRUE)
bg_xy   <- terra::crds(bg_points,   df = TRUE)

pres_df <- data.frame(x = pres_xy[,1], y = pres_xy[,2], resp = 1)
bg_df   <- data.frame(x = bg_xy[,1],   y = bg_xy[,2],   resp = 0)

all_df <- rbind(pres_df, bg_df)

# sf for blockCV, and also terra vect for extraction
all_sf <- st_as_sf(all_df, coords = c("x","y"), crs = target_crs, remove = FALSE)

cat("\n---- Data check ----\n")
cat("Total points:", nrow(all_df),
    "| Presences:", sum(all_df$resp == 1),
    "| Background:", sum(all_df$resp == 0), "\n")

stopifnot(nrow(all_df) == nrow(all_sf))
stopifnot(all(all_df$resp %in% c(0,1)))

# ------------------------------------------------------------------------------
# 9) Spatial block cross-validation folds
# ------------------------------------------------------------------------------
sb <- spatialBlock(
  speciesData     = all_sf,
  theRange        = km_to_m(block_km),
  k               = k_folds,
  selection       = "random",
  iteration       = block_iter,
  biomod2Format   = FALSE,
  progress        = FALSE
)

fold_id <- sb$foldID

# MUST be same length as all_df
if (length(fold_id) != nrow(all_df)) {
  stop("fold_id length != number of points. Something is wrong in spatialBlock output.")
}

cat("\n---- Fold check (counts) ----\n")
print(table(fold_id))
print(table(fold_id, all_df$resp))

# ------------------------------------------------------------------------------
# 10) Metric helpers (robust + NA-safe)
# ------------------------------------------------------------------------------
auc_pr <- function(scores, labels) {
  scores <- as.numeric(scores); labels <- as.integer(labels)
  ok <- is.finite(scores) & !is.na(labels)
  scores <- scores[ok]; labels <- labels[ok]
  if (length(unique(labels)) < 2) return(NA_real_)
  if (sum(labels == 1) < 2 || sum(labels == 0) < 5) return(NA_real_)
  pr <- try(PRROC::pr.curve(
    scores.class0 = scores[labels == 1],
    scores.class1 = scores[labels == 0],
    curve = FALSE
  ), silent = TRUE)
  if (inherits(pr, "try-error")) return(NA_real_)
  as.numeric(pr$auc.integral)
}

tss_best <- function(scores, labels) {
  scores <- as.numeric(scores); labels <- as.integer(labels)
  ok <- is.finite(scores) & !is.na(labels)
  scores <- scores[ok]; labels <- labels[ok]
  if (length(unique(labels)) < 2) return(list(TSS = NA_real_, thr = NA_real_))
  roc_obj <- try(pROC::roc(labels, scores, quiet = TRUE, direction = "<"), silent = TRUE)
  if (inherits(roc_obj, "try-error")) return(list(TSS = NA_real_, thr = NA_real_))
  
  coords_all <- pROC::coords(
    roc_obj, x = "all",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  tss_vals <- coords_all$sensitivity + coords_all$specificity - 1
  j <- which.max(tss_vals)
  list(TSS = as.numeric(tss_vals[j]), thr = as.numeric(coords_all$threshold[j]))
}

boyce_index <- function(pred_pres, pred_bg, nclass = 0) {
  pred_pres <- pred_pres[is.finite(pred_pres)]
  pred_bg   <- pred_bg[is.finite(pred_bg)]
  if (length(pred_pres) < 10 || length(pred_bg) < 50) return(NA_real_)
  b <- try(ecospat.boyce(fit = pred_bg, obs = pred_pres, nclass = nclass, PEplot = FALSE), silent = TRUE)
  if (inherits(b, "try-error")) return(NA_real_)
  as.numeric(b$Spearman.cor)
}

partial_roc_ratio <- function(scores, labels, fpr_max = 0.1) {
  scores <- as.numeric(scores); labels <- as.integer(labels)
  ok <- is.finite(scores) & !is.na(labels)
  scores <- scores[ok]; labels <- labels[ok]
  if (length(unique(labels)) < 2) return(NA_real_)
  
  roc_obj <- try(pROC::roc(labels, scores, quiet = TRUE, direction = "<"), silent = TRUE)
  if (inherits(roc_obj, "try-error")) return(NA_real_)
  
  roc_df <- data.frame(
    tpr = rev(roc_obj$sensitivities),
    fpr = rev(1 - roc_obj$specificities)
  )
  roc_df <- roc_df[is.finite(roc_df$fpr) & is.finite(roc_df$tpr), ]
  roc_df <- roc_df[roc_df$fpr <= fpr_max, ]
  if (nrow(roc_df) < 2) return(NA_real_)
  
  ord <- order(roc_df$fpr)
  x <- roc_df$fpr[ord]
  y <- roc_df$tpr[ord]
  pAUC <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  
  # random expectation in [0, fpr_max] triangle
  pAUC_rand <- (fpr_max^2) / 2
  pAUC / pAUC_rand
}

calibration_summary <- function(scores, labels, n_bins = 10) {
  scores <- as.numeric(scores); labels <- as.integer(labels)
  ok <- is.finite(scores) & !is.na(labels)
  df <- data.frame(score = scores[ok], y = labels[ok])
  if (nrow(df) < 30) return(NULL)
  
  q <- quantile(df$score, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  q <- unique(q)
  if (length(q) < 4) return(NULL)
  
  df$bin <- cut(df$score, breaks = q, include.lowest = TRUE)
  df %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(score, na.rm = TRUE),
      obs_rate  = mean(y, na.rm = TRUE),
      n         = n(),
      .groups   = "drop"
    )
}

plot_calibration <- function(cal_df, title = "Calibration") {
  ggplot(cal_df, aes(x = mean_pred, y = obs_rate)) +
    geom_point() +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(x = "Mean predicted score", y = "Observed presence fraction", title = title) +
    theme_minimal()
}

# ------------------------------------------------------------------------------
# 11) List composite rasters to evaluate
# ------------------------------------------------------------------------------
comp_files <- sort(list.files(composite_dir, pattern = "^composite_.*\\.tif$", full.names = TRUE))
parse_name <- function(f) {
  
  nm <- basename(f)
  
  # extract quarter
  quarter <- stringr::str_extract(nm, "Q[1-4]")
  
  # remove prefix and suffix
  method <- sub("^composite_", "", nm)
  method <- sub("_Q[1-4]\\.tif$", "", method)
  
  # 🔥 REMOVE TRAILING REPLICATE IDS (_1, _2, _3, ...)
  method <- sub("_[0-9]+$", "", method)
  
  data.frame(
    file    = f,
    method  = method,
    quarter = quarter,
    stringsAsFactors = FALSE
  )
}
comp_info <- dplyr::bind_rows(lapply(comp_files, parse_name))
# ------------------------------------------------------------------------------
# 12) Evaluation loop (FIXED & SAFE)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Helper: force metrics to length-1 numeric (never length 0)
# ------------------------------------------------------------------------------
safe_scalar <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  as.numeric(x[1])
}
results <- list()

for (i in seq_len(nrow(comp_info))) {
  
  f       <- comp_info$file[i]
  method  <- comp_info$method[i]
  quarter <- comp_info$quarter[i]
  
  message("Evaluating: ", basename(f))
  
  R <- try(rast(f), silent = TRUE)
  if (inherits(R, "try-error")) next
  
  R <- try(project(R, target_crs), silent = TRUE)
  if (inherits(R, "try-error")) next
  
  scores_all <- try(
    terra::extract(R, terra::vect(all_sf), ID = FALSE)[,1],
    silent = TRUE
  )
  
  labels_all <- all_df$resp
  
  for (k in seq_len(k_folds)) {
    
    idx <- which(fold_id == k)
    
    # defaults (ALWAYS length 1)
    m_aucpr <- m_boyce <- m_tss <- m_thr <- m_proc <- cal_slope <- cal_int <- NA_real_
    
    if (!inherits(scores_all, "try-error") &&
        length(idx) >= min_test_points) {
      
      y <- labels_all[idx]
      s <- scores_all[idx]
      
      if (length(unique(y)) == 2 &&
          sum(y == 1) >= min_pos_fold &&
          sum(y == 0) >= min_bg_fold) {
        
        ok <- is.finite(s)
        
        if (sum(ok) >= min_test_points) {
          
          s2 <- s[ok]
          y2 <- y[ok]
          
          m_aucpr <- safe_scalar(auc_pr(s2, y2))
          
          tss <- tss_best(s2, y2)
          m_tss <- safe_scalar(tss$TSS)
          m_thr <- safe_scalar(tss$thr)
          
          m_boyce <- safe_scalar(
            boyce_index(s2[y2 == 1], s2[y2 == 0])
          )
          
          m_proc <- safe_scalar(
            partial_roc_ratio(s2, y2, fpr_max = fpr_max)
          )
          
          cal <- calibration_summary(s2, y2, n_bins = 10)
          if (!is.null(cal) && nrow(cal) >= 3) {
            fit <- lm(obs_rate ~ mean_pred, data = cal)
            cal_slope <- safe_scalar(coef(fit)[2])
            cal_int   <- safe_scalar(coef(fit)[1])
          }
        }
      }
    }
    
    # ✅ THIS CAN NO LONGER FAIL
    results[[length(results) + 1]] <- data.frame(
      method     = method,
      quarter    = quarter,
      fold       = k,
      AUC_PR     = m_aucpr,
      TSS        = m_tss,
      TSS_thr    = m_thr,
      pROC_ratio = m_proc,
      Cal_slope  = cal_slope,
      Cal_int    = cal_int,
      stringsAsFactors = FALSE
    )
  }
}

res_df <- dplyr::bind_rows(results)
write.csv(res_df, out_fold_csv, row.names = FALSE)

cat("Saved fold metrics to:", out_fold_csv, "\n")



# ------------------------------------------------------------------------------
# 13) Summarise metrics + CI (from spatial CV folds)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Helper: Stand-alone Boyce index (Option B, moving window)
# ------------------------------------------------------------------------------

safe_boyce_optionB <- function(pred_bg, pred_pres,
                               nclass = 20,
                               window.w = 0.10,
                               eps = 1e-6,
                               drop_zeros = TRUE,
                               verbose = FALSE) {
  
  pred_bg   <- pred_bg[is.finite(pred_bg)]
  pred_pres <- pred_pres[is.finite(pred_pres)]
  
  if (length(pred_bg) < 50 || length(pred_pres) < 10) {
    return(list(Boyce = NA_real_, why = "too few points", raw = NULL))
  }
  
  if (drop_zeros) {
    pred_bg2   <- pred_bg[pred_bg > 0]
    pred_pres2 <- pred_pres[pred_pres > 0]
    why0 <- "drop_zeros=TRUE"
  } else {
    pred_bg2   <- pred_bg + eps
    pred_pres2 <- pred_pres + eps
    why0 <- "drop_zeros=FALSE"
  }
  
  if (length(pred_bg2) < 50 || length(pred_pres2) < 10) {
    return(list(Boyce = NA_real_, why = paste("too many zeros;", why0), raw = NULL))
  }
  
  b <- try(
    ecospat::ecospat.boyce(
      fit      = pred_bg2,
      obs      = pred_pres2,
      nclass   = nclass,
      window.w = window.w,
      PEplot   = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(b, "try-error") || is.null(b)) {
    return(list(Boyce = NA_real_, why = "ecospat.boyce error", raw = b))
  }
  
  boyce_val <- NA_real_
  if ("Spearman.cor" %in% names(b)) boyce_val <- as.numeric(b$Spearman.cor)
  if (!is.finite(boyce_val) && "cor" %in% names(b)) boyce_val <- as.numeric(b$cor)
  
  if (!is.finite(boyce_val)) {
    return(list(Boyce = NA_real_, why = "Boyce not computable", raw = b))
  }
  
  list(Boyce = boyce_val, why = "OK", raw = b)
}
library(dplyr)

ci_mean <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(c(NA_real_, NA_real_))
  se <- sd(x) / sqrt(length(x))
  z  <- qnorm((1 + level) / 2)
  c(mean(x) - z * se, mean(x) + z * se)
}

summary_df <- res_df %>%
  group_by(method, quarter) %>%
  summarise(
    AUC_PR_mean = ifelse(all(is.na(AUC_PR)), NA_real_, mean(AUC_PR, na.rm = TRUE)),
    AUC_PR_sd   = ifelse(all(is.na(AUC_PR)), NA_real_, sd(AUC_PR, na.rm = TRUE)),
    AUC_PR_CI_l = ci_mean(AUC_PR)[1],
    AUC_PR_CI_u = ci_mean(AUC_PR)[2],
    
    TSS_mean    = ifelse(all(is.na(TSS)), NA_real_, mean(TSS, na.rm = TRUE)),
    TSS_CI_l    = ci_mean(TSS)[1],
    TSS_CI_u    = ci_mean(TSS)[2],
    
    pROC_mean   = ifelse(all(is.na(pROC_ratio)), NA_real_, mean(pROC_ratio, na.rm = TRUE)),
    pROC_CI_l   = ci_mean(pROC_ratio)[1],
    pROC_CI_u   = ci_mean(pROC_ratio)[2],
    
    Cal_slope_mean = ifelse(all(is.na(Cal_slope)), NA_real_, mean(Cal_slope, na.rm = TRUE)),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 13b) Stand-alone Boyce index (FULL DATA, Option B — NOT per fold)
# ------------------------------------------------------------------------------

library(ecospat)

boyce_results <- list()

for (i in seq_len(nrow(comp_info))) {
  
  f       <- comp_info$file[i]
  method  <- comp_info$method[i]
  quarter <- comp_info$quarter[i]
  
  message("Computing Boyce (full model): ", basename(f))
  
  R <- try(rast(f), silent = TRUE)
  if (inherits(R, "try-error")) next
  
  R <- try(project(R, target_crs), silent = TRUE)
  if (inherits(R, "try-error")) next
  
  pred_pres <- try(terra::extract(R, pts_presence)[,2], silent = TRUE)
  pred_bg   <- try(terra::extract(R, bg_points)[,2], silent = TRUE)
  
  if (inherits(pred_pres, "try-error") || inherits(pred_bg, "try-error")) {
    boyce_results[[length(boyce_results) + 1]] <- data.frame(
      method = method,
      quarter = quarter,
      Boyce = NA_real_,
      Boyce_status = "extraction error",
      stringsAsFactors = FALSE
    )
    next
  }
  
  boy <- safe_boyce_optionB(
    pred_bg   = pred_bg,
    pred_pres = pred_pres,
    nclass    = 20,
    window.w  = 0.10,
    drop_zeros = TRUE,
    verbose   = FALSE
  )
  
  boyce_results[[length(boyce_results) + 1]] <- data.frame(
    method = method,
    quarter = quarter,
    Boyce = boy$Boyce,
    Boyce_status = boy$why,
    stringsAsFactors = FALSE
  )
}

boyce_df <- bind_rows(boyce_results)

# ------------------------------------------------------------------------------
# 13c) Join Boyce + compute FINAL SCORE
# ------------------------------------------------------------------------------

summary_df <- summary_df %>%
  left_join(boyce_df, by = c("method", "quarter")) %>%
  mutate(
    Score = ifelse(
      is.finite(AUC_PR_mean) &
        is.finite(TSS_mean) &
        is.finite(pROC_mean) &
        is.finite(Boyce),
      
      0.20 * AUC_PR_mean +
        0.20 * TSS_mean +
        0.20 * pROC_mean +
        0.40 * Boyce,   # ecological consistency weighted higher
      
      NA_real_
    )
  ) %>%
  arrange(desc(Score))

# ------------------------------------------------------------------------------
# Save FINAL publication table
# ------------------------------------------------------------------------------

write.csv(summary_df, out_summary_csv, row.names = FALSE)

cat("\nSaved FINAL summary table with CI + Boyce to:\n",
    out_summary_csv, "\n")

