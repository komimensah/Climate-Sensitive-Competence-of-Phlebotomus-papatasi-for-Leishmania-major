library(terra)
library(ggplot2)

# ==========================================================
# 0) WORKING DIRECTORY
# ==========================================================
setwd("/Users/kagboka/Desktop/socioeeconomic_modelling_review/Violin plots")

# ==========================================================
# 1) INPUT DIRECTORIES
# ==========================================================
present_dir <- "present"
ssp245_dir  <- "ssp245"
ssp585_dir  <- "ssp585"
quarters <- c("Q1","Q2","Q3","Q4")

out_dir <- "deltas"
dir.create(out_dir, showWarnings = FALSE)

# ==========================================================
# 2) ENDEMIC SITES (CSV → SpatVector)
# ==========================================================
obs_path <- "/Users/kagboka/Desktop/socioeeconomic_modelling_review/Final data DIsease.csv"
obs_df   <- read.csv(obs_path)

sites <- vect(
  obs_df,
  geom = c("lon", "lat"),
  crs  = "EPSG:4326"
)

# ==========================================================
# 3) RASTER LOADERS (EXACT FILENAMES)
# ==========================================================
load_present <- function(q) {
  f <- file.path(present_dir, paste0(q, ".tif"))
  if (!file.exists(f)) stop("Missing present file: ", f)
  rast(f)
}

load_future <- function(dir, q, band = c("mean","lwr","upr")) {
  band <- match.arg(band)
  fname <- switch(
    band,
    mean = paste0(q, ".tif"),
    lwr  = paste0("lwr_", q, ".tif"),
    upr  = paste0("upr_", q, ".tif")
  )
  f <- file.path(dir, fname)
  if (!file.exists(f)) stop("Missing future file: ", f)
  rast(f)
}

# ==========================================================
# 4) MATCHED RANK-BISERIAL EFFECT SIZE (PAIRED)
# ==========================================================
rank_biserial <- function(present, future) {
  d <- future - present
  d <- d[!is.na(d) & d != 0]
  n <- length(d)
  if (n < 10) return(NA_real_)
  r <- rank(abs(d))
  Wp <- sum(r[d > 0])
  T  <- n * (n + 1) / 2
  (2 * Wp - T) / T
}

# ==========================================================
# 5) BOOTSTRAP CI FOR r_rb
# ==========================================================
boot_ci_r_rb <- function(present, future, R = 500, seed = 123) {
  
  keep <- complete.cases(present, future)
  present <- present[keep]
  future  <- future[keep]
  
  n0 <- length(present)
  if (n0 < 30) return(c(NA_real_, NA_real_))
  
  set.seed(seed)
  rb <- numeric(R)
  
  for (i in seq_len(R)) {
    idx <- sample.int(n0, replace = TRUE)
    rb[i] <- rank_biserial(present[idx], future[idx])
  }
  
  rb <- rb[is.finite(rb)]
  quantile(rb, probs = c(0.025, 0.975), na.rm = TRUE)
}

# ==========================================================
# 6) MAIN LOOP
# ==========================================================
results   <- list()
plot_data <- list()

for (q in quarters) {
  
  cat("Processing", q, "\n")
  
  # ---- Load rasters ----
  r_pres <- load_present(q)
  
  r245_m <- load_future(ssp245_dir, q, "mean")
  r245_l <- load_future(ssp245_dir, q, "lwr")
  r245_u <- load_future(ssp245_dir, q, "upr")
  
  r585_m <- load_future(ssp585_dir, q, "mean")
  r585_l <- load_future(ssp585_dir, q, "lwr")
  r585_u <- load_future(ssp585_dir, q, "upr")
  
  # ---- Project sites ----
  sites_p <- project(sites, crs(r_pres))
  
  # ---- Extract values ----
  pres  <- extract(r_pres,  sites_p)[,2]
  
  f245m <- extract(r245_m, sites_p)[,2]
  f245l <- extract(r245_l, sites_p)[,2]
  f245u <- extract(r245_u, sites_p)[,2]
  
  f585m <- extract(r585_m, sites_p)[,2]
  f585l <- extract(r585_l, sites_p)[,2]
  f585u <- extract(r585_u, sites_p)[,2]
  
  # ---- Paired deltas ----
  keep245 <- complete.cases(pres, f245m)
  keep585 <- complete.cases(pres, f585m)
  
  d245m <- f245m[keep245] - pres[keep245]
  d245l <- f245l[keep245] - pres[keep245]
  d245u <- f245u[keep245] - pres[keep245]
  
  d585m <- f585m[keep585] - pres[keep585]
  d585l <- f585l[keep585] - pres[keep585]
  d585u <- f585u[keep585] - pres[keep585]
  
  # ---- Paired Wilcoxon (mean only) ----
  w245 <- wilcox.test(d245m, mu = 0, exact = FALSE)
  w585 <- wilcox.test(d585m, mu = 0, exact = FALSE)
  
  # ---- Effect sizes + CI ----
  rb245 <- rank_biserial(pres[keep245], f245m[keep245])
  rb585 <- rank_biserial(pres[keep585], f585m[keep585])
  
  ci245 <- boot_ci_r_rb(pres[keep245], f245m[keep245])
  ci585 <- boot_ci_r_rb(pres[keep585], f585m[keep585])
  
  # ---- Store results (FIXED COLUMN NAMES) ----
  results[[q]] <- rbind(
    data.frame(
      Quarter=q, Scenario="SSP2-4.5",
      N_sites=length(d245m),
      p_value=w245$p.value,
      r_rb=rb245, CI_low=ci245[1], CI_high=ci245[2]
    ),
    data.frame(
      Quarter=q, Scenario="SSP5-8.5",
      N_sites=length(d585m),
      p_value=w585$p.value,
      r_rb=rb585, CI_low=ci585[1], CI_high=ci585[2]
    )
  )
  
  # ---- Store for plotting ----
  plot_data[[q]] <- rbind(
    data.frame(Quarter=q, Scenario="SSP2-4.5", Delta=d245m),
    data.frame(Quarter=q, Scenario="SSP5-8.5", Delta=d585m)
  )
}

# ==========================================================
# 7) FINAL TABLE
# ==========================================================
summary_table <- do.call(rbind, results)
write.csv(summary_table, "Wilcoxon_DeltaSuitability_ByQuarter.csv", row.names = FALSE)
print(summary_table)


# ==========================================================
# 8) FIGURE — HALF VIOLIN Δ-SUITABILITY (NO gghalves)
# ==========================================================

GeomSplitViolin <- ggproto(
  "GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x))
    
    grp <- data[1, "group"]
    newdata <- data
    gap <- 0.04  # adjust between 0.02–0.06
    
    if (grp %% 2 == 1) {
      newdata$x <- newdata$xminv - gap
    } else {
      newdata$x <- newdata$xmaxv + gap
    }
    
    GeomPolygon$draw_panel(newdata, ...)
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL,
                              stat = "ydensity", position = "identity",
                              ..., draw_quantiles = NULL,
                              trim = TRUE, scale = "area",
                              na.rm = FALSE, show.legend = NA,
                              inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm,
      ...
    )
  )
}
plot_df$Scenario <- factor(
  plot_df$Scenario,
  levels = c("SSP2-4.5", "SSP5-8.5")
)

p <- ggplot(
  plot_df,
  aes(
    x = Quarter,
    y = Delta,
    fill = Scenario
  )
) +
  
  # ---- TRUE split violins ----
geom_split_violin(
  trim = FALSE,
  alpha = 0.7,
  color = "black"
) +
  
  # ---- Central boxplot ----
geom_boxplot(
  width = 0.12,
  outlier.shape = NA,
  alpha = 0.9
) +
  
  # ---- Mean ----
stat_summary(
  fun = mean,
  geom = "point",
  shape = 21,
  size = 2.5,
  fill = "white"
) +
  
  scale_fill_manual(
    values = c(
      "SSP2-4.5" = "#E76F51",
      "SSP5-8.5" = "#2A9D8F"
    )
  ) +
  
  labs(
    x = "Seasonal quarter",
    y = expression(Delta*" suitability (future − present)")
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(
  "Figure4_DeltaSuitability_SplitViolin_TRUE.png",
  p,
  width = 8,
  height = 5,
  dpi = 600
)

print(p)
