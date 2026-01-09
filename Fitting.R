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

# Your data
###################Eggs developmmental rates
df <- data.frame(
  Temperature = c(25.5, 26.5, 29.5, 23, 28),
  SurvivalRate = c(0.164203612, 0.179533214, 0.183150183,
                   0.09596929, 0.1447178)
)

# Fit model: mu(T) = exp(a + b*T + c*T^2)
fit <- nls(
  SurvivalRate ~ exp(a + b * Temperature + c * Temperature^2),
  data = df,
  start = list(a = -5, b = 0.1, c = -0.002),
  control = nls.control(maxiter = 200, warnOnly = TRUE)
)

summary(fit)
coef(fit)
df$Predicted <- predict(fit)
df


# Observed vs predicted
obs <- df$SurvivalRate
pred <- predict(fit)

# Residuals
res <- obs - pred

# ---- RMSE ----
RMSE <- sqrt(mean(res^2))

# ---- MAE ----
MAE <- mean(abs(res))

# ---- R-squared ----
SST <- sum((obs - mean(obs))^2)
SSR <- sum(res^2)
R2 <- 1 - SSR/SST

# ---- AIC and BIC ----
n <- length(obs)
k <- length(coef(fit))
AIC_val <- n * log(SSR/n) + 2 * k
BIC_val <- n * log(SSR/n) + k * log(n)

# Display results
list(
  RMSE = RMSE,
  MAE = MAE,
  R2 = R2,
  AIC = AIC_val,
  BIC = BIC_val
)
#######################Larva developmental rates
df <- data.frame(
  Temperature = c(25.5, 26.5, 29.5),
  SurvivalRate = c(0.039323634, 0.047147572, 0.06779661
                   )
)

fit <- nls(
  SurvivalRate ~ exp(a + b * Temperature),
  data = df,
  start = list(a = -3, b = 0.1),
  control = nls.control(maxiter = 200, warnOnly = TRUE)
)

summary(fit)


summary(fit)
coef(fit)
df$Predicted <- predict(fit)
df

# Observed vs predicted
obs <- df$SurvivalRate
pred <- predict(fit)

# Residuals
res <- obs - pred

# ---- RMSE ----
RMSE <- sqrt(mean(res^2))

# ---- MAE ----
MAE <- mean(abs(res))

# ---- R-squared ----
SST <- sum((obs - mean(obs))^2)
SSR <- sum(res^2)
R2 <- 1 - SSR/SST

# ---- AIC and BIC ----
n <- length(obs)
k <- length(coef(fit))
AIC_val <- n * log(SSR/n) + 2 * k
BIC_val <- n * log(SSR/n) + k * log(n)

# Display results
list(
  RMSE = RMSE,
  MAE = MAE,
  R2 = R2,
  AIC = AIC_val,
  BIC = BIC_val
)
###########################Pupae developmental rates
0.057636888
0.11778563
0.154798762
df <- data.frame(
  Temperature = c(25.5, 26.5, 29.5),
  SurvivalRate = c(0.057636888, 0.11778563, 0.154798762
  )
)

fit <- nls(
  SurvivalRate ~ exp(a + b * Temperature),
  data = df,
  start = list(a = -3, b = 0.1),
  control = nls.control(maxiter = 200, warnOnly = TRUE)
)

summary(fit)


summary(fit)
coef(fit)
df$Predicted <- predict(fit)
df


# Observed vs predicted
obs <- df$SurvivalRate
pred <- predict(fit)

# Residuals
res <- obs - pred

# ---- RMSE ----
RMSE <- sqrt(mean(res^2))

# ---- MAE ----
MAE <- mean(abs(res))

# ---- R-squared ----
SST <- sum((obs - mean(obs))^2)
SSR <- sum(res^2)
R2 <- 1 - SSR/SST

# ---- AIC and BIC ----
n <- length(obs)
k <- length(coef(fit))
AIC_val <- n * log(SSR/n) + 2 * k
BIC_val <- n * log(SSR/n) + k * log(n)

# Display results
list(
  RMSE = RMSE,
  MAE = MAE,
  R2 = R2,
  AIC = AIC_val,
  BIC = BIC_val
)
##################eggs_mortality
# Data
df <- data.frame(
  T = c(15, 18, 20, 25, 28, 32),
  M = c(0.72, 0.40, 0.30, 0.38, 0.28, 0.65)
)


gauss_simple  <- nls(
  M ~ 1 / (1 + exp(-(a + b*T + c*T^2))),
  data = df,
  start = list(a = -3, b = 0.2, c = -0.01),
  algorithm = "port",
  lower = c(a = -10, b = -5, c = -1),
  upper = c(a = 10, b = 5, c = 1)
)

# Model summary
summary(gauss_simple)

# Coefficients
coef(gauss_simple)

# Predicted values
df$Predicted <- predict(gauss_simple)

# Observed vs predicted
obs <- df$M
pred <- df$Predicted

# Residuals
res <- obs - pred

# ---- RMSE ----
RMSE <- sqrt(mean(res^2))

# ---- MAE ----
MAE <- mean(abs(res))

# ---- R-squared ----
SST <- sum((obs - mean(obs))^2)
SSR <- sum(res^2)
R2 <- 1 - SSR / SST

# ---- AIC and BIC ----
n <- length(obs)
k <- length(coef(gauss_simple))
AIC_val <- n * log(SSR / n) + 2 * k
BIC_val <- n * log(SSR / n) + k * log(n)

# Display results
list(
  RMSE = RMSE,
  MAE = MAE,
  R2 = R2,
  AIC = AIC_val,
  BIC = BIC_val
)
###########Larvae mortality rates

# Data
df <- data.frame(
  T = c(15, 18, 20, 25, 28, 32),
  M = c(0.95, 0.85, 0.50, 0.7, 0.30, 0.25)
)


gauss_simple  <- nls(
  M ~ 1 / (1 + exp(-(a + b*T + c*T^2))),
  data = df,
  start = list(a = -3, b = 0.2, c = -0.01),
  algorithm = "port",
  lower = c(a = -10, b = -5, c = -1),
  upper = c(a = 10, b = 5, c = 1)
)

# Model summary
summary(gauss_simple)

# Coefficients
coef(gauss_simple)

# Predicted values
df$Predicted <- predict(gauss_simple)

# Observed vs predicted
obs <- df$M
pred <- df$Predicted

# Residuals
res <- obs - pred

# ---- RMSE ----
RMSE <- sqrt(mean(res^2))

# ---- MAE ----
MAE <- mean(abs(res))

# ---- R-squared ----
SST <- sum((obs - mean(obs))^2)
SSR <- sum(res^2)
R2 <- 1 - SSR / SST

# ---- AIC and BIC ----
n <- length(obs)
k <- length(coef(gauss_simple))
AIC_val <- n * log(SSR / n) + 2 * k
BIC_val <- n * log(SSR / n) + k * log(n)

# Display results
list(
  RMSE = RMSE,
  MAE = MAE,
  R2 = R2,
  AIC = AIC_val,
  BIC = BIC_val
)



######################Pupa_Mortality_rates	

# Data
df <- data.frame(
  T = c(15, 18, 20, 25, 28, 32),
  M = c(0.5, 0.1, 0.15, 0.1, 0.25, 0.3)
)


gauss_simple  <- nls(
  M ~ 1 / (1 + exp(-(a + b*T + c*T^2))),
  data = df,
  start = list(a = -3, b = 0.2, c = -0.01),
  algorithm = "port",
  lower = c(a = -10, b = -5, c = -1),
  upper = c(a = 10, b = 5, c = 1)
)

# Model summary
summary(gauss_simple)

# Coefficients
coef(gauss_simple)

# Predicted values
df$Predicted <- predict(gauss_simple)

# Observed vs predicted
obs <- df$M
pred <- df$Predicted

# Residuals
res <- obs - pred

# ---- RMSE ----
RMSE <- sqrt(mean(res^2))

# ---- MAE ----
MAE <- mean(abs(res))

# ---- R-squared ----
SST <- sum((obs - mean(obs))^2)
SSR <- sum(res^2)
R2 <- 1 - SSR / SST

# ---- AIC and BIC ----
n <- length(obs)
k <- length(coef(gauss_simple))
AIC_val <- n * log(SSR / n) + 2 * k
BIC_val <- n * log(SSR / n) + k * log(n)

# Display results
list(
  RMSE = RMSE,
  MAE = MAE,
  R2 = R2,
  AIC = AIC_val,
  BIC = BIC_val
)

############################ Adults mortality rates
# Data
df <- data.frame(
  T = c(15, 18, 20, 25, 28, 32),
  M = c(0.052521008, 0.085910653, 0.10080645, 0.129533679, 0.099206349, 0.173611111)
)


gauss_simple  <- nls(
  M ~ 1 / (1 + exp(-(a + b*T + c*T^2))),
  data = df,
  start = list(a = -3, b = 0.2, c = -0.01),
  algorithm = "port",
  lower = c(a = -10, b = -5, c = -1),
  upper = c(a = 10, b = 5, c = 1)
)

# Model summary
summary(gauss_simple)

# Coefficients
coef(gauss_simple)

# Predicted values
df$Predicted <- predict(gauss_simple)

# Observed vs predicted
obs <- df$M
pred <- df$Predicted

# Residuals
res <- obs - pred

# ---- RMSE ----
RMSE <- sqrt(mean(res^2))

# ---- MAE ----
MAE <- mean(abs(res))

# ---- R-squared ----
SST <- sum((obs - mean(obs))^2)
SSR <- sum(res^2)
R2 <- 1 - SSR / SST

# ---- AIC and BIC ----
n <- length(obs)
k <- length(coef(gauss_simple))
AIC_val <- n * log(SSR / n) + 2 * k
BIC_val <- n * log(SSR / n) + k * log(n)

# Display results
list(
  RMSE = RMSE,
  MAE = MAE,
  R2 = R2,
  AIC = AIC_val,
  BIC = BIC_val
)

############################
# Fecundity (Gaussian model)
############################

# Input data
df_fec <- data.frame(
  T = c(25.5, 26.5, 29.5),
  F = c(42.49, 42.71, 41.30)
)

# Fit Gaussian thermal fecundity curve:
# F(T) = A * exp(- (T - mu)^2 / (2 * sigma^2))
gauss_fec <- nls(
  F ~ A * exp(-((T - mu)^2) / (2 * sigma^2)),
  data = df_fec,
  start = list(A = 45, mu = 27, sigma = 2),
  algorithm = "port",
  lower = c(A = 0,  mu = 20, sigma = 0.1),
  upper = c(A = 100, mu = 35, sigma = 10),
  control = nls.control(maxiter = 200, warnOnly = TRUE)
)

# Model summary
summary(gauss_fec)

# Estimated parameters
coef(gauss_fec)

# Predicted fecundity
df_fec$Predicted <- predict(gauss_fec)
df_fec
# Observed vs predicted
obs  <- df_fec$F
pred <- df_fec$Predicted

# Residuals
res <- obs - pred

# RMSE
RMSE <- sqrt(mean(res^2))

# MAE
MAE <- mean(abs(res))

# R-squared
SST <- sum((obs - mean(obs))^2)
SSR <- sum(res^2)
R2  <- 1 - SSR / SST

# Information criteria
n <- length(obs)
k <- length(coef(gauss_fec))
AIC_val <- n * log(SSR / n) + 2 * k
BIC_val <- n * log(SSR / n) + k * log(n)

# Output diagnostics
list(
  RMSE = RMSE,
  MAE  = MAE,
  R2   = R2,
  AIC  = AIC_val,
  BIC  = BIC_val
)

