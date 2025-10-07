library(tidyverse)
library(mgcv)
library(patchwork)

set.seed(123)
N <- 100
x_coord  <- runif(N) * 100
y_coord  <- runif(N) * 100
XY <- cbind(x_coord, y_coord)

# --- helper: Matérn covariance (ν = 1.5; smooth but not too smooth) -------------
matern32 <- function(d, rho) {
  # ν = 3/2 Matérn: (1 + sqrt(3) d/rho) exp(-sqrt(3) d/rho)
  a <- sqrt(3) * d / rho
  (1 + a) * exp(-a)
}

# --- build covariance matrices for two fields (different ranges) ----------------
D <- as.matrix(dist(XY))
sigma2_g <- 1.0   # process variance for green
sigma2_a <- 1.0   # process variance for albedo
tau2_g   <- 0.05  # nugget (measurement noise)
tau2_a   <- 0.10
rho_g    <- 40    # spatial range (larger => smoother over space)
rho_a    <- 20

Cg <- sigma2_g * matern32(D, rho_g) + diag(tau2_g, N)
Ca <- sigma2_a * matern32(D, rho_a) + diag(tau2_a, N)

# --- simulate Gaussian random fields via Cholesky --------------------------------
# add a tiny jitter for numerical stability if needed
Cg <- Cg + diag(1e-10, N)
Ca <- Ca + diag(1e-10, N)

Lg <- chol(Cg)
La <- chol(Ca)

zg <- as.numeric(t(Lg) %*% rnorm(N))
za <- as.numeric(t(La) %*% rnorm(N))

# map to [0,1] (optional; keeps the "feel" of your original scales)
green  <- (zg - min(zg)) / diff(range(zg))
albedo <- (za - min(za)) / diff(range(za))

# response with *predictor* effects + noise (still non-spatial if you want)
# also can change these if you want one to have a stronger or weaker effect
beta_green <- 2
beta_albedo <- 1.5
air_temp <- beta_green * green + beta_albedo * albedo + rnorm(N, sd = 0.5)

# make output df
df <- tibble(air_temp, green, albedo, x_coord, y_coord)




