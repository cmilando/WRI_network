#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' This script creates the simulated data
#' with some minor covariance
#' 
#' change the beta coefficients to obeserve a different effect later
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
library(tidyverse)
library(mgcv)
library(patchwork)

set.seed(123)
N <- 100             # number of monitors
beta_green  <- 2     # beta coefficient for greenness
beta_albedo <- 1.5   # beta coefficient for albedo

N_daymet    <- 20    # number of days of daily meteorology
beta_daymet <- 10    # beta coefficient for the meteorology parameter

air_temp_sd <- 2     # st dev for random error, increase for more noise

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
x_coord  <- runif(N) * 100
y_coord  <- runif(N) * 100
XY <- cbind(x_coord, y_coord)
monitor_points <- data.frame(XY)
monitor_points$id <- 1:N
dim(monitor_points)

# --- helper: Matérn covariance (ν = 1.5; smooth but not too smooth) ----------
matern32 <- function(d, rho) {
  # ν = 3/2 Matérn: (1 + sqrt(3) d/rho) exp(-sqrt(3) d/rho)
  a <- sqrt(3) * d / rho
  (1 + a) * exp(-a)
}

# --- build covariance matrices for two fields (different ranges) -------------
D <- as.matrix(dist(XY))
dim(D)

# process variance for green
sigma2_g <- 1.0   
sigma2_a <- 1.0   
sigma2_d <- 1.0   

# nugget (measurement noise)
tau2_g   <- 0.05  
tau2_a   <- 0.10
tau2_d   <- 0.75

# spatial range (larger => smoother over space)
rho_g    <- 40    
rho_a    <- 20
rho_d    <- 30

# build the covariance matrix
# add a tiny jitter for numerical stability if needed
Cg <- sigma2_g * matern32(D, rho_g) + diag(tau2_g, N) + diag(1e-10, N)
Ca <- sigma2_a * matern32(D, rho_a) + diag(tau2_a, N) + diag(1e-10, N)
Cd_l <- vector("list", N_daymet) 
for(i in 1:N_daymet) {
  Cd_l[[i]] <- sigma2_d * matern32(D, rho_d) + diag(tau2_d, N) + diag(1e-10, N)
}


# --- simulate Gaussian random fields via Cholesky ----------------------------
Lg <- chol(Cg)
La <- chol(Ca)
Ld_l <- lapply(Cd_l, chol)
length(La)        
#
zg <- as.numeric(t(Lg) %*% rnorm(N))
za <- as.numeric(t(La) %*% rnorm(N))
length(za)
zd_l <- lapply(Ld_l, \(l) as.numeric(t(l) %*% rnorm(N)))

# map to [0,1] (optional; keeps the "feel" of your original scales)
green  <- (zg - min(zg)) / diff(range(zg))
albedo <- (za - min(za)) / diff(range(za))
daymet_l <- lapply(zd_l, \(z) (z - min(z)) / diff(range(z)))
length(daymet_l[[1]])

# repeat N * daymet times
green  <- rep(green, times = N_daymet)
albedo <- rep(albedo, times = N_daymet)
id <- rep(1:N, times = N_daymet)
length(green)
length(albedo)

# also add in daymet which
daymet <- do.call(c, daymet_l)
length(daymet)

X_mat <- cbind(daymet, green, albedo)
dim(X_mat)
beta_mat <- c(beta_daymet, beta_green, beta_albedo)

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
# response with *predictor* effects + noise (still non-spatial if you want)
# also can change these if you want one to have a stronger or weaker effect
air_temp <- X_mat %*% beta_mat + rnorm(N * daymet, sd = air_temp_sd)
air_temp
dim(air_temp)

# make output df
df <- tibble(air_temp = as.vector(air_temp), green, albedo, daymet, id)
df <- left_join(df, monitor_points, by = join_by(id))
df$day_id <- rep(1:N_daymet, each = N)

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
# --- plot --------------------------------
# quick look at spatial signal
p1 <- ggplot(df %>% filter(day_id == 1), 
             aes(x_coord, y_coord, fill = green))  + 
  geom_point(size = 2, shape = 21) + 
  scale_fill_viridis_c() + ggtitle("a. Measured Green") + 
  theme(legend.position = 'bottom')

p2 <- ggplot(df %>% filter(day_id == 1), 
             aes(x_coord, y_coord, fill = albedo)) + 
  geom_point(size = 2, shape = 21) + 
  scale_fill_viridis_c() + ggtitle("b. Measured Albedo") + 
  theme(legend.position = 'bottom')

p3 <- ggplot(df %>% filter(day_id == 1), 
             aes(x_coord, y_coord, fill = daymet)) + 
  geom_point(size = 2, shape = 21) + 
  scale_fill_viridis_c() + ggtitle("c. DayMet (day = 1)") + 
  theme(legend.position = 'bottom')

p4 <- ggplot(df %>% filter(day_id == 1), 
             aes(x_coord, y_coord, fill = air_temp)) + 
  geom_point(size = 2, shape = 21) + 
  scale_fill_viridis_c() + ggtitle("d. Measured Air temp (day = 1)") + 
  theme(legend.position = 'bottom')

(p1 + p2) / (p3 + p4) 
ggsave("img/fig1.png", width = 8, height = 7, dpi = 300)


