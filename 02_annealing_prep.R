#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Prepare for the simulated annealing process
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
library(tidyverse)

S <- rep(0, N)  # the blank matrix
k <- 5          # the first number of monitors to start with

# importantly you can now compare this output with the fortran output for score
system("R CMD SHLIB simann.f90")

dyn.load("simann.so")

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# probably the way to do this in FORTRAN is to convert df to wide first
# because then its just rows and columns
# and probably becomes easier to subset
df_wide <- df %>% pivot_wider(id_cols = 'monitor_id', names_from = 'day_id',
                              values_from = 'air_temp') %>% 
  select(-monitor_id) %>% as.matrix()

# get true values by day
# and this is for the response
# daily NETWORK 50th percentile, 5th percentile, and 95th percentile
# choosing these and not min and max so its a little more robust

#' @param x a matrix
#' @param q a quantile
get_metric <- function(x, q) apply(x, 2, function(a) quantile(a, q))

# check that this works against the Fortran version
check_r1 <- get_metric(df_wide, 0.50)

check_f1 <- .Fortran('get_metric', 
                     x = df_wide, 
                     nrow = as.integer(N), 
                     ncol = as.integer(N_daymet), 
                     q = 0.5, 
                     xout = vector("numeric", N_daymet))

stopifnot(mean(abs(check_r1 - check_f1$xout)) < 0.001)

# now wrap into a bigger function
get_metrics <- function(df_sub) {
  stat1 <- get_metric(df_sub, 0.05)
  stat2 <- get_metric(df_sub, 0.50)
  stat3 <- get_metric(df_sub, 0.95)
  return(cbind(stat1, stat2, stat3))
}

# get the true values to compare against later
df_best <- get_metrics(df_wide)

# function for MSE
get_mse <- function(a, b) mean((a - b)^2)

# check that its comparable
# get_mse(a, b, n, mse)
set.seed(123)
a = rnorm(10)
b = rnorm(10)
check_f3 <- .Fortran('get_mse', 
                     a = a, 
                     b = b, 
                     n = as.integer(10), 
                     mse = 0.)
stopifnot(mean(abs(check_f3$mse - get_mse(a, b))) < 0.001)


# This is the R version of your score matrix
# you can use it to check your math
# for now, lets say is the sum of squared errors for 
# daily NETWORK mean, 5th percentile, and 95th percentile
get_score <- function(S_local) {
  
  # S_local = S
  
  S_ones <- which(S_local == 1)
  
  # filter the data frame
  df_sub <- df_wide[S_ones, ]

  # get just the points for this sample
  yy <- get_metrics(df_sub)

  # get day-wise errors
  zz = vector("numeric", ncol(df_best))

  for(i in 1:(ncol(df_best))) {
    a = unlist(unname(as.vector(df_best[, i])))
    b = unlist(unname(as.vector(yy[, i])))
    zz[i] = get_mse(a, b)
  }

  # get sum
  # since we want this to be low, you don't have to make it negative
  sum(zz)

}

# get the first set
S[c(1:k)] <- 1

# get the score
check_r2 <- get_score(S)

# compare
# S, df_wide, magic_n, nsites, ndays, SCORE
# get_score(S, df_wide, magic_n, nsites, ndays, score_cols, SCORE)
check_f2 <- .Fortran('get_score', 
                     S = as.integer(S),
                     df_wide = df_wide,
                     magic_n = as.integer(k),
                     nsites = as.integer(N),
                     ndays = as.integer(N_daymet),
                     score_cols = as.integer(3),
                     SCORE = 0.)


stopifnot(mean(abs(check_r2 - check_f2$SCORE)) < 0.001)

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////








