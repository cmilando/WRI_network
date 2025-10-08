#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Prepare for the simulated annealing process
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
library(tidyverse)

S <- rep(0, N)  # the blank matrix
k <- 5          # the first number of monitors to start with

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# probably the way to do this in FORTRAN is to convert df to wide first
# because then its just rows and columns
# and probably becomes easier to subset
df_wide <- df %>% pivot_wider(id_cols = 'id', names_from = 'day_id',
                              values_from = 'air_temp') %>% 
  select(-id) %>% as.matrix()

# get true values by day
# and this is for the response
# daily NETWORK 50th percentile, 5th percentile, and 95th percentile
# choosing these and not min and max so its a little more robust

#' @param x a matrix
#' @param q a quantile
get_metric <- function(x, q) apply(x, 2, function(a) quantile(a, q))

get_metrics <- function(df_sub) {
  stat1 <- get_metric(df_sub, 0.50)
  stat2 <- get_metric(df_sub, 0.05)
  stat3 <- get_metric(df_sub, 0.95)
  return(cbind(stat1, stat2, stat3))
}

# get the true values to compare against later
df_best <- get_metrics(df_wide)

# function for MSE
get_mse <- function(a, b) mean((a - b)^2)

# This is the R version of your score matrix
# you can use it to check your math
# for now, lets say is the sum of squared errors for 
# daily NETWORK mean, 5th percentile, and 95th percentile
get_score <- function(S_local) {

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
  -1*sum(zz)

}

# get the first set
S[c(1:k)] <- 1

# get the score
get_score(S)

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# importantly you can now compare this output with the fortran output for score
system("R CMD SHLIB simann.f90")
dyn.unload("simann.so")
dyn.load("simann.so")






