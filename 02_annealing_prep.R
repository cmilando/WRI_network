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

# system("rm *.so")
# system("rm *.o")
# 
dyn.unload("simann.so")
dyn.load("simann.so")


# now for the linear model section
predictor_cols <- c('daymet', 'green', 'albedo')
monitor_id_col <- c('monitor_id')
response_col   <- c('air_temp')

beta_vector <- c(beta_daymet, beta_green, beta_albedo) # same order as above

# this assigns the weights to the different components of the composite score
# position 1 is for predictors similarity
# position 2 is for the output of the linear model
lambda_vec <- c(1, 1)


#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# probably the way to do this in FORTRAN is to convert df to wide first
# because then its just rows and columns
# and probably becomes easier to subset
df_wide <- df %>% pivot_wider(id_cols = 'monitor_id', names_from = 'day_id',
                              values_from = 'air_temp') %>% 
  select(-monitor_id) %>% as.matrix()

# also have to do some treatment to the co-variates
N_predictors <- length(predictor_cols)
Y_matrix    <- as.matrix(df[, response_col])
X_matrix    <- as.matrix(df[, predictor_cols])
ID_vector   <- as.integer(unlist(df[, monitor_id_col]))


#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
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

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# now wrap into a bigger function
get_metrics <- function(df_sub) {
  stat1 <- get_metric(df_sub, 0.05)
  stat2 <- get_metric(df_sub, 0.50)
  stat3 <- get_metric(df_sub, 0.95)
  return(cbind(stat1, stat2, stat3))
}

# get the true values to compare against later
df_best <- get_metrics(df_wide)

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# function for MSE
get_rmse <- function(a, b) sqrt(mean((a - b)^2))

# check that its comparable
# get_rmse(a, b, n, mse)
set.seed(123)

a = rnorm(10)
b = rnorm(10)

check_f3 <- .Fortran('get_rmse', 
                     a = a, 
                     b = b, 
                     n = as.integer(10), 
                     mse = 0.)

stopifnot(mean(abs(check_f3$mse - get_rmse(a, b))) < 0.001)

check_f3 <- .Fortran('get_rmse', 
                     a = a, 
                     b = a, 
                     n = as.integer(10), 
                     rmse = 0.)

stopifnot(check_f3$rmse < 0.001)


#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# check OLS decomposition


set.seed(1)
n <- 100L; p <- 3L
X <- matrix(rnorm(n*p), n, p)
Y <- rnorm(n)

beta <- double(p); 
info <- integer(1)

## SVD version (robust)
rank <- integer(1)
out2 <- .Fortran("ols_svd",
                 X     = as.double(X),
                 Y     = as.double(Y),
                 n     = as.integer(n),
                 p     = as.integer(p),
                 rcond = as.double(1e-12),
                 beta  = beta,
                 rank  = rank,
                 info  = info)
out2$beta
coef(lm(Y ~ X + 0))

# nice looks good


#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# This is the R version of your score matrix
# you can use it to check your math
# for now, lets say is the sum of squared errors for 
# daily NETWORK mean, 5th percentile, and 95th percentile
get_score <- function(S_local) {
  
  # S_local = S
  
  S_ones <- which(S_local == 1)
  
  # ---------------
  # -- Part 1 -----
  # ---------------
  
  if (lambda_vec[1] != 0) {
    
    # filter the data frame
    df_sub <- df_wide[S_ones, ]
  
    # get just the points for this sample
    yy <- get_metrics(df_sub)
  
    # get day-wise errors
    zz = vector("numeric", ncol(df_best))
  
    for(i in 1:(ncol(df_best))) {
      a = unlist(unname(as.vector(df_best[, i])))
      b = unlist(unname(as.vector(yy[, i])))
      zz[i] = get_rmse(a, b)
    }
  
    # get sum
    # since we want this to be low, you don't have to make it negative
    z1 <- sum(zz)
    
  } else {
    
    z1 <- 0
    
  }
  
  # ---------------
  # -- Part 2 -----
  # ---------------
  
  if (lambda_vec[2] != 0) {
    
    m_sub_rows <- which(ID_vector %in% S_ones)
    
    Y_sub <- matrix(Y_matrix[m_sub_rows, ], ncol = 1)
    
    X_matrix_sub <- X_matrix[m_sub_rows, ]
    
    # classic OLM invervsion lets go
    # β = (XTX)−1XTy
    beta_vector <- MASS::ginv(t(X_matrix_sub) %*% X_matrix_sub) %*% 
      t(X_matrix_sub) %*% Y_sub
    
    print(beta_vector)
    
    # now predict on the full set
    pred_sub <- X_matrix %*% beta_vector
    
    # and, you guessed it get_rmse
    z2 <- get_rmse(Y_matrix, pred_sub)
    
  } else {
    
    z2 <- 0
    
  }
    
  # ---------------
  # -- Part 3 -----
  # ---------------
  
  print(z1)
  print(z2)
  
  print(z1 * lambda_vec[1] + z2 * lambda_vec[2])
  z1 * lambda_vec[1] + z2 * lambda_vec[2]
}

# get the first set
S[c(1:k)] <- 1
S

# get the score
check_r2 <- get_score(S)
check_r2

# compare
# S, df_wide, magic_n, nsites, ndays, SCORE
# get_score(S, df_wide, magic_n, nsites, ndays, score_cols, 
# X_matrix, n_predictors, ID_vector, Y_matrix, lambda,  & 
#   &                   SCORE_z1, SCORE_z2
check_f2 <- .Fortran('get_score', 
                     S = as.integer(S),
                     df_wide = df_wide,
                     magic_n = as.integer(k),
                     nsites = as.integer(N),
                     ndays = as.integer(N_daymet),
                     score_cols = as.integer(3),
                     X_matrix = X_matrix,
                     n_predictors = as.integer(N_predictors),
                     ID_vector = as.integer(ID_vector),
                     Y_matrix = Y_matrix,
                     labmda = lambda_vec,
                     SCORE_z1 = 0.,
                     SCORE_z2 = 0.,
                     SCORE = 0.)
check_f2$SCORE

stopifnot(mean(abs(check_r2 - check_f2$SCORE)) < 0.001)

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////

# another check, the mse should be 0 if the k = nsites
# but you should never be able to run this in sim-ann because 
# this breaks the algorithm

# get the first set
S[c(1:N)] <- 1

# get the score
check_r2 <- get_score(S)
check_r2

check_zero <- .Fortran('get_score', 
                     S = as.integer(S),
                     df_wide = df_wide,
                     magic_n = as.integer(N),
                     nsites = as.integer(N),
                     ndays = as.integer(N_daymet),
                     score_cols = as.integer(3),
                     SCORE = 0.)
check_zero

stopifnot(mean(abs(check_r2 - check_zero$SCORE)) < 0.001)






