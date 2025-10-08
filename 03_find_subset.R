#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Test with one K-monitors
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
source("00_create_simulated_dataset.R")
source("01_model_airtemp.R")
source("02_annealing_prep.R")

system("R CMD SHLIB simann.f90")

dyn.unload("simann.so")

dyn.load("simann.so")

find_optimal_set <- function(k_in, rep_in, verbose = 0) {
  
  write.table(NA, file = paste0("tmp/k",k_in, "_rep", rep_in))
  
  dyn.load("simann.so")
  
  # initialize S to be relative to this k_in
  S <- rep(0, N)  
  S[c(1:k_in)] <- 1
  
  oo <- .Fortran("simann",
                 S = as.integer(S),
                 df_wide = df_wide,
                 magic_n = as.integer(k_in),
                 nsites = as.integer(N),
                 ndays = as.integer(N_daymet),
                 score_cols = as.integer(3),
                 SCORE = 0.,
                 cooling_rate = 0.95,
                 verbose = as.integer(verbose)) 
  
  return(list(SCORE = oo$SCORE, S = oo$S, k = k_in, rep = rep_in))
}

oo <- find_optimal_set(k_in = 50, rep_in = 1, verbose = 1)

# confirming math
get_score(oo$S)
oo$SCORE

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Now in parallel, computer for 1:N monitors, and maybe like 10 times per K
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

library(future)
library(future.apply)
plan(multisession)

# so what I'm learning here is that it probably makes sense
# to think of reasonable bounds for k before you start
# specifically, you can't do k = N, k must always be < N

test_grid <- expand_grid(k = seq(10, 100, by = 10), rep = 1:5)
test_grid
system("rm -r tmp/*")

# run in parallel, takes  
total_oo <- future_lapply(1:nrow(test_grid), \(i) {
  
  # you can't do k = N, k must always be < N
  if (test_grid$k[i] == N) return(NULL)
  
  find_optimal_set(k_in = test_grid$k[i], 
                   rep_in = test_grid$rep[i],
                   verbose = 0)
  
})

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Make the curve 
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

plot_df <- do.call(rbind, lapply(total_oo, \(l) data.frame(l$SCORE, l$k, l$rep)))

ggplot(plot_df) +
  geom_boxplot(aes(x = l.k, y = l.SCORE / l.k, group = l.k))

ggsave("img/fig2.png", width = 8, height = 4, dpi = 500)
