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
  
  # k_in = 5
  # verbose = 1
  
  write.table(NA, file = paste0("tmp/k",k_in, "_rep", rep_in))
  
  dyn.load("simann.so")
  
  # initialize S to be relative to this k_in
  S <- rep(0, N)  
  
  # randomly get k starting positions
  init <- sample(1:N, size = k_in, replace = F)
  
  S[init] <- 1

  oo <- .Fortran("simann",
                 S = as.integer(S),
                 df_wide = df_wide,
                 magic_n = as.integer(k_in),
                 nsites = as.integer(N),
                 ndays = as.integer(N_daymet),
                 q = assess_quantiles,
                 score_cols = as.integer(length(assess_quantiles)),
                 X_matrix = X_matrix,
                 n_predictors = as.integer(N_predictors),
                 ID_vector = as.integer(ID_vector),
                 Y_matrix = Y_matrix,
                 labmda_vec = lambda_vec,
                 SCORE_z1 = 0.,
                 SCORE_z2 = 0.,
                 SCORE = 0.,
                 cooling_rate = 0.95,
                 verbose = as.integer(verbose)) 
  
  return(list(SCORE = oo$SCORE, SCORE_z1 = oo$SCORE_z1, SCORE_z2 = oo$SCORE_z2,
              S = oo$S, k = k_in, rep = rep_in))
}

oo <- find_optimal_set(k_in = 50, rep_in = 1, verbose = 1)
oo

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

test_grid <- expand_grid(k = seq(70, 90, by = 4), rep = 1:5)
test_grid
system("mkdir tmp")
system("rm -r tmp/*")


# run in parallel, takes  
total_oo <- future_lapply(1:nrow(test_grid), \(i) {
  
  # you can't do k = N, k must always be < N
  # just for sim-ann, because the algorithm
  # you can get a score, which you checked earlier
  if (test_grid$k[i] == N) return(NULL)
  
  find_optimal_set(k_in = test_grid$k[i], 
                   rep_in = test_grid$rep[i],
                   verbose = 0)
  
}, future.seed = T)

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Make the curve 
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

plot_df <- do.call(rbind, lapply(total_oo, \(l) 
                                 data.frame(l$SCORE, l$SCORE_z1, l$SCORE_z2,
                                            l$k, l$rep)))

plot_df_long <- plot_df %>% 
  pivot_longer(cols = c(l.SCORE, l.SCORE_z1, l.SCORE_z2))

plot_med <- plot_df_long %>%
  group_by(l.k, name) %>% summarize(.groups = 'keep', med = median(value))


p1 <- ggplot(plot_df_long %>% filter(name == 'l.SCORE')) +
  geom_line(data = plot_med %>% filter(name == 'l.SCORE'),
            aes(x = l.k / N * 100, y = med / l.k),
            color = 'blue', linewidth = 1, linetype = '11') + 
  geom_boxplot(aes(x = l.k / N * 100, y = value / l.k, group = l.k),
               outlier.shape = NA) +
  xlab("% of monitors included") + 
  ylab("Composite error score per monitor")

p2 <- ggplot(plot_df_long %>% filter(name != 'l.SCORE')) +
  geom_line(data = plot_med %>% filter(name != 'l.SCORE'),
            aes(x = l.k / N * 100, y = med / l.k, color = name,
                group = name),
            linewidth = 1, linetype = '11') + 
  geom_boxplot(aes(x = l.k / N * 100, y = value / l.k, group = paste0(name,l.k),
                   color = name),
               outlier.shape = NA, position = 'identity') +
  scale_color_discrete(labels = c('feature coverage', 'response fit'),
                       name = 'Error type') +
  xlab("% of monitors included") + 
  ylab("Composite error score per monitor")  

p1 + ggtitle('a. Combined score') + p2 + ggtitle('b. Individual components')

ggsave("img/fig2_122341233.png", width = 8.9, height = 4.2, dpi = 500)
