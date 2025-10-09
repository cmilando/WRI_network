#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' How to run with a new dataset
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

library(tidyverse)
library(patchwork)

library(future)
library(future.apply)
plan(multisession)

system("R CMD SHLIB simann.f90")
dyn.load("simann.so")

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Dataset
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

# Load the dataset
df <- readRDS("test_df.RDS")

# df needs coordinates cols
coords_cols   <- c('x_coord', 'y_coord')

# names of the predictor columns
predictor_cols <- c('daymet', 'green', 'albedo')

# df needs a integer monitor_id
monitor_id_col <- c('monitor_id')

# df needs a response variable
response_col   <- c('air_temp')

# df needs a day_id
day_id_col <- c('day_id')

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Optimization params
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

# quantiles to assess feature similarity
# moved in from the edges to be slightly more robust

# assess_quantiles <- c(0.05, 0.50, 0.95)
assess_quantiles <- c(0.50)

# this assigns the weights to the different components of the composite score
# position 1 is for predictors similarity
# position 2 is for the output of the linear model
lambda_vec <- c(1, 1)

# this roughly determines how quickly the simulated annealing converges
# 0.95 is a good starting point
cooling_rate <- 0.95

# number of monitoring groups to assess at
# number of repetitions at each group
# note: you many want to change the test_grid below if these values
#       are really large, I tested this with a dataset of 100 Monitors
#       but if you have 1000 that will be probably not great
N_groups     <- 10
N_group_reps <- 5

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Validation of df columns - Fortran really cares about types
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

# names of the coordinates
stopifnot(all(coords_cols %in% colnames(df)))
stopifnot(all(typeof(unlist(df[,coords_cols])) == 'double'))

# names of the linear model predictors
stopifnot(all(predictor_cols %in% colnames(df)))
stopifnot(all(typeof(unlist(df[,predictor_cols])) == 'double'))

# now set N_predictors
N_predictors <- length(predictor_cols)

# name of the monitor ID column
stopifnot(all(monitor_id_col %in% colnames(df)))
stopifnot(typeof(unlist(df[,monitor_id_col])) == 'integer')

# from this set N
N <- length(unlist(unique(df[, monitor_id_col])))
stopifnot(all(unlist(unique(df[, monitor_id_col])) %in% 1:N))

# name of the response column
stopifnot(all(response_col %in% colnames(df)))
stopifnot(typeof(unlist(df[,response_col])) == 'double')

# name of the day_id
stopifnot(all(day_id_col %in% colnames(df)))
stopifnot(typeof(unlist(df[,day_id_col])) == 'integer')

# from this set N days
n_days <- length(unlist(unique(df[, day_id_col])))
stopifnot(all(unlist(unique(df[, day_id_col])) %in% 1:n_days))

# check that df is N x Ndays
stopifnot(nrow(df) == N * n_days)

# check lambda
stopifnot(all(lambda_vec >= 0))

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Transformed data  
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

# probably the way to do this in FORTRAN is to convert df to wide first
# because then its just rows and columns
# and probably becomes easier to subset
df_wide <- df %>% pivot_wider(id_cols = !!monitor_id_col, 
                              names_from = !!day_id_col,
                              values_from = !!response_col) %>% 
  select(-!!monitor_id_col) %>% as.matrix()

stopifnot(all(dim(df_wide) == c(N, n_days)))

# also have to do some treatment to the co-variates
Y_matrix    <- as.matrix(df[, response_col])
X_matrix    <- as.matrix(df[, predictor_cols])
ID_vector   <- as.integer(unlist(df[, monitor_id_col]))

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Function for simulated annealing with one test run  
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

find_optimal_set <- function(k_in, rep_in, verbose = 0) {
  
  # ***********
  # k_in = 5
  # verbose = 1
  # ***********
  
  write.table(NA, file = paste0("tmp/k", k_in, "_rep", rep_in))
  
  dyn.load("simann.so")
  
  # initialize S to be relative to this k_in
  S <- rep(0, N)  
  
  # randomly get k starting positions
  init <- sample(1:N, size = k_in, replace = F)
  
  S[init] <- 1
  
  oo <- .Fortran("simann",
                 S            = as.integer(S),
                 df_wide      = df_wide,
                 magic_n      = as.integer(k_in),
                 nsites       = as.integer(N),
                 ndays        = as.integer(n_days),
                 q            = assess_quantiles,
                 score_cols   = as.integer(length(assess_quantiles)),
                 X_matrix     = X_matrix,
                 n_predictors = as.integer(N_predictors),
                 ID_vector    = as.integer(ID_vector),
                 Y_matrix     = Y_matrix,
                 labmda_vec   = lambda_vec,
                 SCORE_z1     = 0.,
                 SCORE_z2     = 0.,
                 SCORE        = 0.,
                 cooling_rate = cooling_rate,
                 verbose      = as.integer(verbose)) 
  
  return(list(SCORE = oo$SCORE, SCORE_z1 = oo$SCORE_z1, SCORE_z2 = oo$SCORE_z2,
              S = oo$S, k = k_in, rep = rep_in))
}

oo <- find_optimal_set(k_in = round(0.5 * N), rep_in = 1, verbose = 1)

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Now in parallel, computer for 1:N monitors, and maybe like 10 times per K
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

# so what I'm learning here is that it probably makes sense
# to think of reasonable bounds for k before you start
# specifically, you can't do k = N, k must always be < N

k_vec <- round(seq(0, N, length.out = (N_groups+1))[2:(N_groups+1)])

test_grid <- expand_grid(k = k_vec, rep = 1:N_group_reps)
test_grid

system("mkdir tmp")
system("rm -r tmp/*")

# run in parallel, takes a few minutes total
# can watch progress in /tmp
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

ggsave("img/curves.png", width = 8.9, height = 4.2, dpi = 600)

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' Make the spatial plot 
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================

S_plot_df <- do.call(rbind, lapply(total_oo, \(l) {
  if(any(is.null(l))) return(NULL)
  x1 <- data.frame(S = l$S)
  x1$monitor_id = 1:nrow(x1)
  colnames(x1)[2] = monitor_id_col
  x1 <- left_join(x1, unique(df[, c(monitor_id_col, coords_cols)]))
  x1$k <- l$k
  x1$rep <- l$rep
  x1
}))

head(S_plot_df)

##
S_plot_df_agg <- S_plot_df %>% 
  group_by(x_coord, y_coord, monitor_id, k) %>%
  summarize(
    .groups = 'keep',
    S_sum = sum(S)
  ) %>%
  mutate(S_prob = S_sum / N_group_reps)

pt_size = 3
pt_shape = 21

p1 <- ggplot(S_plot_df_agg, 
       aes(x_coord, y_coord, fill = S_prob))  + 
  geom_point(size = pt_size, shape = pt_shape, show.legend = F) + 
  scale_fill_viridis_c() + 
  ggtitle("Probability of selection - by k") +
  facet_wrap(~k)

##
S_plot_df_agg2 <- S_plot_df %>% 
  group_by(x_coord, y_coord, monitor_id) %>%
  summarize(
    .groups = 'keep',
    n = n(),
    S_sum = sum(S)
  ) %>%
  mutate(S_prob = S_sum / n)

S_plot_df_agg2

p2 <- ggplot(S_plot_df_agg2, 
       aes(x_coord, y_coord, fill = S_prob))  + 
  geom_point(size = pt_size, shape = pt_shape, show.legend = F) + 
  scale_fill_viridis_c() + coord_equal() +
  ggtitle("Probability of selection - overall")

p1 + p2 + plot_layout(widths = c(0.8, 0.4))

ggsave("img/P_by_k.png", width = 10.29167, height = 6.37500, dpi = 600)
