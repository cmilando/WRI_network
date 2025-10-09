#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' This script creates the simulated data
#' with some minor covariance
#' 
#' change the beta coefficients to obeserve a different effect later
#' a little help from ChatGPT to get the spatial covariance part done
#' 
#' ////////////////////////////////////////////////////////////////////////////
#' ============================================================================
library(tidyverse)
library(mgcv)
library(patchwork)

# set.seed(123)
set.seed(122341233)

N <- 50             # number of monitors

beta_green  <- 2     # beta coefficient for greenness
beta_albedo <- 1.5   # beta coefficient for albedo

N_daymet    <- 100    # number of days of daily meteorology
beta_daymet <- 10    # beta coefficient for the meteorology parameter

air_temp_sd <- 2     # st dev for random error, increase for more noise

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
#' using the above variables, get the simulated data, 
#' details here don't really matter
source("functions.R")
df <- get_simulated_df()
head(df)

saveRDS(df, "test_df.RDS")

#' ============================================================================
#' ////////////////////////////////////////////////////////////////////////////
# plot
pt_size = 3
pt_shape = 21

p1 <- ggplot(df %>% filter(day_id == 1), 
             aes(x_coord, y_coord, fill = green))  + 
  geom_point(size = pt_size, shape = pt_shape, show.legend = F) + 
  scale_fill_viridis_c() + ggtitle("a. Measured Green") + 
  theme(legend.position = 'bottom')

p2 <- ggplot(df %>% filter(day_id == 1), 
             aes(x_coord, y_coord, fill = albedo)) + 
  geom_point(size = pt_size, shape = pt_shape, show.legend = F) + 
  scale_fill_viridis_c() + ggtitle("b. Measured Albedo") + 
  theme(legend.position = 'bottom')

p3 <- ggplot(df %>% filter(day_id == 1), 
             aes(x_coord, y_coord, fill = daymet)) + 
  geom_point(size = pt_size, shape = pt_shape, show.legend = F) + 
  scale_fill_viridis_c() + ggtitle("c. DayMet (day = 1)") + 
  theme(legend.position = 'bottom')

p4 <- ggplot(df %>% filter(day_id == 1), 
             aes(x_coord, y_coord, fill = air_temp)) + 
  geom_point(size = pt_size, shape = pt_shape, show.legend = F) + 
  scale_fill_viridis_c() + ggtitle("d. Measured Air temp (day = 1)") + 
  theme(legend.position = 'bottom')

(p1 + p2) / (p3 + p4) 
ggsave("img/features.png", width = 8, height = 6, dpi = 300)


