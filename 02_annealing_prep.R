library(tidyverse)
library(sf)
library(pbapply)
load("all_objects.RData")

# ok now set up for simann

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  #geom_sf(data = highway_lines_sf, color = "orange", size = 0.1, alpha = 0.5) +
  geom_sf(data = filtered_points, color = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")

# the SimAnn question -- what selection of 5 red points 
# with a set radius of Y
# gives you the most green points

# first pass - you could do great circle distance, but 
# for a first pass just assume flat

RADIUS <- 20 * 1609.34
potential_sites <- st_coordinates(filtered_points)
population_dots <- st_coordinates(st_transform(population_points_sf, crs = 3857))

head(potential_sites)
head(population_dots)

# distance formula
# ^ you only need to calculate this once
# and can be done outside of the loop

# and you could also do this in FORTRAN again, could do great-circle distance
# if you really wanted to

site_pairs <- expand_grid(site_id = 1:nrow(potential_sites),
                          pop_id = 1:nrow(population_dots))

site_dist <- pbsapply(1:nrow(site_pairs), function(i) {
  site_i = site_pairs$site_id[i]
  pop_i = site_pairs$pop_id[i]
  sqrt((potential_sites[site_i,1] - population_dots[pop_i, 1])^2 +
         (potential_sites[site_i,2] - population_dots[pop_i, 2])^2)
})

summary(site_dist)
which.min(site_dist)
site_pairs[692220, ]
length(population_points_sf)

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  geom_sf(data = population_points_sf[760, ], color = "forestgreen", size = 0.5) +
  #geom_sf(data = highway_lines_sf, color = "orange", size = 0.1, alpha = 0.5) +
  geom_sf(data = filtered_points[771, ], color = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")

# so first you can subset to only ones within the radius
site_pairs$dist <- site_dist
dim(site_pairs)

site_pairs_sub <- site_pairs %>% filter(site_dist < RADIUS)
dim(site_pairs_sub)

# yeah this is now only 93,920 sites-pop sites
length(unique(site_pairs_sub$site_id))
nrow(potential_sites)

# so i think the way to do this is to get clever with your indices

# all these variables should be range01 for additivity
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

set.seed(123)
site_pairs_sub$population <- 0.5

metric_df <- data.frame(pop_id = unique(site_pairs_sub$pop_id))

metric_df$metric1 <- range01(rnorm(1:nrow(metric_df)))
site_pairs_sub <- left_join(site_pairs_sub, metric_df)
dim(site_pairs_sub)

hist(site_pairs_sub$metric1)

head(site_pairs_sub)
tail(site_pairs_sub)
summary(site_pairs_sub)

beta_weights <- c(10, 4)

score = as.matrix(site_pairs_sub[, c('population', 'metric1')]) %*% beta_weights

site_pairs_sub$score = as.vector(score)

head(site_pairs_sub)

# make sure its arranged
stopifnot(site_pairs_sub$site_id == sort(site_pairs_sub$site_id))

rlepX <- rle(site_pairs_sub$site_id)

# Compute the indices of the start and end of each run
end <- cumsum(rlepX$lengths)
start <- c(1, lag(end)[-1] + 1)

xx <- data.frame(site_id = unique(site_pairs_sub$site_id), start, end)
head(xx)

#### WHAT ARE YOU PASSING IN
head(xx)

head(site_pairs_sub)

S <- rep(0, length(unique(site_pairs_sub$site_id)))
S[c(1:5)] <- 1

penalty <- 1


# the tricky part is the cancelling out.
# you know the point to site matrix though so ..

# (1) choose a set of 5
set.seed(123)


#S_ones = sample(xx$site_id, size = length(), replace = F)

# (2) update the indicator column
# this has to be a group multiplier, not sequential

#
get_score <- function(S_local) {

  S_ones <- which(S_local == 1)

  # get just the points for this sample
  yy <- do.call(rbind,lapply(1:length(S_ones), function(i) {
    this_block <- xx[S_ones[i],]
    site_pairs_sub[this_block$start:this_block$end, ]
  }))

  ## now make a frequency table of pop_id
  ## you could probably do this differently
  ## in a way that only requires going through the list once
  freq <- data.frame(table(yy$pop_id), stringsAsFactors = F)
  names(freq) <- c("pop_id", 'freq')
  freq$pop_id <- as.integer(as.character(freq$pop_id))

  yy <- left_join(yy, freq, by = join_by(pop_id))

  yy$score2 <- yy$score * 1 / (penalty * yy$freq)

  ## get sum
  -1*sum(yy$score2)

}

get_score(S)

aa <- as.matrix(site_pairs_sub, 
                nrow = nrow(site_pairs_sub),
                ncol = 6)
head(aa)
dim(aa)

bb = as.matrix(xx, 
               nrow = nrow(xx),
               ncol = ncol(xx))
head(bb)
dim(bb)



