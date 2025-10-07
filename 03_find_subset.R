source("03_prep_simann.R")

system("R CMD SHLIB simann.f90")
dyn.unload("simann.so")
dyn.load("simann.so")
oo <- .Fortran("simann",
               S = as.integer(S),
               magic_n = as.integer(length(which(S == 1))),
               np = as.integer(nrow(site_pairs_sub)),
               nsites = as.integer(length(unique(site_pairs_sub$site_id))),
               site_pairs_sub = aa,
               row_lookup = bb,
               penalty = 1,
               SCORE = 0.,
               cooling_rate = 0.95,
               verbose = as.integer(1)) # the closer to 1 the slower the cooling. 0.92 is fast

# confirming math
get_score(oo$S)
oo$SCORE

pop_pts <- site_pairs_sub %>% filter(site_id %in% which(oo$S == 1)) %>% pull(pop_id)

final_buffers <- st_buffer(filtered_points[which(oo$S == 1), ], RADIUS <- 20 * 1609.34)
# plot(final_buffers)

ggplot() +
  geom_sf(data = state_map_sf, fill = "white", color = "black") +
  geom_sf(data = population_points_sf, color = "forestgreen", size = 0.5) +
  geom_sf(data = highway_lines_sf, color = "grey", size = 0.1, alpha = 0.5) +
  geom_sf(data = filtered_points[which(oo$S == 1), ], color = "red", size = 2) +
  geom_sf(data = final_buffers, color = "red", size = 5, fill = NA) +
  geom_sf(data = population_points_sf[pop_pts, ], color = "orange", size = 1) +
  #geom_sf(data = buffers[673,], color = 'red', fill = NA) +
  theme_minimal() +
  labs(title = "Population Points and Major Highways in Selected Northeastern States")
