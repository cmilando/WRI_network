# --- fit 2D thin-plate spline smooths (low-rank TPS) ----------------------------
# scale coords for numerics
rngx <- range(df$x_coord); rngy <- range(df$y_coord)
df <- df %>%
  mutate(xs = (x_coord - rngx[1]) / diff(rngx),
         ys = (y_coord - rngy[1]) / diff(rngy))

m_green  <- gam(green  ~ s(xs, ys, bs = "tp", k = 60), data = df, method = "REML")
m_albedo <- gam(albedo ~ s(xs, ys, bs = "tp", k = 60), data = df, method = "REML")

# if "k' index" says k too low, increase k and refit
gam.check(m_green)   
gam.check(m_albedo)
# need to understand this a little more, but seems ok

# --- predict to a regular grid ---------------------------------------------------
x_pts <- seq(0, 100, by = 1)
y_pts <- seq(0, 100, by = 1)
grid  <- expand.grid(x_coord = x_pts, y_coord = y_pts) %>%
  mutate(xs = (x_coord - rngx[1]) / diff(rngx),
         ys = (y_coord - rngy[1]) / diff(rngy))

grid$green_hat  <- predict(m_green,  newdata = grid)
grid$albedo_hat <- predict(m_albedo, newdata = grid)

# plots
p4 <- ggplot(grid, aes(x_coord, y_coord)) +
  geom_point( shape = 3, size = 0.02) + ggtitle("a. Points")

p5 <- ggplot(grid, aes(x_coord, y_coord, fill = green_hat)) +
  geom_raster() +  
  geom_point(data = df, aes( x= x_coord, y =y_coord, fill = green), shape = 21) + 
  scale_fill_viridis_c() + 
  labs(title = "b. Green", fill = "pred")+ 
  theme(legend.position = 'bottom')

p6 <- ggplot(grid, aes(x_coord, y_coord, fill = albedo_hat)) +
  scale_fill_viridis_c() + 
  geom_raster() + 
  geom_point(data = df, aes( x= x_coord, y =y_coord, fill = albedo), shape = 21) + 
  labs(title = "c. Albedo", fill = "pred")+ 
  theme(legend.position = 'bottom')

p4 + p5 + p6