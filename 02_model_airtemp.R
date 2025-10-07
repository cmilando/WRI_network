
# create a linear model to recover air_temp prediction
lm1 <- glm(air_temp ~ green + albedo, data = df)

summary(lm1)
# we'll 


